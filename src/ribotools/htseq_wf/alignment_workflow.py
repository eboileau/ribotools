#! /usr/bin/env python3

"""Provide wrapper for general alignment workflow.

(1) Trim adapters using Flexbar.
(2) Remove rRNA, tRNA, etc. using Bowtie 2.
(3) Align reads to the genome (default) using STAR.

Note* Bowtie 2 and STAR indices must already be available.
      For mapping with STAR, annotations are used on the fly.

      This script is a generalisation of create-base-genome-profile
      from the rpbp package.
"""

import os
import sys
import argparse
import logging
import shlex
import yaml

import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.utils as utils
import pbio.misc.slurm as slurm

import pbio.utils.bam_utils as bam_utils
import pbio.utils.fastx_utils as fastx_utils
import pbio.utils.pgrm_utils as pgrm_utils

import pbio.ribo.ribo_utils as ribo_utils
import pbio.ribo.ribo_filenames as filenames

import pproc.utils.cl_utils as clu

from pproc.defaults import default_num_cpus, default_mem, default_gff_spec,\
    star_executable, star_options, flexbar_options, metagene_options

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper script for a general RNA- or 
        Ribo-Seq workflow: flexbar -> bowtie2 -> STAR. File names and directory structures 
        follow the conventions used in pbio.ribo.riboseq-utils to be compatible with Rp-Bp.""")

    parser.add_argument('seq', choices=['rna', 'ribo'])

    parser.add_argument('raw_data', help="The raw fastq[.gz] file.")

    parser.add_argument('config', help="The yaml config file.")

    parser.add_argument('name', help="The name of the dataset.")

    parser.add_argument('--trim-rna-to-max-fragment-size', help="""Flag: trim RNA post 
        adapter removal using max fragment size from the matching Ribo-seq sample. Note* At least
        the "periodic-offsets" file must be available. The config file must also include 
        "matching_samples" and the path to the Ribo-seq config must be given [--ribo-config]).
        If the [--post-trim-length] option was passed via [flexbar-options], it will 
        silently override this option.""", action='store_true')

    parser.add_argument('--ribo-config', help="""Optional argument: the Ribo-seq config file
        if using [--trim-rna-to-max-fragment-size].""",
                        required='--trim-rna-to-max-fragment-size' in sys.argv, type=str)

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    pgrm_utils.add_flexbar_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[alignment-workflow]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    # if using slurm, submit the script
    if args.use_slurm:
        cmd = "{}".format(' '.join(shlex.quote(s) for s in sys.argv))
        slurm.check_sbatch(cmd, args=args)
        return

    # check that all of the necessary programs are callable
    programs = [
        'flexbar',
        'bowtie2',
        args.star_executable,
        'samtools']
    shell_utils.check_programs_exist(programs)

    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    note = config.get('note', None)

    base_keys = ['ribosomal_index',
                 'star_index',
                 'genome_base_path',
                 'genome_name',
                 'fasta',
                 'gtf']

    # only for RNA if using [trim-rna-to-max-fragment-size]
    filename_length = None

    if args.seq == 'rna':
        config_keys = base_keys + ['rnaseq_data']
        utils.check_keys_exist(config, config_keys)
        seq_data = config['rnaseq_data']
        adapter_seq_str = utils.get_config_argument(config, 'rna_adapter_sequence', 'adapter-seq')
        adapter_file_str = utils.get_config_argument(config, 'rna_adapter_file', 'adapters')
        star_output_prefix = filenames.get_rnaseq_bam_base(seq_data,
                                                           args.name,
                                                           note=note)

        if args.trim_rna_to_max_fragment_size:

            ribo_config = yaml.load(open(args.ribo_config), Loader=yaml.FullLoader)
            is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)

            matching_ribo_sample = config['matching_samples'][args.name]

            # get the lengths, we don't need the offsets
            lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
                ribo_config,
                matching_ribo_sample,
                is_unique=is_unique_ribo,
                default_params=metagene_options)

            if len(lengths) == 0:
                msg = """No periodic read lengths were found, but the
                                [trim-rna-to-max-fragment-size] option was given!"""
                logger.critical(msg)
                return

            max_length = max([int(l) for l in lengths])
            filename_length = max_length

            # ask flexbar to trim to specified read length from 3' end after removal
            # if called from 'run-htseq-worflow', we are fine, but if called independently
            # and ['post-trim-length'] is passed via [flexbar-options], then it will
            # override this value silently
            flexbar_options['post-trim-length'] = max_length

    else:
        config_keys = base_keys + ['riboseq_data']
        utils.check_keys_exist(config, config_keys)
        seq_data = config['riboseq_data']
        adapter_seq_str = utils.get_config_argument(config, 'adapter_sequence', 'adapter-seq')
        adapter_file_str = utils.get_config_argument(config, 'adapter_file', 'adapters')
        star_output_prefix = filenames.get_riboseq_bam_base(seq_data,
                                                            args.name,
                                                            note=note)

    # (1) Run flexbar to remove adapter sequences.
    # -------------------------------------------
    raw_data = args.raw_data

    flexbar_target = filenames.get_without_adapters_base(seq_data,
                                                         args.name,
                                                         note=note)

    without_adapters = filenames.get_without_adapters_fastq(seq_data,
                                                            args.name,
                                                            note=note)

    # get all options, command line options override defaults
    flexbar_option_str = pgrm_utils.get_final_args(flexbar_options, args.flexbar_options)

    cmd = "flexbar -r {} -t {} {} {} {} -n {}".format(raw_data,
                                                      flexbar_target,
                                                      adapter_seq_str,
                                                      adapter_file_str,
                                                      flexbar_option_str,
                                                      args.num_cpus)
    in_files = [raw_data]
    out_files = [without_adapters]
    file_checkers = {
        without_adapters: fastx_utils.check_fastq_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite, call=call)

    # (2) Run Bowtie2 to remove rRNA, tRNA, etc. alignments.
    # ------------------------------------------------------
    out = utils.abspath("dev", "null")  # these are discarded anyway

    without_rrna = filenames.get_without_rrna_fastq(seq_data,
                                                    args.name,
                                                    note=note)

    with_rrna = filenames.get_with_rrna_fastq(seq_data,
                                              args.name,
                                              note=note)

    cmd = "bowtie2 -p {} --very-fast -x {} -U {} -S {} --un-gz {} --al-gz {}".format(
        args.num_cpus,
        config['ribosomal_index'],
        without_adapters,
        out,
        without_rrna,
        with_rrna)

    in_files = [without_adapters]
    in_files.extend(pgrm_utils.get_bowtie2_index_files(config['ribosomal_index']))
    out_files = [without_rrna, with_rrna]
    to_delete = [without_adapters, with_rrna]
    file_checkers = {
        without_rrna: fastx_utils.check_fastq_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite, call=call,
                                   keep_delete_files=keep_delete_files, to_delete=to_delete)

    # (3) Run STAR to align rRNA-depleted reads to the genome.
    # -------------------------------------------------------
    star_genome_bam = "{}{}".format(star_output_prefix, "Aligned.sortedByCoord.out.bam")

    # get all options, command line options override defaults

    mem_bytes = utils.human2bytes(args.mem)
    star_options['limitBAMsortRAM'] = mem_bytes

    if args.tmp is not None:
        star_tmp_name = str(args.name + "_STARtmp")
        star_tmp_dir = pgrm_utils.create_star_tmp(args.tmp, star_tmp_name)
        star_options['outTmpDir'] = star_tmp_dir

    star_option_str = pgrm_utils.get_final_args(star_options, args.star_options)

    # If gff3 spec, then we need to inform STAR.
    sjdb_gtf_tag_str = ""
    gff_spec = config.get('gff_spec', default_gff_spec)
    if gff_spec.lower() == 'gff3':
        sjdb_gtf_tag_str = "--sjdbGTFtagExonParentTranscript Parent"

    cmd = ("{} --runThreadN {} --genomeDir {} --sjdbGTFfile {} {} --readFilesIn {} "
           "{} --outFileNamePrefix {}".format(args.star_executable,
                                              args.num_cpus,
                                              config['star_index'],
                                              config['gtf'],
                                              sjdb_gtf_tag_str,
                                              without_rrna,
                                              star_option_str,
                                              star_output_prefix))
    in_files = [without_rrna]
    in_files.extend(pgrm_utils.get_star_index_files(config['star_index']))
    out_files = [star_genome_bam]
    file_checkers = {star_genome_bam: bam_utils.check_bam_file}
    to_delete = [without_rrna]

    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite, call=call,
                                   keep_delete_files=keep_delete_files, to_delete=to_delete)

    # symlink the (genome) STAR output
    indexed_genome_bam = filenames.get_seq_bam(args.seq,
                                               seq_data,
                                               args.name,
                                               length=filename_length,
                                               note=note)
    if os.path.exists(star_genome_bam):
        shell_utils.create_symlink(star_genome_bam, indexed_genome_bam, remove=False, call=call)
    else:
        msg = """Could not find the STAR genome bam alignment file. Unless the
              [--do-not-call] option was given, this is a problem!"""
        logger.warning(msg)

    # create the bamtools index if it does not already exists
    bam_utils.index_bam_file(indexed_genome_bam, args)

    # check if we want to keep multimappers
    if ((args.seq == "rna" and "keep_rnaseq_multimappers" in config) or
            (args.seq == "ribo" and "keep_riboseq_multimappers" in config)):
        return

    # remove multi mappers from the genome alignments
    is_unique = True
    keep_genome_bam = filenames.get_seq_bam(args.seq,
                                            seq_data,
                                            args.name,
                                            is_unique=is_unique,
                                            length=filename_length,
                                            note=note)
    in_files = [indexed_genome_bam]
    out_files = [keep_genome_bam]
    to_delete = [indexed_genome_bam]
    file_checkers = {
        keep_genome_bam: bam_utils.check_bam_file
    }
    utils.call_func_if_not_exists(bam_utils.remove_multimappers, out_files,
                                  in_files=in_files, overwrite=args.overwrite,
                                  call=call, file_checkers=file_checkers, to_delete=to_delete,
                                  keep_delete_files=keep_delete_files)


if __name__ == '__main__':
    main()
