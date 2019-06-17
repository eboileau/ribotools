#! /usr/bin/env python3

"""Call htseq-count on a list of samples

Functions:
    get_count_table
"""

import sys
import os
import argparse
import logging
import yaml
import shlex

import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.utils as utils
import pbio.misc.slurm as slurm

import pbio.utils.pgrm_utils as pgrm_utils

import pbio.ribo.ribo_filenames as filenames
import pbio.ribo.ribo_utils as ribo_utils

import pproc.utils.cl_utils as clu

from pproc.defaults import default_num_cpus, default_mem, htseq_options, \
    htseq_executable, metagene_options

logger = logging.getLogger(__name__)


def get_count_table(seq_base, name, length=None, is_unique=False, note=None):

    unique_str = filenames.get_unique_string(is_unique)
    length_str = filenames.get_length_string(length)
    note_str = filenames.get_note_string(note)

    fn = ''.join([name, note_str, unique_str, length_str, '.tsv'])

    return os.path.join(seq_base, 'count-tables', fn)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Call htseq-count and extract read counts.""")

    parser.add_argument('seq', choices=['rna', 'ribo'])

    parser.add_argument('config', help="The yaml config file.")

    parser.add_argument('name', help="The name of the dataset.")

    parser.add_argument('--not-periodic', help="""Flag: non-periodic read 
            lengths are NOT filtered out; this is to get the right file names 
            (length-). For Ribo-seq only.""", action='store_true')

    parser.add_argument('--ribo-config', help="""Optional argument: the Ribo-seq config file
            when seq is RNA and RNA reads have been trimmed to max fragment size from 
            the matching Ribo-seq sample. In addition, the RNA config file must include 
            "matching_samples". If not given, then the 'normal' RNA alignment files will
            be used.""", type=str, default=None)

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)
    clu.add_htseq_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[call-htseq-count]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    # if using slurm, submit the script
    if args.use_slurm:
        cmd = "{}".format(' '.join(shlex.quote(s) for s in sys.argv))
        slurm.check_sbatch(cmd, args=args)
        return

    # check that all of the necessary programs are callable
    programs = [htseq_executable]
    shell_utils.check_programs_exist(programs)

    if ((args.seq.lower() == 'ribo' and args.ribo_config) or
            (args.seq.lower() == 'rna' and args.not_periodic)):
        msg = """seq is Ribo and [--ribo-config] was passed or seq is 
            RNA and [--not-periodic] was passed. These options will
            be ignored."""
        logger.warning(msg)

    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    sample_name_map = ribo_utils.get_sample_name_map(config)
    note = config.get('note', None)

    keep_key = 'keep_' + str(args.seq) + 'seq_multimappers'
    data_key = str(args.seq) + 'seq_data'
    is_unique = not (keep_key in config)

    if args.seq.lower() == 'ribo' and not args.not_periodic:

        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            args.name,
            is_unique=is_unique,
            default_params=metagene_options
        )

        if len(lengths) == 0:
            msg = "No periodic read lengths and offsets were found!"
            logger.critical(msg)
            return

    elif args.seq.lower() == 'rna' and args.ribo_config:

        config_keys = ['matching_samples']
        utils.check_keys_exist(config, config_keys)
        matching_ribo_sample = config['matching_samples'][args.name]

        ribo_config = yaml.load(open(args.ribo_config), Loader=yaml.FullLoader)
        is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)

        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            ribo_config,
            matching_ribo_sample,
            is_unique=is_unique_ribo,
            default_params=metagene_options
        )

        if len(lengths) == 0:
            msg = """No periodic read lengths and offsets were found!, but the
                                [--ribo-config] option was given!"""
            logger.critical(msg)
            return

        lengths = str(max([int(l) for l in lengths]))

    else:
        lengths = None

    bam_file = filenames.get_seq_bam(
        args.seq,
        config[data_key],
        args.name,
        is_unique=is_unique,
        length=lengths,
        note=note
    )

    count_file = get_count_table(config[data_key],
                                 sample_name_map[args.name],
                                 is_unique=is_unique,
                                 length=lengths,
                                 note=note)

    # all options to htseq-count
    htseq_options_str = pgrm_utils.get_final_args(htseq_options, args.htseq_options)

    cmd = "{} {} {} {} > {}".format(
        htseq_executable,
        htseq_options_str,
        bam_file,
        config['gtf'],
        count_file)

    in_files = [bam_file, config['gtf']]
    out_files = [count_file]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call,
                                   keep_delete_files=keep_delete_files)


if __name__ == '__main__':
    main()
