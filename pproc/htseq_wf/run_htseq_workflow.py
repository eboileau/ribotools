#! /usr/bin/env python3

"""Helper script to submit a set of RNA- or Ribo-seq
samples for alignment (Flexbar, Bowtie 2, STAR), and count
reads mapping to selected features (HTSeq).

Note* Example workflow: submit Ribo-seq samples for
      mapping and periodicity estimation (or else use
      data available from the Rp-Bp pipeline). When
      completed, submit RNA-seq samples, trimming
      reads to max periodic fragment length of matching
      Ribo-seq sample.

      Ribo-seq: Estimate read length periodicity,
      and filter out non-periodic fragments from the
      final BAM files (default). If periodicity estimates
      and/or mapped reads were previously obtained by running
      rpbp, they must be available. Reads are NOT P-site
      offset-shifted.

      Ribo-seq: Genome (unique) alignments are used to
      construct the profiles, which may include anti-sense
      reads.

      RNA-seq: Reads can be trimmed to max periodic fragment
      length from the matching Ribo-seq data before mapping,
      using the same default mapping options. Ribo-seq periodicity
      estimates must be available. For general long read RNA-seq
      mapping, these default options must be overridden
      via command line using 'star-options' (and/or 'flexbar-options').

      HTSeq: Deals with library strandedness. All options are passed
      via this script.
"""

import os
import sys
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

from pproc.defaults import default_num_cpus, default_mem, star_executable

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Submit a set of samples to RNA- or 
        Ribo-seq alignment pipelines. The pipeline is called for every sample in the 
        configuration file.""")

    parser.add_argument('seq', choices=['rna', 'ribo'])

    parser.add_argument('config', help="The yaml configuration file.")

    parser.add_argument('-skip', '--skip-periodicity-estimation', help="""Flag: by default, 
        estimate periodicity and filter non-periodic read lengths from the final alignment 
        files. If this flag is passed, then only mapping is performed. For Ribo-seq only.""",
                        action='store_true')

    parser.add_argument('--run-all', help="""Flag: map and count Ribo-seq and RNA-seq, one
        after the other, provided that ALL Ribo-seq jobs completed successfully.""",
                        action='store_true')

    parser.add_argument('--trim-rna-to-max-fragment-size', help="""Flag: trim RNA post 
        adapter removal using max fragment size from the matching Ribo-seq sample. Note* At least
        the "periodic-offsets" file must be available. The config file must also include 
        "matching_samples" and the path to the Ribo-seq config must be given [--ribo-config])""",
        action='store_true')

    parser.add_argument('--ribo-config', help="""Optional argument: the Ribo-seq config file
        if using [--trim-rna-to-max-fragment-size].""",
        required='--trim-rna-to-max-fragment-size' in sys.argv, type=str)

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    pgrm_utils.add_flexbar_options(parser)
    clu.add_htseq_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # check that all of the necessary programs are callable
    programs = [
        'flexbar',
        'bowtie2',
        args.star_executable,
        'samtools',
        'htseq-workflow',
        'get-ribo-periodic',
        'keep-ribo-periodic'
    ]
    shell_utils.check_programs_exist(programs)

    # handle all option strings to call the pipeline script
    logging_str = logging_utils.get_logging_options_string(args)
    file_str = clu.get_file_options_string(args)
    slurm_str = slurm.get_slurm_options_string(args)
    star_str = pgrm_utils.get_star_options_string(args)
    flexbar_str = pgrm_utils.get_flexbar_options_string(args)
    htseq_str = clu.get_htseq_options_string(args)
    tmp_str = ""

    # handle do_not_call so that we do call the pipeline script, but that it does not run anything
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"
    args.do_not_call = False

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    base_keys = ['ribosomal_index',
                 'star_index',
                 'genome_base_path',
                 'genome_name',
                 'fasta',
                 'gtf']

    all_ribo_jobs = None

    # now process samples
    if args.seq == 'ribo':

        config_keys = base_keys + ['riboseq_data', 'riboseq_samples']
        utils.check_keys_exist(config, config_keys)

        # first submit all alignment jobs, then estimate periodicity
        # for each one
        job_ids_mapping = {}
        for sample_name, data in config['riboseq_samples'].items():

            if args.tmp is not None:
                tmp = os.path.join(args.tmp, "{}_htseq_wf_ribo".format(sample_name))
                tmp_str = "--tmp {}".format(shlex.quote(tmp))

            cmd = "alignment-workflow {} {} {} {} {} {} {} {} {} {} {}".format(
                args.seq,
                data,
                args.config,
                sample_name,
                do_not_call_str,
                logging_str,
                tmp_str,
                file_str,
                slurm_str,
                star_str,
                flexbar_str)

            job_id = slurm.check_sbatch(cmd, args=args)
            job_ids_mapping[sample_name] = job_id

        job_ids_periodic = []
        if not args.skip_periodicity_estimation:

            filter_non_periodic_str = '--filter-non-periodic'

            for sample_name in config['riboseq_samples'].keys():

                cmd = "get-ribo-periodic {} {} {} {} {} {} {}".format(
                    args.config,
                    sample_name,
                    filter_non_periodic_str,
                    do_not_call_str,
                    logging_str,
                    file_str,
                    slurm_str)

                job_id = [job_ids_mapping[sample_name]]
                job_id_periodic = slurm.check_sbatch(cmd, args=args, dependencies=job_id)
                job_ids_periodic.append(job_id_periodic)

        # collect all jobs and submit to htseq-count
        all_ribo_jobs = list(job_ids_mapping.values())
        all_ribo_jobs.extend(job_ids_periodic)

        cdm = "".format()

        slurm.check_sbatch(cmd, args=args, dependencies=all_ribo_jobs)

        # we are done with ribo, however if we run all, we need to continue
        if args.run_all:
            args.seq = 'rna'

    if args.seq == 'rna':

        config_keys = base_keys + ['rnaseq_data', 'rnaseq_samples']

        trim_str = ''
        ribo_cfg_str = ''

        if args.trim_rna_to_max_fragment_size:

            config_keys.extend(['matching_samples'])

            if '--post-trim-length' in flexbar_str:
                msg = ("Flexbar option [--post-trim-length] already given..."
                   "will ignore [trim-rna-to-max-fragment-size].")
                logger.critical(msg)
            else:
                trim_str = '--trim-rna-to-max-fragment-size'
                ribo_cfg_str = '--ribo-config {}'.format(shlex.quote(args.ribo_config))

        if not args.trim_rna_to_max_fragment_size and args.ribo_config:
            msg = """The [--ribo-config] option is passed without
            [trim-rna-to-max-fragment-size], it will be ignored."""
            logger.warning(msg)

        utils.check_keys_exist(config, config_keys)

        all_rna_jobs = []
        for sample_name, data in config['rnaseq_samples'].items():

            if args.tmp is not None:
                tmp = os.path.join(args.tmp, "{}_htseq_wf_rna".format(sample_name))
                tmp_str = "--tmp {}".format(shlex.quote(tmp))

            cmd = "alignment-workflow {} {} {} {} {} {} {} {} {} " \
                  "{} {} {} {}".format(args.seq,
                                       data,
                                       args.config,
                                       sample_name,
                                       do_not_call_str,
                                       logging_str,
                                       tmp_str,
                                       file_str,
                                       slurm_str,
                                       star_str,
                                       flexbar_str,
                                       trim_str,
                                       ribo_cfg_str)

            job_id = slurm.check_sbatch(cmd, args=args, dependencies=all_ribo_jobs)
            all_rna_jobs.append(job_id)

        # collect all jobs and submit to htseq-count

        cdm = "".format()

        slurm.check_sbatch(cmd, args=args, dependencies=all_rna_jobs)


if __name__ == '__main__':
    main()
