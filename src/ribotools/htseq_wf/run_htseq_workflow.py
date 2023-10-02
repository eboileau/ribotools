#! /usr/bin/env python3

"""Helper script to submit a set of RNA- and/or Ribo-seq
samples for alignment and abundance quantification (using HTSeq).

Note* Example workflow: Submit Ribo-seq samples for
      mapping and periodicity estimation (or else use
      data available from the Rp-Bp pipeline). When
      completed, submit RNA-seq samples, trimming
      reads to max periodic fragment length of matching
      Ribo-seq samples.

      Ribo-seq: Estimate read length periodicity,
      and filter out non-periodic fragments from the
      final BAM files (default). If periodicity estimates
      and/or mapped reads were previously obtained by running
      Rp-Bp, they must be available. Reads are NOT P-site
      offset-shifted.

      RNA-seq: Reads can be trimmed to max periodic fragment
      length from the matching Ribo-seq data before mapping,
      using the same default mapping options. Ribo-seq periodicity
      estimates must be available.

      HTSeq: Deals with library strandedness. All options are passed
      via this script.
"""

import os
import sys
import argparse
import logging
import yaml
import shlex
import re

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils
import pbiotools.misc.slurm as slurm

import pbiotools.utils.pgrm_utils as pgrm_utils

import ribotools.utils.cl_utils as clu

from rpbp.defaults import default_num_cpus, default_mem, star_executable
from ribotools.defaults import htseq_options

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Submit a set of samples to RNA- or
        Ribo-seq alignment pipelines. The pipeline is called for every sample in the
        configuration file.""",
    )

    parser.add_argument("seq", choices=["rna", "ribo"])

    parser.add_argument("config", help="The yaml configuration file.")

    parser.add_argument(
        "--skip-periodicity-estimation",
        help="""Flag: by default,
        estimate periodicity and filter non-periodic read lengths from the final alignment
        files. For Ribo-seq only.""",
        action="store_true",
    )

    parser.add_argument(
        "--run-all",
        help="""Flag: map and count Ribo-seq and RNA-seq, one
        after the other, provided that ALL Ribo-seq jobs completed successfully. The same
        general options are used, including options used for Flexbar, STAR and Bowtie2.
        For Ribo-seq only.""",
        action="store_true",
        required="--stranded" in sys.argv,
    )

    parser.add_argument(
        "--stranded",
        help="""Optional argument: library strandedness
        for RNA-seq. This option is passed to htseq-count and overrides the same option
        passed via [--htseq-options] and used for Ribo-seq. Unless given, the default
        value will be used.""",
        choices=["yes", "reverse", "no"],
        default="no",
    )

    parser.add_argument(
        "--trim-rna-to-max-fragment-size",
        help="""Flag: trim RNA post
        adapter removal using max fragment size from the matching Ribo-seq sample. Note* At least
        the "periodic-offsets" file must be available. The config file must also include
        "matching_samples" and the path to the Ribo-seq config must be given [--ribo-config])""",
        action="store_true",
    )

    parser.add_argument(
        "--ribo-config",
        help="""Optional argument: the Ribo-seq config file
        if using [--trim-rna-to-max-fragment-size].""",
        required="--trim-rna-to-max-fragment-size" in sys.argv,
        type=str,
    )

    parser.add_argument(
        "--rna-config",
        help="""Optional argument: the RNA-seq config file
            if using [--run-all].""",
        required="--run-all" in sys.argv,
        type=str,
    )

    parser.add_argument(
        "--gtf",
        help="""A different GTF file for abundance estimation, e.g. the output
        of get-gtf-from-predictions (Ribo-seq ORFs). This is passed to
        htseq-count and overrides the GTF file from the config.""",
        type=str,
        dest="htseq_gtf",
    )

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
        "flexbar",
        "bowtie2",
        args.star_executable,
        "samtools",
        "alignment-workflow",
        "get-ribo-periodic",
        "keep-ribo-periodic",
        "call-htseq-count",
    ]
    shell_utils.check_programs_exist(programs)

    # handle all option strings to call the pipeline script
    logging_str = logging_utils.get_logging_options_string(args)
    file_str = clu.get_file_options_string(args)
    slurm_str = slurm.get_slurm_options_string(args)
    star_str = pgrm_utils.get_star_options_string(args)
    flexbar_str = pgrm_utils.get_flexbar_options_string(args)
    htseq_str = clu.get_htseq_options_string(args)
    htseq_gtf_str = ""
    if args.htseq_gtf:
        htseq_gtf_str = "--gtf {}".format(args.htseq_gtf)
    tmp_str = ""

    # handle do_not_call so that we do call the pipeline script, but that it does not run anything
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"
    args.do_not_call = False

    if args.seq.lower() == "rna" and (args.skip_periodicity_estimation or args.run_all):
        msg = """seq is RNA and either [--skip-periodicity-estimation] or [--run-all]
            were given. These options will be ignored."""
        logger.warning(msg)

    if not args.run_all and args.rna_config:
        msg = """The [--rna-config] option is passed without
        [--run-all], it will be ignored."""
        logger.warning(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    base_keys = [
        "ribosomal_index",
        "star_index",
        "genome_base_path",
        "genome_name",
        "fasta",
        "gtf",
    ]

    all_ribo_jobs = None

    # now process samples
    if args.seq.lower() == "ribo":

        config_keys = base_keys + ["riboseq_data", "riboseq_samples"]
        utils.check_keys_exist(config, config_keys)

        # first submit all alignment jobs, then estimate periodicity
        # for each one
        job_ids_mapping = {}
        for sample_name, data in config["riboseq_samples"].items():

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
                flexbar_str,
            )

            job_id = slurm.check_sbatch(cmd, args=args)
            job_ids_mapping[sample_name] = job_id

        job_ids_periodic = {}
        not_periodic_str = "--not-periodic"

        if not args.skip_periodicity_estimation:

            not_periodic_str = ""
            filter_non_periodic_str = "--filter-non-periodic"

            for sample_name in config["riboseq_samples"].keys():

                cmd = "get-ribo-periodic {} {} {} {} {} {} {}".format(
                    args.config,
                    sample_name,
                    filter_non_periodic_str,
                    do_not_call_str,
                    logging_str,
                    file_str,
                    slurm_str,
                )

                job_id = [job_ids_mapping[sample_name]]
                job_id_periodic = slurm.check_sbatch(
                    cmd, args=args, dependencies=job_id
                )
                job_ids_periodic[sample_name] = job_id_periodic

        # run htseq-count
        for sample_name in config["riboseq_samples"].keys():

            cmd = "call-htseq-count {} {} {} {} {} {} {} {} {}".format(
                args.seq,
                args.config,
                sample_name,
                logging_str,
                file_str,
                slurm_str,
                htseq_str,
                not_periodic_str,
                htseq_gtf_str,
            )

            job_id = [
                job_ids_periodic[sample_name]
                if job_ids_periodic.get(sample_name, None)
                else job_ids_mapping[sample_name]
            ]
            slurm.check_sbatch(cmd, args=args, dependencies=job_id)

        # we are done with ribo, however if we run all,
        # we need to collect all jobs and check htseq-options for stranded
        if args.run_all:

            all_ribo_jobs = list(job_ids_mapping.values())
            all_ribo_jobs.extend(list(job_ids_periodic.values()))

            if args.htseq_options is not None:
                # replace option used for ribo
                args.htseq_options = [
                    opt
                    if "--stranded" not in opt
                    else re.sub("yes|no|reverse", args.stranded, opt)
                    for opt in args.htseq_options
                ]
            else:
                # use default if no arguments are passed
                args.htseq_options = ["--stranded {}".format(args.stranded)]

            htseq_str = clu.get_htseq_options_string(args)

            # and pass the right config file
            args.config = args.rna_config
            config = yaml.load(open(args.config), Loader=yaml.FullLoader)

            args.seq = "rna"

    if args.seq.lower() == "rna":

        config_keys = base_keys + ["rnaseq_data", "rnaseq_samples"]

        trim_str = ""
        ribo_cfg_str = ""

        if args.trim_rna_to_max_fragment_size:

            config_keys.extend(["matching_samples"])

            if "--post-trim-length" in flexbar_str:
                msg = (
                    "Flexbar option [--post-trim-length] already given..."
                    "will ignore [trim-rna-to-max-fragment-size]."
                )
                logger.critical(msg)
            else:
                trim_str = "--trim-rna-to-max-fragment-size"
                ribo_cfg_str = "--ribo-config {}".format(shlex.quote(args.ribo_config))

        if not args.trim_rna_to_max_fragment_size and args.ribo_config:
            msg = """The [--ribo-config] option is passed without
            [--trim-rna-to-max-fragment-size], it will be ignored."""
            logger.warning(msg)

        utils.check_keys_exist(config, config_keys)

        job_ids_mapping = {}
        for sample_name, data in config["rnaseq_samples"].items():

            if args.tmp is not None:
                tmp = os.path.join(args.tmp, "{}_htseq_wf_rna".format(sample_name))
                tmp_str = "--tmp {}".format(shlex.quote(tmp))

            cmd = "alignment-workflow {} {} {} {} {} {} {} {} {} " "{} {} {} {}".format(
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
                flexbar_str,
                trim_str,
                ribo_cfg_str,
            )

            job_id = slurm.check_sbatch(cmd, args=args, dependencies=all_ribo_jobs)
            job_ids_mapping[sample_name] = job_id

        # run htseq-count
        for sample_name in config["rnaseq_samples"].keys():

            cmd = "call-htseq-count {} {} {} {} {} {} {} {} {}".format(
                args.seq,
                args.config,
                sample_name,
                logging_str,
                file_str,
                slurm_str,
                htseq_str,
                ribo_cfg_str,
                htseq_gtf_str,
            )

            job_id = [job_ids_mapping[sample_name]]
            slurm.check_sbatch(cmd, args=args, dependencies=job_id)


if __name__ == "__main__":
    main()
