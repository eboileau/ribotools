#! /usr/bin/env python3

"""Wrapper for 'keep-ribo-periodic', to extract periodic
periodic reads from existing BAM files (Rp-Bp pipeline).
See 'run-htseq-workflow' for more details. No slurm option."""


import argparse
import logging
import yaml

import pbiotool.misc.logging_utils as logging_utils
import pbiotool.misc.shell_utils as shell_utils

import pbiotool.utils.bam_utils as bam_utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

import pproc.utils.cl_utils as clu

from rpbp.defaults import metagene_options

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract periodic reads for a set of
        existing BAM files. The periodic lengths and offsets file must be available.""",
    )

    parser.add_argument("config", help="The yaml configuration file.")

    parser.add_argument(
        "--do-not-call",
        help="""If this flag is present, then the program
        will not be executed (but all commands will be printed).""",
        action="store_true",
    )

    clu.add_file_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs = ["samtools", "keep-ribo-periodic"]
    shell_utils.check_programs_exist(programs)

    # handle all option strings to call the pipeline script
    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # get config optional arguments
    note = config.get("note", None)
    is_unique = not ("keep_riboseq_multimappers" in config)

    for sample_name in config["riboseq_samples"].keys():

        bam_to_filter = filenames.get_riboseq_bam(
            config["riboseq_data"], sample_name, is_unique=is_unique, note=note
        )

        # get the lengths but we don't use the offsets
        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            config, sample_name, is_unique=is_unique, default_params=metagene_options
        )

        if len(lengths) == 0:
            msg = (
                "No periodic read lengths and offsets were found! Are the "
                "periodic lengths and offsets files available?"
            )
            logger.critical(msg)
            return

        lengths_str = " ".join(lengths)

        filtered_bam = filenames.get_riboseq_bam(
            config["riboseq_data"],
            sample_name,
            is_unique=is_unique,
            length=lengths,
            note=note,
        )

        cmd = "keep-ribo-periodic {} {} --lengths {} {} ".format(
            bam_to_filter, filtered_bam, lengths_str, logging_str
        )
        in_files = [bam_to_filter]
        out_files = [filtered_bam]
        file_checkers = {filtered_bam: bam_utils.check_bam_file}
        shell_utils.call_if_not_exists(
            cmd,
            out_files,
            in_files=in_files,
            file_checkers=file_checkers,
            overwrite=args.overwrite,
            call=call,
        )


if __name__ == "__main__":
    main()
