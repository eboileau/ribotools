#! /usr/bin/env python3

"""Provide wrapper for Ribo-seq workflow.

(1) Extract metagene profiles.
(2) Estimate metagene profiles Bayes factors.
(3) Select periodic fragments and offsets.
(4) Optionally, filter non-periodic read lengths from alignment file (BAM)

Note* This is similar to create-orf-profiles from the Rp-Bp
      pipeline, except that it starts from the alignment BAM files,
      and it does not create the ORF profiles, and possibly filters
      the final alignment BAM file.
"""

import sys
import argparse
import logging
import yaml
import shlex

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils
import pbiotools.misc.slurm as slurm

import pbiotools.utils.bam_utils as bam_utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

import ribotools.utils.cl_utils as clu

from rpbp.defaults import default_num_cpus, default_mem, metagene_options

logger = logging.getLogger(__name__)

default_models_base = filenames.get_default_models_base()


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Wrapper script for a Ribo-seq workflow
        starting from alignment files, based on Rp-Bp. File names and directory structures
        follow the conventions used in Rp-Bp.""",
    )

    parser.add_argument("config", help="The yaml config file.")

    parser.add_argument("name", help="The name of the dataset.")

    parser.add_argument(
        "--filter-non-periodic",
        help="""Flag: if this flag is passed,
        non-periodic read lengths will be filtered out of the final (BAM) file available.""",
        action="store_true",
    )

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[get-ribo-periodic]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    # if using slurm, submit the script
    if args.use_slurm:
        cmd = "{}".format(" ".join(shlex.quote(s) for s in sys.argv))
        slurm.check_sbatch(cmd, args=args)
        return

    # check that all of the necessary programs are callable
    programs = [
        "extract-metagene-profiles",
        "estimate-metagene-profile-bayes-factors",
        "select-periodic-offsets",
    ]
    shell_utils.check_programs_exist(programs)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    required_keys = ["riboseq_data", "genome_base_path", "genome_name"]
    utils.check_keys_exist(config, required_keys)

    models_base = config.get("models_base", default_models_base)

    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    # set some general argument strings
    logging_str = logging_utils.get_logging_options_string(args)

    # get config optional arguments
    note = config.get("note", None)
    is_unique = not config.get("keep_riboseq_multimappers", False)

    #  create the metagene profiles
    start_upstream_str = utils.get_config_argument(
        config,
        "metagene_start_upstream",
        "start-upstream",
        default=metagene_options["metagene_start_upstream"],
    )
    start_downstream_str = utils.get_config_argument(
        config,
        "metagene_start_downstream",
        "start-downstream",
        default=metagene_options["metagene_start_downstream"],
    )
    end_upstream_str = utils.get_config_argument(
        config,
        "metagene_end_upstream",
        "end-upstream",
        default=metagene_options["metagene_end_upstream"],
    )
    end_downstream_str = utils.get_config_argument(
        config,
        "metagene_end_downstream",
        "end-downstream",
        default=metagene_options["metagene_end_downstream"],
    )

    riboseq_bam_filename = filenames.get_riboseq_bam(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    metagene_profiles = filenames.get_metagene_profiles(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    # use the canonical transcripts for extracting the metagene profiles
    transcript_bed = filenames.get_bed(
        config["genome_base_path"], config["genome_name"], is_annotated=True
    )

    cmd = "extract-metagene-profiles {} {} {} --num-cpus {} {} {} {} {} {}".format(
        riboseq_bam_filename,
        transcript_bed,
        metagene_profiles,
        args.num_cpus,
        logging_str,
        start_upstream_str,
        start_downstream_str,
        end_upstream_str,
        end_downstream_str,
    )

    in_files = [riboseq_bam_filename, transcript_bed]
    out_files = [metagene_profiles]
    file_checkers = {metagene_profiles: utils.check_gzip_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
        keep_delete_files=keep_delete_files,
    )

    # estimate the periodicity for each offset for all read lengths
    metagene_profile_bayes_factors = filenames.get_metagene_profiles_bayes_factors(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    periodic_models = filenames.get_models(models_base, "periodic")
    non_periodic_models = filenames.get_models(models_base, "nonperiodic")

    periodic_models_str = " ".join(periodic_models)
    non_periodic_models_str = " ".join(non_periodic_models)

    periodic_models_str = "--periodic-models {}".format(periodic_models_str)
    non_periodic_models_str = "--nonperiodic-models {}".format(non_periodic_models_str)

    periodic_offset_start_str = utils.get_config_argument(
        config,
        "periodic_offset_start",
        default=metagene_options["periodic_offset_start"],
    )
    periodic_offset_end_str = utils.get_config_argument(
        config, "periodic_offset_end", default=metagene_options["periodic_offset_end"]
    )
    metagene_profile_length_str = utils.get_config_argument(
        config,
        "metagene_profile_length",
        default=metagene_options["metagene_profile_length"],
    )
    seed_str = utils.get_config_argument(
        config, "seed", default=metagene_options["seed"]
    )
    chains_str = utils.get_config_argument(
        config, "chains", default=metagene_options["chains"]
    )
    iterations_str = utils.get_config_argument(
        config,
        "metagene_iterations",
        "iterations",
        default=metagene_options["metagene_iterations"],
    )

    cmd = (
        "estimate-metagene-profile-bayes-factors {} {} --num-cpus {} {} {} "
        "{} {} {} {} {} {} {}".format(
            metagene_profiles,
            metagene_profile_bayes_factors,
            args.num_cpus,
            periodic_models_str,
            non_periodic_models_str,
            periodic_offset_start_str,
            periodic_offset_end_str,
            metagene_profile_length_str,
            seed_str,
            chains_str,
            iterations_str,
            logging_str,
        )
    )

    in_files = [metagene_profiles]
    in_files.extend(periodic_models)
    in_files.extend(non_periodic_models)
    out_files = [metagene_profile_bayes_factors]
    file_checkers = {metagene_profile_bayes_factors: utils.check_gzip_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
        keep_delete_files=keep_delete_files,
    )

    # select the best read lengths for constructing the signal
    periodic_offsets = filenames.get_periodic_offsets(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    cmd = "select-periodic-offsets {} {}".format(
        metagene_profile_bayes_factors, periodic_offsets
    )

    in_files = [metagene_profile_bayes_factors]
    out_files = [periodic_offsets]
    file_checkers = {periodic_offsets: utils.check_gzip_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
        keep_delete_files=keep_delete_files,
    )

    # filter the last BAM file to retain periodic fragments only
    if not args.filter_non_periodic:
        return

    # get the lengths but we don't use the offsets
    lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
        config, args.name, is_unique=is_unique, default_params=metagene_options
    )

    if len(lengths) == 0:
        msg = (
            "No periodic read lengths and offsets were found. Try relaxing "
            "min_metagene_profile_count, min_metagene_bf_mean, "
            "max_metagene_bf_var, and/or min_metagene_bf_likelihood. Quitting."
        )
        logger.critical(msg)
        return

    lengths_str = " ".join(lengths)

    bam_to_filter = riboseq_bam_filename

    filtered_bam = filenames.get_riboseq_bam(
        config["riboseq_data"],
        args.name,
        is_unique=is_unique,
        length=lengths,
        note=note,
    )

    cmd = "keep-ribo-periodic {} {} --lengths {} ".format(
        bam_to_filter, filtered_bam, lengths_str
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
        keep_delete_files=keep_delete_files,
    )


if __name__ == "__main__":
    main()
