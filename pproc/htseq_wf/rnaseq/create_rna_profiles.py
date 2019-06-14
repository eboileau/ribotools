#! /usr/bin/env python3

""" Get 'RNA profiles', i.e. map reads
to their 5' ends and find their position
with respect to annotated feature (CDS).

"""

import argparse
import logging
import yaml

import misc.logging_utils as logging_utils
import misc.utils as utils
import misc.slurm as slurm

import bio_utils.bam_utils as bam_utils

import riboutils.ribo_filenames as filenames
import riboutils.ribo_utils as ribo_utils

import btea.utils.cl_utils as clu

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The rna yaml config file.")
    parser.add_argument('ribo_config', help="The ribo config file for rna trimmed data.")
    parser.add_argument('bed', help="""The annotated transcripts file in BED format.""")

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser)
    clu.add_isoform_strategy(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))

    config_keys = ['rnaseq_data', 'rnaseq_samples', 'matching_samples']
    utils.check_keys_exist(config, config_keys)

    is_unique = not ('keep_rnaseq_multimappers' in config)
    note = config.get('note', None)

    ribo_config = yaml.load(open(args.ribo_config))
    is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)

    for sample_name in config['rnaseq_samples'].keys():

        matching_ribo_sample = config['matching_samples'][sample_name]
        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            ribo_config,
            matching_ribo_sample,
            is_unique=is_unique_ribo,
            isoform_strategy=args.isoform_strategy
        )

        if len(lengths) == 0:
            msg = "No periodic read lengths and offsets were found!"
            logger.critical(msg)
            return

        lengths = str(max([int(l) for l in lengths]))

        genome_bam = filenames.get_rnaseq_bam(
            config['rnaseq_data'],
            sample_name,
            is_unique=is_unique,
            length=lengths,
            isoform_strategy=args.isoform_strategy,
            note=note
        )

        # first index file
        # bam_utils.index_bam_file(genome_bam, args)

        # create the metagene profiles, use riboseq nomenclature
        metagene_profiles = filenames.get_riboseq_base(config['rnaseq_data'], sample_name,
            'metagene-profiles', length=lengths, is_unique=is_unique,
            isoform_strategy=args.isoform_strategy, note=note)
        metagene_profiles = metagene_profiles + ".metagene-profile.csv.gz"

        start_upstream_str = utils.get_config_argument(config,
            'metagene_profile_start_upstream', 'start-upstream')
        start_downstream_str = utils.get_config_argument(config,
            'metagene_profile_start_downstream', 'start-downstream')
        end_upstream_str = utils.get_config_argument(config,
            'metagene_profile_end_upstream', 'end-upstream')
        end_downstream_str = utils.get_config_argument(config,
            'metagene_profile_end_downstream', 'end-downstream')

        # cmd = ("extract-metagene-profiles {} {} {} --num-cpus {} {} {} {} {} {}".format(
        #     genome_bam, args.bed, metagene_profiles, args.num_cpus, logging_str,
        #     start_upstream_str, start_downstream_str, end_upstream_str, end_downstream_str))

        cmd = ("/home/eboileau/devc/ext/extract_rna_profiles.py {} {} {} --num-cpus {} {} {} {} {} {}".format(
            genome_bam, args.bed, metagene_profiles, args.num_cpus, logging_str,
            start_upstream_str, start_downstream_str, end_upstream_str, end_downstream_str))

        slurm.check_sbatch(cmd, args=args)

if __name__ == '__main__':
    main()
