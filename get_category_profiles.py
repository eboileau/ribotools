#! /usr/bin/env python3

""" Get 'metagene profiles' for all ribo samples,
but from different groups of transcripts, these
are given in a file.
"""

import os
import argparse
import logging
import yaml

import pandas as pd

import misc.logging_utils as logging_utils
import misc.utils as utils
import misc.slurm as slurm

import riboutils.ribo_filenames as filenames

import btea.utils.cl_utils as clu

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The yaml config file.")
    parser.add_argument('bed_list', help="""List of annotated transcript files
        to extract metagene profiles for each category (could be only CDS, or else
        uORF, etc). Full path to each file must be given, as well as a name for a
        subdirectory where profiles will be saved (tab separated, no header). 
        Files must be in BED12 format.""")

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))

    config_keys = ['riboseq_data', 'riboseq_samples']
    utils.check_keys_exist(config, config_keys)

    is_unique = not ('keep_riboseq_multimappers' in config)
    note = config.get('note', None)

    # read categories
    bed_list = pd.read_csv(args.bed_list, sep='\t', header=None, names=['bed_file', 'dirloc'])

    # base path
    base_path = os.path.join(config['riboseq_data'], 'metagene-category-profiles')
    os.makedirs(base_path, exist_ok=True)

    for row in bed_list.itertuples():

        bed = row.bed_file
        dirloc = row.dirloc

        msg = 'Creating profiles for {}'.format(bed)
        logger.info(msg)

        # check output path
        output_path = os.path.join(base_path, dirloc)
        os.makedirs(output_path, exist_ok=True)

        for sample_name in config['riboseq_samples'].keys():

            riboseq_bam_filename = filenames.get_riboseq_bam(config['riboseq_data'],
                                                             sample_name,
                                                             is_unique=is_unique,
                                                             note=note)

            # create the metagene profiles, use riboseq nomenclature
            # to specific locations
            metagene_profiles = filenames.get_riboseq_base(config['riboseq_data'], sample_name,
                output_path, is_unique=is_unique, note=note)
            metagene_profiles = metagene_profiles + ".metagene-profile.csv.gz"

            start_upstream_str = utils.get_config_argument(config,
                'metagene_profile_start_upstream', 'start-upstream')
            start_downstream_str = utils.get_config_argument(config,
                'metagene_profile_start_downstream', 'start-downstream')
            end_upstream_str = utils.get_config_argument(config,
                'metagene_profile_end_upstream', 'end-upstream')
            end_downstream_str = utils.get_config_argument(config,
                'metagene_profile_end_downstream', 'end-downstream')

            cmd = ("extract-metagene-profiles {} {} {} --num-cpus {} {} {} {} {} {}"
                   .format(riboseq_bam_filename, bed, metagene_profiles,
                           args.num_cpus, logging_str, start_upstream_str,
                           start_downstream_str, end_upstream_str, end_downstream_str))


            slurm.check_sbatch(cmd, args=args)



if __name__ == '__main__':
    main()
