#! /usr/bin/env python3

""" Get 'RNA profiles', i.e. map reads
to their 5' ends and find their position
with respect to annotated feature (CDS).
"""

import sys
import argparse
import logging
import yaml

import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.utils as utils
import pbio.misc.slurm as slurm

import pbio.utils.bam_utils as bam_utils

import pbio.ribo.ribo_filenames as filenames
import pbio.ribo.ribo_utils as ribo_utils

from pproc.defaults import default_num_cpus, default_mem, metagene_options

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The rna yaml config file.")
    parser.add_argument('ribo_config', help="The ribo config file for rna trimmed data.")
    parser.add_argument('bed', help="""The annotated transcripts file in BED format.""")

    parser.add_argument('-t', '--tmp', help="""Optional argument: where to write 
                        temporary files. If not specified, programs-specific tmp will be used.""", 
                        default=None)

    parser.add_argument('--overwrite', help="""Flag: overwrite existing files.""", 
                        action='store_true')

    parser.add_argument('-k', '--keep-intermediate-files', help="""Flag: unless this flag is 
                        given, all intermediate files (such as discarded reads) will be deleted, 
                        unless the [--do-not-call] option is also given.""", 
                        action='store_true')
    
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    logging_str = logging_utils.get_logging_options_string(args)
    
    msg = "[create-rna-profiles]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    config_keys = ['rnaseq_data', 'rnaseq_samples', 'matching_samples']
    utils.check_keys_exist(config, config_keys)

    is_unique = not ('keep_rnaseq_multimappers' in config)
    note = config.get('note', None)

    ribo_config = yaml.load(open(args.ribo_config), Loader=yaml.FullLoader)
    is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)
    
    # check that all of the necessary programs are callable
    programs = ['extract-rna-profiles']
    shell_utils.check_programs_exist(programs)
    

    for sample_name in config['rnaseq_samples'].keys():

        matching_ribo_sample = config['matching_samples'][sample_name]
        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            ribo_config,
            matching_ribo_sample,
            is_unique=is_unique_ribo,
            default_params=metagene_options
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
            note=note
        )

        # create the metagene profiles, use riboseq nomenclature
        metagene_profiles = filenames.get_riboseq_base(config['rnaseq_data'], 
                                                       sample_name,
                                                       'metagene-profiles', 
                                                       length=lengths, 
                                                       is_unique=is_unique,
                                                       note=note)
        metagene_profiles = metagene_profiles + ".metagene-profile.csv.gz"

        start_upstream_str = utils.get_config_argument(config,
                                                       'metagene_start_upstream',
                                                       'start-upstream',
                                                       default=metagene_options['metagene_start_upstream'])
        start_downstream_str = utils.get_config_argument(config,
                                                         'metagene_start_downstream',
                                                         'start-downstream',
                                                         default=metagene_options['metagene_start_downstream'])
        end_upstream_str = utils.get_config_argument(config,
                                                     'metagene_end_upstream',
                                                     'end-upstream',
                                                     default=metagene_options['metagene_end_upstream'])
        end_downstream_str = utils.get_config_argument(config,
                                                       'metagene_end_downstream',
                                                       'end-downstream',
                                                       default=metagene_options['metagene_end_downstream'])

        cmd = ("extract-rna-profiles {} {} {} --num-cpus {} {} {} {} {} {}".format(genome_bam, 
                                                                                   args.bed, 
                                                                                   metagene_profiles, 
                                                                                   args.num_cpus, 
                                                                                   logging_str,
                                                                                   start_upstream_str, 
                                                                                   start_downstream_str, 
                                                                                   end_upstream_str, 
                                                                                   end_downstream_str))

        slurm.check_sbatch(cmd, args=args)


if __name__ == '__main__':
    main()
