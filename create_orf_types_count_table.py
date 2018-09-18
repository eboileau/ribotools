#! /usr/bin/env python3

import argparse
import logging
import os
import yaml
import csv

import pandas as pd

import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.utils as utils
import misc.pandas_utils

import bio_utils.bed_utils as bed_utils

import riboutils.ribo_filenames as filenames
import riboutils.ribo_utils as ribo_utils

logger = logging.getLogger(__name__)


def get_orf_type_counts(name, is_single_sample, name_map, config, args):
    
    note_str = config.get('note', None)
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)
    is_unique = not ('keep_riboseq_multimappers' in config)
    
    if is_single_sample:
        try:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
                config, 
                name,
                is_unique=is_unique
            )
        except FileNotFoundError:
            msg = "Could not find metagene periodicity file. Skipping. name: {}".format(name,)
            logger.warning(msg)
            return None
    else:
        lengths, offsets = None, None
        
    predicted_orfs = filenames.get_riboseq_predicted_orfs(
        config['riboseq_data'], 
        name, 
        length=lengths, 
        offset=offsets, 
        is_unique=is_unique, 
        note=note_str, 
        fraction=fraction, 
        reweighting_iterations=reweighting_iterations,
        is_filtered=True, 
        is_chisq=False
    )
    
    if not os.path.exists(predicted_orfs):
        msg = "Could not find predicted ORFs. name: {}. file: {}".format(name, predicted_orfs)
        logger.warning(msg)
        return None
    
    bed = bed_utils.read_bed(predicted_orfs)

    if args.use_labels:
        bed['orf_type_group'] = bed['orf_type'].map(
            ribo_utils.orf_type_labels_reverse_mapping)
        orf_type_counts = bed.groupby(['orf_type_group', 'strand']).size()
        orf_type_counts = orf_type_counts.reset_index(name="count")
        orf_type_counts['category'] = orf_type_counts['orf_type_group'].map(
            ribo_utils.orf_type_labels_display_name_map)
    else:
        orf_type_counts = bed.groupby(['orf_type', 'strand']).size()
        orf_type_counts = orf_type_counts.reset_index(name="count")
        orf_type_counts['category'] = orf_type_counts['orf_type'].map(
            ribo_utils.orf_type_display_name_map)

    orf_type_counts['source'] = name_map[name]

    return orf_type_counts


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Parse predictions and create ORF types count table.")

    parser.add_argument('config', help="The (yaml) config file.")
    parser.add_argument('out', help='''The output file (gz), overwritten by default.
        If path to file does not exist, it will be created.''')

    parser.add_argument('--use-labels', help='''If this flag is present, then 
        group similar ORF types using label mappings defined in "ribo_utils".''',
                        action='store_true')
    parser.add_argument('--use-conditions', help='''If this flag is present, then the 
        count table includes results from the merged replicates in addition to those 
        from each of the individual samples, else only results from single samples 
        are used.''', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    sample_name_map = ribo_utils.get_sample_name_map(config)
    condition_name_map = ribo_utils.get_condition_name_map(config)

    msg = 'Parsing predictions for replicates.'
    logger.info(msg)

    is_single_sample = True
    sample_orf_types = parallel.apply_iter_simple(
        config['riboseq_samples'].keys(),
        get_orf_type_counts,
        is_single_sample,
        sample_name_map,
        config,
        args
        )
    sample_orf_types = utils.remove_nones(sample_orf_types)
    
    if args.use_conditions:
        msg = 'Parsing predictions for conditions.'
        logger.info(msg)
        is_single_sample = False
        merged_sample_orf_types = parallel.apply_iter_simple(
            ribo_utils.get_riboseq_replicates(config),
            get_orf_type_counts,
            is_single_sample,
            condition_name_map,
            config,
            args
        )
        merged_sample_orf_types = utils.remove_nones(merged_sample_orf_types)
        sample_orf_types = sample_orf_types + merged_sample_orf_types
    
    sample_orf_types_df = pd.concat(sample_orf_types)

    msg = "Writing output to: {}".format(args.out)
    logger.info(msg)

    misc.pandas_utils.write_df(sample_orf_types_df, args.out, create_path=True,
                               index=False, sep='\t', header=True,
                               do_not_compress=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
