#! /usr/bin/env python3

import os
import argparse
import logging
import pandas as pd
import yaml
import csv

import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.pandas_utils as pandas_utils

import bio_utils.bed_utils as bed_utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as filenames

logger = logging.getLogger(__name__)

selected_fields = bed_utils.bed12_field_names + ['orf_num', 'orf_len', 'orf_type',
                                                 'bayes_factor_mean', 'bayes_factor_var',
                                                 'x_1_sum', 'x_2_sum', 'x_3_sum']

selected_fields_name_map = {'bayes_factor_mean': 'BF_mean', 'bayes_factor_var': 'BF_var',
                            'x_1_sum': 'frame1_count', 'x_2_sum': 'frame2_count',
                            'x_3_sum': 'frame3_count'}


def get_predictions_file(name, is_sample, config):

    note_str = config.get('note', None)
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)
    is_unique = not ('keep_riboseq_multimappers' in config)

    if is_sample:
        try:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
                config, name, is_unique=is_unique)
        except FileNotFoundError:
            msg = "Could not find periodic lengths and offset file. Skipping. name: {}".format(name)
            logger.warning(msg)
            return None
    else:
        lengths, offsets = None, None

    predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], name,
                                                          length=lengths, offset=offsets,
                                                          is_unique=is_unique, note=note_str,
                                                          fraction=fraction,
                                                          reweighting_iterations=reweighting_iterations,
                                                          is_filtered=True)

    if not os.path.exists(predicted_orfs):
        msg = "Could not find predicted ORFs. name: {}. file: {}".format(name, predicted_orfs)
        raise FileNotFoundError(msg)

    return predicted_orfs


def add_data(name, sample_name_map, is_sample, config):

    orfs_file = get_predictions_file(name, is_sample, config)
    orfs = bed_utils.read_bed(orfs_file)
    orfs = orfs[selected_fields]
    orfs['orf_type'] = orfs['orf_type'].map(
        ribo_utils.orf_type_display_name_map)
    orfs['source'] = sample_name_map[name]

    return orfs


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''Parse predictions and create unique output 
                                     table with selected fields.''')

    parser.add_argument('config', help="The (yaml) config file.")
    parser.add_argument('out', help="The BED12+ output file name (gz).")

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    sample_name_map = ribo_utils.get_sample_name_map(config)
    condition_name_map = ribo_utils.get_condition_name_map(config)

    msg = 'Parsing predictions for replicates.'
    logger.info(msg)

    is_sample = True
    all_predictions = parallel.apply_iter_simple(
        config['riboseq_samples'].keys(),
        add_data,
        sample_name_map,
        is_sample,
        config
    )

    if 'riboseq_biological_replicates' in config:
        msg = 'Parsing predictions for conditions.'
        logger.info(msg)
        is_sample = False
        condition_predictions = parallel.apply_iter_simple(
            ribo_utils.get_riboseq_replicates(config).keys(),
            add_data,
            condition_name_map,
            is_sample,
            config
        )
        all_predictions = all_predictions + condition_predictions

    all_predictions_df = pd.concat(all_predictions)
    all_predictions_df.rename(columns=selected_fields_name_map, inplace=True)

    msg = "Writing output to: {}".format(args.out, )
    logger.info(msg)

    header = ['#{}'.format(c) for c in all_predictions_df.columns]
    pandas_utils.write_df(all_predictions_df, args.out, create_path=True,
                          index=False, sep='\t',
                          header=header, do_not_compress=False,
                          quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
