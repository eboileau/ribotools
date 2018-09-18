#! /usr/bin/env python3

"""Collect all read length distribution .csv.gz files for all replicates
into a large table (unique counts only). The files must be
available (they are created as part of the pre-processing report.)
"""

import argparse
import logging
import yaml
import csv

import pandas as pd

import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.pandas_utils as pandas_utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as ribo_filenames


logger = logging.getLogger(__name__)


def add_data(sample_name, config):

    note = config.get('note', None)
    read_length_distribution_file = ribo_filenames.get_riboseq_read_length_distribution(
        config['riboseq_data'], sample_name, note=note)
    read_length_distribution = pd.read_csv(read_length_distribution_file)
    # keep only unique reads
    read_length_distribution = read_length_distribution[
        read_length_distribution['basename'].str.contains('unique')]
    ret = pd.pivot_table(read_length_distribution,
                         values='count',
                         columns='length',
                         index='basename')

    return ret


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The yaml config file.")
    parser.add_argument('out', help='''The output file, overwritten by default.
                If path to file does not exist, it will be created.''')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    sample_name_map = ribo_utils.get_sample_name_map(config)

    msg = 'Parsing predictions for replicates.'
    logger.info(msg)

    all_length_distributions = parallel.apply_iter_simple(
        config['riboseq_samples'].keys(),
        add_data,
        config
    )

    all_length_distributions_df = pd.concat(all_length_distributions).fillna(0)
    all_length_distributions_df['name'] = all_length_distributions_df.index
    all_length_distributions_df['name'] = all_length_distributions_df['name'].str.replace('-unique', '')
    all_length_distributions_df['name'] = all_length_distributions_df['name'].apply(lambda x: sample_name_map[x])
    all_length_distributions_df.set_index('name', inplace=True)
    all_length_distributions_df.index.name = None

    msg = "Writing output to: {}".format(args.out)
    logger.info(msg)

    pandas_utils.write_df(all_length_distributions_df, args.out, create_path=True,
                          index=True, sep='\t', header=True,
                          do_not_compress=True, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
