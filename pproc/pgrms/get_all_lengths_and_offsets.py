#! /usr/bin/env python3

"""Get all lengths/offsets for all Ribo-seq samples
in the configuration file and output these to a dataframe.

Functions:
    get_lengths_and_offsets_call
"""

import os
import argparse
import logging
import yaml
import csv

import pandas as pd

import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.pandas_utils as pandas_utils

import pbio.ribo.ribo_utils as ribo_utils

from rpbp.defaults import metagene_options

logger = logging.getLogger(__name__)


def get_lengths_and_offsets_call(sample_name, config, args):

    is_unique = not ('keep_riboseq_multimappers' in config)

    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
        config,
        sample_name,
        is_unique=is_unique,
        default_params=metagene_options
    )

    if len(lengths) == 0:
        msg = "No periodic read lengths and offsets were found!"
        logger.critical(msg)
        return

    sample_name_map = ribo_utils.get_sample_name_map(config)

    lengths = [int(l) for l in lengths]
    offsets = [int(o) for o in offsets]
    data = {'condition': sample_name_map[sample_name], 'lengths': lengths, 'offsets': offsets}

    return pd.DataFrame(data, columns=['condition', 'lengths', 'offsets'])


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The yaml config file.")

    parser.add_argument('out', help='''The (.csv.gz) output file complete path.''')

    parser.add_argument('--overwrite', help='''Overwrites output.''', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    msg = 'Parsing predictions for replicates.'
    logger.info(msg)

    all_lengths_and_offsets = parallel.apply_iter_simple(
        config['riboseq_samples'].keys(),
        get_lengths_and_offsets_call,
        config,
        args
    )

    all_lengths_and_offsets_df = pd.concat(all_lengths_and_offsets)

    msg = "Writing output to: {}".format(args.out)
    logger.info(msg)

    if os.path.exists(args.out) and not args.overwrite:
        msg = "Output file {} already exists. Skipping.".format(args.out)
        logger.warning(msg)
    else:
        pandas_utils.write_df(all_lengths_and_offsets_df, args.out, create_path=True,
                              index=False, sep=',', header=True,
                              do_not_compress=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
