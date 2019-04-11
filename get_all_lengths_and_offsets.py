#! /usr/bin/env python3

"""Get all lengths/offsets for all ribo-seq samples
in the configuration file and output these to a dataframe.
"""

import os
import argparse
import logging
import yaml
import csv

import pandas as pd

import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.pandas_utils

import riboutils.ribo_utils as ribo_utils

import btea.utils.cl_utils as clu

logger = logging.getLogger(__name__)


def get_lengths_and_offsets_call(sample_name, args):

    config = yaml.load(open(args.config))
    is_unique = not ('keep_riboseq_multimappers' in config)

    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
        config,
        sample_name,
        is_unique=is_unique,
        isoform_strategy=args.isoform_strategy
    )

    if len(lengths) == 0:
        msg = "No periodic read lengths and offsets were found!"
        logger.critical(msg)
        return

    sample_name_map = ribo_utils.get_sample_name_map(config)

    lengths = [int(l) for l in lengths]
    offsets = [int(o) for o in offsets]
    data = {'name': sample_name_map[sample_name], 'lengths': lengths, 'offsets': offsets}
    return pd.DataFrame(data, columns=['name', 'lengths', 'offsets'])


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The yaml config file.")
    parser.add_argument('out', help='''The output file, overwritten by default.''')

    parser.add_argument('--note', default=None)

    clu.add_isoform_strategy(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    msg = 'Parsing predictions for replicates.'
    logger.info(msg)

    all_lengths_and_offsets = parallel.apply_iter_simple(
        config['riboseq_samples'].keys(),
        get_lengths_and_offsets_call,
        args
    )

    all_lengths_and_offsets_df = pd.concat(all_lengths_and_offsets)

    msg = "Writing output to: {}".format(args.out)
    logger.info(msg)

    misc.pandas_utils.write_df(all_lengths_and_offsets_df, args.out, create_path=True,
                               index=False, sep='\t', header=True,
                               do_not_compress=True, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
