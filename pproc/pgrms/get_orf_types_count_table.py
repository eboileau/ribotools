#! /usr/bin/env python3

"""Get ORF count table from the merged predictions.

Functions:
    get_orf_type_counts
"""

import os
import argparse
import logging
import csv

import pbio.misc.logging_utils as logging_utils
import pbio.misc.pandas_utils as pandas_utils

import pbio.utils.bed_utils as bed_utils

import pbio.ribo.ribo_utils as ribo_utils

logger = logging.getLogger(__name__)


def get_orf_type_counts(condition):

    orf_type_counts = condition.groupby(['orf_type', 'strand']).size()
    orf_type_counts = orf_type_counts.reset_index(name="count")

    return orf_type_counts


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Parse predictions and create ORF types count table.")

    parser.add_argument('orfs', help="The combined (long format) ORF predictions complete path")

    parser.add_argument('out', help='''The (.csv.gz) output file complete path.''')

    parser.add_argument('--group-labels', help='''Types are substituted by the broader
        categories defined in ribo_utils.''', action='store_true')

    parser.add_argument('--overwrite', help='''Overwrites output.''', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    
    orfs = bed_utils.read_bed(args.orfs)

    if args.group_labels:
        orfs['orf_type'] = orfs['orf_type'].map(ribo_utils.orf_type_labels_reverse_mapping)

    orf_types = orfs.groupby('condition').apply(get_orf_type_counts)
    orf_types.reset_index(inplace=True)
    orf_types.rename(columns={'orf_type': 'ORF category'}, inplace=True)
    orf_types.drop(columns='level_1', inplace=True)

    if os.path.exists(args.out) and not args.overwrite:
        msg = "Output file {} already exists. Skipping.".format(args.out)
        logger.warning(msg)
    else:
        pandas_utils.write_df(orf_types, args.out, create_path=True,
                              index=False, sep=',', header=True,
                              do_not_compress=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
