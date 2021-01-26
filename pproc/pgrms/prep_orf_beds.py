#! /usr/bin/env python3

"""Prepare the Rp-Bp ORF predictions BED files to be used 
with trackhub-utils (get-bed2bigBed).

Note: The header must match those of bed_utils.bed12_field_names

Functions:
    _get_bed
"""

import sys
import os
import argparse
import logging
import yaml
import csv
import pandas as pd

import pbio.misc.utils as utils
import pbio.misc.pandas_utils as pandas_utils
import pbio.misc.logging_utils as logging_utils

import pbio.utils.bed_utils as bed_utils

import pbio.ribo.ribo_utils as ribo_utils
import pbio.ribo.ribo_filenames as filenames

from rpbp.defaults import metagene_options

logger = logging.getLogger(__name__)



ORF_FIELDS = bed_utils.bed12_field_names + ['orf_num', 'orf_len', 'orf_type',
                                            'bayes_factor_mean', 'bayes_factor_var',
                                            'x_1_sum', 'x_2_sum', 'x_3_sum']

# use the base labels, but for display we need to adjust
display_name_map = ribo_utils.orf_type_labels_display_name_map
display_name_map['canonical_variant'] = 'Variant'

color_mapping = {
    'Canonical': '31,119,181',  # #1f77b4
    'Variant': '158,218,229',  # #9edae5
    'ncORF': '197,176,213',  # #c5b0d5
    'uORF': '152,223,138',  # #98df8a
    'dORF': '219,219,141',  # #dbdb8d
    'Novel': '255,187,120',  # #ffbb78
    'Other': '255,0,0'
}



def _get_bed(input_filename, output_filename, args):

    """Get BED12+ file and adjust features. The fields must match those
    defined in bed_utils.bed12_field_names and extra fields from rpbp.
    """

    bed_df = bed_utils.get_bed_df(input_filename)
    
    msg = "No. of unique features by 'id': {}".format(len(bed_df['id'].unique()))
    logger.info(msg)
    
    # first keep only selected ORFs/features
    if args.id_list is not None:
        bed_df = bed_df[bed_df['id'].isin(args.id_list)]
        
    # filter by p-sites: final list may include ORFs satisfying the criteria,
    # but not for a given replicate (e.g. >=10 p-sites in 3 replicates, but not
    # all replicates)
    if args.filter_by_psites:
        bed_df = bed_df[bed_df['x_1_sum'] >= args.psite_threshold]

    # add display label, make sure orf_category is also in fields
    for orf_type, labels in ribo_utils.orf_type_labels_mapping.items():
        bed_df.loc[bed_df['orf_type'].isin(labels), 'orf_category'] = display_name_map[orf_type]

    # just to make sure...
    remove_m = bed_df['orf_type'].isna()
    if not args.keep_other:
        other_m = bed_df['orf_category'] == 'Other'
        remove_m = remove_m | other_m
    bed_df = bed_df[~remove_m]

    msg = "No. of unique features (after filtering): {}".format(len(bed_df['id'].unique()))
    logger.info(msg)
    
    # convert counts to int
    bed_df = bed_df.astype({"x_1_sum": int, "x_2_sum": int, "x_3_sum": int})
    
    if not args.no_color:
        for label, color in color_mapping.items():
            label_m = bed_df['orf_category'] == label
            bed_df.loc[label_m, 'color'] = color
    else:
        pass

    # Sort on the chrom field, and then on the chromStart field.
    bed_df.sort_values(['seqname', 'start'], ascending=[True, True], inplace=True)

    # remove unused fields, and get order
    bed_df = bed_df[ORF_FIELDS]

    # Writes bed file to output directory
    # TODO: add option to write uncompressed w/o header 
    output_filename = '{}.bed.gz'.format(output_filename)
    output_filename = os.path.join(args.dirloc, output_filename)
    if os.path.exists(output_filename) and not args.overwrite:
        msg = "Output file {} already exists. Skipping.".format(output_filename)
        logger.warning(msg)
    else:
        pandas_utils.write_df(bed_df,
                              output_filename,
                              index=False,
                              sep='\t',
                              header=True,
                              do_not_compress=False,
                              quoting=csv.QUOTE_NONE)
        
    return


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Prepare one or more ORF BED files e.g.
        to be used with trackhub-utils. Standard 'rpbp/pbio' defaults are used.""")

    parser.add_argument('dirloc', help="""The output directory. All ORF BED files are 
        written to this location.""")

    parser.add_argument('--config', help="""The (yaml) config file for the Rp-Bp predictions.
        If not given, then [--input-list] with [--input-type f] must be used.""", type=str)

    parser.add_argument('--input-list', help="""A space-delimited list of input files, sample 
        names or conditions, each quoted separately. They must be either all files, or 
        sample/condition names (in which case the [--config] option must also be used), and
        this must be specified with the [--input-type]. Only these will be processed.""",
        nargs='*', required='--input-type' in sys.argv, type=str)

    parser.add_argument('--output-list', help="""A space-delimited list of output file base names,
        each quoted separately, WITHOUT extension (required if [--input-type f], else ignored.)""",
        nargs='*', type=str)

    parser.add_argument('--input-type', help="""The 'type' of [--input-list], either f (files),
        s (samples) or c (conditions).""", required='--input-list' in sys.argv, type=str,
        choices=['f', 's', 'c'])

    parser.add_argument('--no-color', help="""If this flag is present then custom color field 
        is NOT added for the ORF categories.""", action='store_true')

    parser.add_argument('-k', '--keep-other', help="""If this flag is present then ORFs labeled
        as "Other" will be included. They are discarded by default.""", action='store_true')

    parser.add_argument('-fid', '--filter-by-id', help="""Full path to a list of ORF ids, one per 
        line without header, to keep in the final set. This can also be a BED12 file, with the
        "id" field.""", type=str)

    parser.add_argument('-fps', '--filter-by-psites', help="""If this flag is present then keep ORFs
        with at least [--psite-threshold] in-frame p-sites.""", action='store_true')
    
    parser.add_argument('--psite-threshold', help="""Filter ORF predictions: keep ORFs
        with at least [--psite-threshold] in-frame p-sites. Silently ignored if [--filter-by-psites]
        is not set.""", type=int, default=10)

    parser.add_argument('-a', '--all-replicates', help="""If this flag is present then BED files
        are created for all replicates in the config file, in addition to the merged replicates or
        conditions. By default, only the latter are created, unless the option [no-merged] is set.""", 
        action='store_true')

    parser.add_argument('--no-merged', help="""If this flag is present then predictions from merged
        replicates are ignored.""", action='store_true')
    
    parser.add_argument('--overwrite', help='''If this flag is present, then existing files
        will be overwritten.''', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    if args.input_list and not args.config and args.input_type in ['s', 'c']:
        logger.critical("Missing [--config]")
        return
    if args.input_list and not args.output_list and args.input_type == 'f':
        logger.critical("Missing [--output-list]")
        return
    if args.no_merged and not args.all_replicates:
        logger.warning("""Option [--no-merged] is set without [--all-replicates], setting
            [--all-replicates] to True.""")
        args.all_replicates = True

    msg = "[prep-orf-beds]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    if args.config:
        config = yaml.load(open(args.config), Loader=yaml.FullLoader)
        required_keys = [
            'riboseq_data',
            'riboseq_samples',
            'riboseq_biological_replicates'
        ]
        utils.check_keys_exist(config, required_keys)

        # TODO add check for options when not all orf fields are used
        sample_name_map = ribo_utils.get_sample_name_map(config)
        condition_name_map = ribo_utils.get_condition_name_map(config)
        
    # check if we filter the list of ORFs
    args.id_list = None
    if args.filter_by_id:
        # check if bed or text file
        id_list = bed_utils.read_bed(args.filter_by_id)
        if len(id_list.columns) > 1:
            try:
                id_list = id_list.id.unique()
            except:
                msg = 'Using [--filter-by-id], but BED12+ file does not contain id field!'
                logger.critical(msg)
        else:
            id_list = set(open(args.filter_by_id).read().split())
        args.id_list = id_list
        
    files_only = False
    sample_names = {}
    condition_names = {}
    if args.input_list:
        if args.input_type == 's':
            logger.warning("""Using --input-type s, setting [--all-replicates] to True, and 
                ignoring merged replicates.""")
            args.all_replicates = True
            args.no_merged = True
            sample_names = {name: [name] for name in args.input_list}
        elif args.input_type == 'c':
            logger.warning("""Using --input-type c, setting [--no-merged] to False, and 
                ignoring replicates.""")
            args.all_replicates = False
            args.no_merged = False
            condition_names = {name: [name] for name in args.input_list}
        else:
            logger.warning("""Using --input-type f, ignoring replicate options.""")
            files_only = True
    else:
        sample_names = config['riboseq_samples']
        condition_names = ribo_utils.get_riboseq_replicates(config)

    if not files_only and args.output_list:
        msg = "[--output-list] will be ignored"
        logger.warning(msg)

    # check output path
    if os.path.exists(args.dirloc):
        args.dirloc = os.path.join(args.dirloc, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.dirloc)
        raise OSError(msg)

    if files_only:
        for bed_file, output_file in zip(args.input_list, args.output_list):
            if not os.path.exists(bed_file):
                msg = "Could not find the bed file: {}. Terminating.".format(bed_file)
                raise FileNotFoundError(msg)
            _get_bed(bed_file, output_file, args)


        return

    note_str = config.get('note', None)
    is_unique = not ('keep_riboseq_multimappers' in config)
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    if args.all_replicates:
        logger.info("Processing replicates.")
        for name in sorted(sample_names.keys()):
            msg = "Processing sample: {}".format(name)
            logger.info(msg)

            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config,
                                                                           name,
                                                                           is_unique=is_unique,
                                                                           default_params=metagene_options)

            predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'],
                                                                  name,
                                                                  length=lengths,
                                                                  offset=offsets,
                                                                  is_unique=is_unique,
                                                                  note=note_str,
                                                                  fraction=fraction,
                                                                  reweighting_iterations=reweighting_iterations,
                                                                  is_filtered=True)

            if not os.path.exists(predicted_orfs):
                msg = "Name: {}, file: {}".format(name, predicted_orfs)
                raise FileNotFoundError(msg)

            pretty_name = sample_name_map[name]
            _get_bed(predicted_orfs, pretty_name, args)

    # the merged replicates or conditions are always created, unless no_merged is set
    if args.no_merged:
        return
    
    logger.info("Processing merged replicates.")
    lengths = None
    offsets = None
    for name in sorted(condition_names.keys()):
        msg = "Processing condition: {}".format(name)
        logger.info(msg)

        predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'],
                                                              name,
                                                              length=lengths,
                                                              offset=offsets,
                                                              is_unique=is_unique,
                                                              note=note_str,
                                                              fraction=fraction,
                                                              reweighting_iterations=reweighting_iterations,
                                                              is_filtered=True)

        if not os.path.exists(predicted_orfs):
            msg = "Name: {}, file: {}".format(name, predicted_orfs)
            raise FileNotFoundError(msg)

        pretty_name = condition_name_map[name]
        _get_bed(predicted_orfs, pretty_name, args)


if __name__ == '__main__':
    main()
