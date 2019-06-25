#! /usr/bin/env python3

"""Collect all read length distribution .csv.gz files for all replicates
into a large table (unique counts only). If the files are not already available,
then they are created.

Note* Does not consider the filtered Ribo-seq BAM files with
      periodic fragment lengths (i.e. uses all fragments,
      they can easily be filtered afterwards); however, for RNA-seq, if
      mapping was done on trimmed reads, additional arguments must be
      given to find the right files.

Functions:
    add_data
"""

import os
import argparse
import logging
import yaml
import csv
import re

import pandas as pd

import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.pandas_utils as pandas_utils
import pbio.misc.utils as utils
import pbio.misc.shell_utils as shell_utils

import pbio.ribo.ribo_utils as ribo_utils
import pbio.ribo.ribo_filenames as filenames

from rpbp.defaults import metagene_options

logger = logging.getLogger(__name__)

data_get_length_distribution = {
    'rna': filenames.get_rnaseq_read_length_distribution,
    'ribo': filenames.get_riboseq_read_length_distribution,
}

data_get_bam = {
    'rna': filenames.get_rnaseq_bam,
    'ribo': filenames.get_riboseq_bam,
}


def add_data(sample_name, get_length_distribution, data, note, args):

    read_length_distribution_file = get_length_distribution(data,
                                                            sample_name,
                                                            note=note)

    read_length_distribution = pd.read_csv(read_length_distribution_file)

    ret = pd.pivot_table(read_length_distribution,
                         values='count',
                         columns='length',
                         index='basename')

    return ret


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', help="The yaml configuration file.")

    parser.add_argument('seq', choices=['rna', 'ribo'])

    parser.add_argument('out', help='''The output file complete path without extension. 
        If relevant, two files are created, one for all reads, and another for unique mappers.''')

    parser.add_argument('--overwrite', help='''Overwrites output file. This will NOT 
        overwrite the read length distribution files.''', action='store_true')

    parser.add_argument('--ribo-config', help="""Optional argument: the Ribo config file
        when seq is rna and rna reads have been trimmed to max ribo fragment lengths.
        If reads are trimmed then this needs to be given, otherwise the program will not find 
        the alignment files. In addition, the rna config file must include 'matching_samples'.""",
                        type=str)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    note = config.get('note', None)

    if args.ribo_config:
        ribo_config = yaml.load(open(args.ribo_config), Loader=yaml.FullLoader)
        is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)

    keep_key = 'keep_' + str(args.seq) + 'seq_multimappers'
    is_unique = not (keep_key in config)

    data_key = str(args.seq) + 'seq_data'
    rep_key = str(args.seq) + 'seq_samples'
    replicates = config[rep_key].keys()

    get_length_distribution = data_get_length_distribution[args.seq]
    get_bam = data_get_bam[args.seq]

    sample_name_map = ribo_utils.get_sample_name_map(config)

    # first check if the files exist, if not generate the read length distributions
    for replicate in replicates:

        read_length_distribution = get_length_distribution(config[data_key],
                                                           replicate,
                                                           note=note)

        ret = utils.check_files_exist([read_length_distribution], raise_on_error=False)
        if ret:
            continue

        msg = 'Creating {}'.format(read_length_distribution)
        logger.info(msg)

        if args.seq == 'rna' and args.ribo_config:
            config_keys = ['matching_samples']
            utils.check_keys_exist(config, config_keys)
            matching_ribo_sample = config['matching_samples'][replicate]

            # get the lengths, we don't need the offsets
            lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(ribo_config,
                                                                     matching_ribo_sample,
                                                                     is_unique=is_unique_ribo,
                                                                     default_params=metagene_options)

            if len(lengths) == 0:
                msg = "No periodic read lengths and offsets were found!"
                logger.error(msg)

            lengths = str(max([int(l) for l in lengths]))

        else:
            lengths = None

        #  all aligned reads
        genome_bam = get_bam(config[data_key],
                             replicate,
                             is_unique=False,
                             length=lengths,
                             note=note)

        if is_unique:
            # uniquely aligned reads
            unique_bam = get_bam(config[data_key],
                                 replicate,
                                 is_unique=is_unique,
                                 length=lengths,
                                 note=note)

            in_files = [genome_bam, unique_bam]
            cmd = "get-read-length-distribution {} {} --out {} {}".format(
                genome_bam,
                unique_bam,
                read_length_distribution,
                logging_str
            )
        else:
            in_files = [genome_bam]
            cmd = "get-read-length-distribution {} --out {} {}".format(
                genome_bam,
                read_length_distribution,
                logging_str
            )

        out_files = [read_length_distribution]
        shell_utils.call_if_not_exists(cmd,
                                       out_files,
                                       in_files=in_files,
                                       call=True)

    msg = 'Parsing read length distribution files'
    logger.info(msg)

    all_length_distributions = parallel.apply_iter_simple(
        replicates,
        add_data,
        get_length_distribution,
        config[data_key],
        note,
        args
    )

    all_length_distributions_df = pd.concat(all_length_distributions).fillna(0)
    all_length_distributions_df['name'] = all_length_distributions_df.index

    # adjust the names
    if args.seq == 'rna' and args.ribo_config:
        repl = re.compile('.length-\d{2}')
        all_length_distributions_df['name'] = all_length_distributions_df['name'].str.replace(repl, '')
    if note is not None:
        repl = '.{}'.format(note)
        all_length_distributions_df['name'] = all_length_distributions_df['name'].str.replace(repl, '')

    # split into all and unique
    m_unique = False
    out_files = [args.out + '.csv.gz']
    out_df = []
    if is_unique:
        m_unique = all_length_distributions_df['name'].str.contains('unique')
        unique_read_length_distribution_df = all_length_distributions_df[m_unique].copy()
        unique_read_length_distribution_df['name'] = unique_read_length_distribution_df['name'].str.replace('-unique',
                                                                                                            '')
        unique_read_length_distribution_df['name'] = unique_read_length_distribution_df['name'].apply(
            lambda x: sample_name_map[x])
        unique_read_length_distribution_df.set_index('name', inplace=True)
        unique_read_length_distribution_df.index.name = 'condition'
        out_files.append(args.out + '-unique.csv.gz')
        out_df.append(unique_read_length_distribution_df)

    all_length_distributions_df = all_length_distributions_df[~m_unique]
    all_length_distributions_df['name'] = all_length_distributions_df['name'].apply(lambda x: sample_name_map[x])
    all_length_distributions_df.set_index('name', inplace=True)
    all_length_distributions_df.index.name = 'condition'
    out_df.append(all_length_distributions_df)
    out_df.reverse()

    msg = "Writing output to: {}".format(','.join(out_files))
    logger.info(msg)

    for file, df in zip(out_files, out_df):
        if os.path.exists(file) and not args.overwrite:
            msg = "Output file {} already exists. Skipping.".format(file)
            logger.warning(msg)
        else:
            pandas_utils.write_df(df, file, create_path=True,
                                  index=True, sep=',', header=True,
                                  do_not_compress=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
