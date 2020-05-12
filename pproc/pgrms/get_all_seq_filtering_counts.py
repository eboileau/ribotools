#! /usr/bin/env python3

""" To replace eventually 'get-all-read-filtering-counts' in
Rp-Bp pipeline, adjusted to count reads from the final genome
alignment files (already filtered for periodic lengths and presumably
for library strandedness), and also include option to count
reads for the RNA pipeline. To be completed.
"""

import sys
import argparse
import yaml
import logging
import pandas as pd
import numpy as np

import pbio.ribo.ribo_filenames as ribo_filenames
import pbio.ribo.ribo_utils as ribo_utils

import pbio.utils.bam_utils as bam_utils
import pbio.utils.fastx_utils as fastx_utils

import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.shell_utils as shell_utils
import pbio.misc.pandas_utils as pandas_utils

logger = logging.getLogger(__name__)

default_num_cpus = 2


def get_counts(name_data, seq_key, sample_name_map, config, ribo_config, args):

    name, data = name_data
    msg = "processing {}...".format(name)
    logger.info(msg)

    data_key = seq_key + 'seq_data'
    keep_key = 'keep_' + seq_key + 'seq_multimappers'
    is_unique = not (keep_key in config)
    note = config.get('note', None)

    # get the periodic fragment lengths if ribo, or else if rna is trimmed
    trim_length = None
    lengths = None
    if seq_key == 'ribo':
        # get the lengths, we don't need the offsets
        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            name,
            is_unique=is_unique,
            isoform_strategy=args.isoform_strategy
        )

        if len(lengths) == 0:
            msg = "No periodic read lengths and offsets were found!"
            logger.critical(msg)
            return
    elif args.trimmed_rna:
        is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)
        matching_ribo_sample = config['matching_samples'][name]
        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            ribo_config,
            matching_ribo_sample,
            is_unique=is_unique_ribo,
            isoform_strategy=args.isoform_strategy
        )

        if len(lengths) == 0:
            msg = "No periodic read lengths and offsets were found!"
            logger.critical(msg)
            return

        lengths = str(max([int(l) for l in lengths]))
        trim_length = lengths

    # first, get the file names
    raw_data = data
    without_adapters = ribo_filenames.get_without_adapters_fastq(
        config[data_key], name, note=note)
    with_rrna = ribo_filenames.get_with_rrna_fastq(
        config[data_key], name, note=note)
    without_rrna = ribo_filenames.get_without_rrna_fastq(
        config[data_key], name, note=note)

    genome_bam = ribo_filenames.get_seq_bam(
        seq_key, config[data_key], name, isoform_strategy=args.isoform_strategy,
        length=trim_length, note=note)
    unique_bam = ribo_filenames.get_seq_bam(
        seq_key, config[data_key], name, isoform_strategy=args.isoform_strategy,
        length=trim_length, is_unique=is_unique, note=note)

    # now count the reads of each type
    msg = "{}: collecting read counts".format(name)
    logger.info(msg)
    
    # get the read counts
    msg = "{}: counting reads in raw data".format(name)
    logger.info(msg)
    raw_data_count = fastx_utils.get_read_count(raw_data, is_fasta=False)

    msg = "{}: counting reads without adapters".format(name)
    logger.info(msg)
    without_adapters_count = fastx_utils.get_read_count(without_adapters, is_fasta=False)
    
    msg = "{}: counting reads with rrna".format(name)
    logger.info(msg)
    with_rrna_count = fastx_utils.get_read_count(with_rrna, is_fasta=False)
    
    msg = "{}: counting reads without rrna".format(name)
    logger.info(msg)
    without_rrna_count = fastx_utils.get_read_count(without_rrna, is_fasta=False)
    
    msg = "{}: counting genome-aligned reads".format(name)
    logger.info(msg)
    genome_count = bam_utils.count_aligned_reads(genome_bam)

    msg = "{}: counting uniquely-aligned reads".format(name)
    logger.info(msg)
    unique_count = bam_utils.count_aligned_reads(unique_bam)

    # this is where the pipelines differ, whether we are counting reads from the final
    # alignment files (transcriptome-filtered genome alignments, adjusted for library
    # strandedness and periodicity if ribo), or else if we just count the periodic reads
    # from the uniquely-aligned reads for ribo, or if this is complete for rna.

    ret = {
        'note': sample_name_map[name],
        'raw_data_count': raw_data_count,
        'without_adapters_count': without_adapters_count,
        'without_rrna_count': without_rrna_count,
        'genome_count': genome_count,
        'unique_count': unique_count
    }

    if False:#args.read_count_pipeline:
        pass
        ## get transcriptome-filtered genome alignment file
        #filtered_genome_bam = ribo_filenames.get_seq_bam(
            #seq_key, config[data_key], name, isoform_strategy=args.isoform_strategy,
            #length=lengths, is_unique=is_unique, stranded=args.strandedness, note=note)
        #msg = "{}: counting transcriptome-filtered genome aligned reads".format(name)
        #logger.info(msg)
        #length_count = bam_utils.count_aligned_reads(filtered_genome_bam)
        #ret['length_count'] = length_count
    else:
        if seq_key == 'ribo':
            lengths_str = ','.join(lengths)
            length_counts = bam_utils.get_length_distribution(unique_bam)
            lengths = set([int(l) for l in lengths])
            m_lengths = length_counts['length'].isin(lengths)
            length_count = np.sum(length_counts.loc[m_lengths, 'count'])
            ret['length_count'] = length_count

            msg = ("{}: found the following periodic lengths: {}. The number of reads "
                   "of these lengths: {}".format(name, lengths_str, length_count))
            logger.info(msg)

    return pd.Series(ret)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script collects counts of ribo or rna seq reads filtered at "
        each step in the Rp-Bp or counting pipeline. ** to be completed ** """)

    parser.add_argument('config', help="The yaml config file.") # same ribo and rna

    parser.add_argument('out', help="The output csv file with the counts")
    parser.add_argument('-p', '--num-cpus', help="The number of processors to use", 
        type=int, default=default_num_cpus)
    parser.add_argument('--overwrite', action='store_true')

    parser.add_argument('--rna', help="""Optional argument: ribo-seq by default, unless "
        this flag is given, in which case reads will be collected for rna-seq. The correct
        configuration file must be used.""", action='store_true')

    parser.add_argument('--use-pretty-names', help="If this flag is given, then will use the names"
        "in 'ribo/rnaseq_sample_name_map' if present.", action='store_true')

    #subparser = parser.add_subparsers(help="""Optional arguments if collecting reads from 
        #the counting pipeline.""", dest='read_count_pipeline')
    #parser_read_count_pipeline = subparser.add_parser('read-count-pipeline', help=""" If this
        #flag is given, then additional options can be passed to collect the reads from
        #the counting pipeline. In particular, the 'usable reads' are taken from the final
        #genome alignments (see 'run-seq-align-pipeline').""")
    #parser_read_count_pipeline.add_argument('--strandedness', help="Library strandedness.",
        #choices=['fr', 'rf'])
    #parser_read_count_pipeline.add_argument('--isoform-strategy', help="""See b-tea.cl_utils.
        #Currently we only use 'merged'""", default=None)
    #parser_read_count_pipeline.add_argument('--trimmed-rna', help="""If rna reads are trimmed
        #to max fragment size from corresponding Ribo-seq data. N.B. At least the 
        #"periodic-offsets" file must be available. The config file must include 
        #"matching_samples" and the path to the Ribo-seq config must also be given using
         #[--ribo-config]). This is only relevant with [--rna].""",
        #action='store_true')
    #parser_read_count_pipeline.add_argument('--ribo-config', help="""Optional argument: the Ribo config file
            #if using [--trimmed-rna].""", required='--trimmed-rna' in sys.argv, type=str)


    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    programs = ['samtools']
    shell_utils.check_programs_exist(programs)

    #if args.read_count_pipeline and not args.strandedness:
        #msg = """[--strandedness] option is missing, falling back to default pipeline.
            #Other options passed with [read-count-pipeline] will still be parsed."""
        #logger.info(msg)
        #args.read_count_pipeline = None
    
    args.isoform_strategy = None

    config = yaml.load(open(args.config))

    ribo_config = None
    seq_key = 'ribo'
    if args.rna:
        seq_key = 'rna'
        #if args.trimmed_rna:
        args.trimmed_rna = True
        #    ribo_config = yaml.load(open(args.ribo_config))
        ribo_config = config
    sample_key = seq_key + 'seq_samples'

    if args.use_pretty_names:
        sample_name_map = ribo_utils.get_sample_name_map(config)
    else:
        sample_name_map = {name: [name] for name in config[sample_key].keys()}

    res = parallel.apply_parallel_iter(config[sample_key].items(),
                                       args.num_cpus,
                                       get_counts,
                                       seq_key, sample_name_map, config, ribo_config, args)

    res = [r for r in res if r is not None]
    res_df = pd.DataFrame(res)

    pandas_utils.write_df(res_df, args.out, index=False)
    
if __name__ == '__main__':
    main()

