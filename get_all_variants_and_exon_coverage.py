#! /usr/bin/env python3

"""Post-processing pipeline for Rp-Bp ORF predictions: group ORF
predictions (conditions only) by genomic loci. The latter are determined
based on exon overlap only (using annotations, annotated and/or de novo
in BED format). Additionally, add exon coverage using the ORF profiles.
** The merging and part of the processing could eventually be split into
separate scripts/programs.
"""

import os
import logging
import argparse
import yaml
import csv

import itertools
import scipy.io
import pandas as pd
import numpy as np
from collections import defaultdict

import bio_utils.bed_utils as bed_utils

import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.utils as utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as filenames

logger = logging.getLogger(__name__)

# using standard names, however we still need the Rp-Bp-related nomenclature...
std_names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
             'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
names = bed_utils.bed12_field_names + ['orf_num', 'orf_len', 'orf_type',
                                       'bayes_factor_mean', 'bayes_factor_var']
names_to_std_names = {name: std_name for name, std_name in zip(bed_utils.bed12_field_names, std_names)}
std_names_to_names = {std_name: name for name, std_name in names_to_std_names.items()}

loc_format = 'TLOC_{0:0>{1}}'


def _get_intervals(row, intron=False):
    """Pull start and end genomic coordinates
    out of the BED block fields, either exons
    or introns.

    Args:
        row (named tuple): BED entry default fields
    Returns:
        starts, ends: numpy array of matching starts, ends
    """

    chromStart = row.chromStart
    blockSizes = np.array(row.blockSizes.split(","), dtype=int)
    blockStarts = np.array(row.blockStarts.split(","), dtype=int)
    starts = chromStart + blockStarts
    ends = starts + blockSizes
    if intron:
        exon_starts = starts
        starts = ends[:-1]
        ends = exon_starts[1:]

    return starts, ends


def _get_all_intervals(df, **kwargs):
    """Get all intervals in dataframe as
    single lists, using name as info."""

    all_starts, all_ends, all_info = [], [], []
    for row in df.itertuples():
        row_starts, row_ends = _get_intervals(row, **kwargs)
        all_starts.append(row_starts)
        all_ends.append(row_ends)
        all_info.append([row.name] * len(row_starts))
    all_starts = np.concatenate(all_starts)
    all_ends = np.concatenate(all_ends)
    all_info = [info for row_info in all_info for info in row_info]

    return all_starts, all_ends, all_info


def _merge_intervals(interval_starts, interval_ends, interval_info=None):
    """Merge a list of intervals based on overlaps. This function is
    almost exactly the same as bed_utils.merge_intervals, except that
    the info is always added (when the cache is not empty), otherwise we
    don't know which intervals are being merged together. The end interval
    is also included when extending, but this should not make a difference.
    See the documentation fo this function."""

    if len(interval_starts) == 0:
        return [], [], []

    interval_starts = np.array(interval_starts)
    interval_ends = np.array(interval_ends)

    sorted_interval_start_indices = np.argsort(interval_starts)

    interval_starts = interval_starts[sorted_interval_start_indices]
    interval_ends = interval_ends[sorted_interval_start_indices]

    if interval_info is not None:
        interval_info = np.array(interval_info)
        interval_info = interval_info[sorted_interval_start_indices]
        next_interval_info = interval_info[0]

    else:
        next_interval_info = None

    next_interval = 0
    num_intervals = len(interval_starts)
    next_interval_start = interval_starts[0]
    next_interval_end = interval_ends[0]

    cache = []

    merged_starts = []
    merged_ends = []
    merged_info = []

    cur_interval_start = None
    cur_interval_end = None
    cur_interval_info = None

    while next_interval < num_intervals:
        # first, remove everything from the cache that ended before this begins
        cache = [c for c in cache if interval_ends[c] > next_interval_start]

        # if we cleared the cache, then we just finished an interval
        if len(cache) == 0:
            if cur_interval_start is not None:
                # add it to the output list
                merged_starts.append(cur_interval_start)
                merged_ends.append(cur_interval_end)
                merged_info.append(cur_interval_info)

            # and start a new interval
            cur_interval_start = next_interval_start
            cur_interval_end = next_interval_end
            cur_interval_info = [next_interval_info]

        else:
            # otherwise, extend the previous interval
            if next_interval_end >= cur_interval_end:
                cur_interval_end = next_interval_end
            # always add info, unless cache is empty
            cur_interval_info.append(next_interval_info)

        # add the next interval to the current cache
        cache.append(next_interval)

        # and advance
        next_interval += 1
        if next_interval < num_intervals:

            next_interval_start = interval_starts[next_interval]
            next_interval_end = interval_ends[next_interval]

            if interval_info is not None:
                next_interval_info = interval_info[next_interval]

    # add the last merged interval to the output list
    merged_starts.append(cur_interval_start)
    merged_ends.append(cur_interval_end)
    merged_info.append(cur_interval_info)

    merged_starts = np.array(merged_starts)
    merged_ends = np.array(merged_ends)

    return merged_starts, merged_ends, merged_info


def _collapse_all(locus):
    """Collapse all transcripts under a given locus, keeping
    the largest exonic span, i.e. there is no container and/or
    containment restrictions. Coordinates only reveal the
    location (trace or imprint) of all merged exons."""

    # get largest container, but only for assigning a name
    locus['span'] = locus['chromEnd'] - locus['chromStart']
    locus.sort_values('span', ascending=False, inplace=True)
    name = locus['name'].iloc[0]

    # get all exon starts and ends
    starts, ends = [], []
    collapsed = []
    for row in locus.itertuples():
        collapsed.append(row.name)
        row_starts, row_ends = _get_intervals(row)
        starts.append(row_starts)
        ends.append(row_ends)
    collapsed.remove(name)

    starts = np.concatenate(starts)
    ends = np.concatenate(ends)

    merged_intervals = bed_utils.merge_intervals(starts, ends)
    merged_starts, merged_ends, _ = merged_intervals

    # convert back to BED
    chromStart = np.min(merged_starts)
    chromEnd = np.max(merged_ends)

    blockCount = len(merged_starts)
    blockSizes = ','.join(str(size) for size in merged_ends - merged_starts)
    blockStarts = ','.join(str(start) for start in merged_starts - chromStart)
    dupinfo = ','.join(str(n) for n in collapsed)

    collapsed_locus = pd.DataFrame(index=[0])  # pandas version?!
    collapsed_locus['chrom'] = locus['chrom'].iloc[0]
    collapsed_locus['chromStart'] = chromStart
    collapsed_locus['chromEnd'] = chromEnd
    collapsed_locus['name'] = name
    collapsed_locus['score'] = '0'
    collapsed_locus['strand'] = locus['strand'].iloc[0]
    collapsed_locus['thickStart'] = '-1'
    collapsed_locus['thickEnd'] = '-1'
    collapsed_locus['itemRgb'] = '0'
    collapsed_locus['blockCount'] = blockCount
    collapsed_locus['blockSizes'] = blockSizes
    collapsed_locus['blockStarts'] = blockStarts
    collapsed_locus['locus'] = locus['locus'].iloc[0]
    collapsed_locus['dupInfo'] = dupinfo

    return collapsed_locus


def _count_exon_coverage_from_profiles(orfs, config):
    """Count periodic coverage for each exon block"""

    is_unique = not ('keep_riboseq_multimappers' in config)
    note_str = config.get('note', None)

    condition = orfs['source'].iloc[0]

    mtx = filenames.get_riboseq_profiles(
        config['riboseq_data'],
        condition,
        length=None,
        offset=None,
        is_unique=is_unique,
        note=note_str
    )

    mtx = scipy.io.mmread(mtx)

    exon_coverage = pd.DataFrame(columns=['orfBlockCoverage', 'source'])

    for row in orfs.itertuples():
        dense = utils.to_dense(mtx, row.orf_num, float, length=row.orf_len)
        if row.strand == '-':
            dense = dense[::-1]

        frames = []
        for in_frame in range(3):
            idx = [i for i in range(row.orf_len) if i % 3 != in_frame]
            frames.append(dense.copy())
            frames[in_frame][idx] = 0

        block_sizes = np.array(row.blockSizes.split(","), dtype=int)
        block_sizes = np.insert(block_sizes, 0, 0)
        exon_blocks_loc = block_sizes.cumsum()

        exon_block_frame_counts = []
        for exon_block_start, exon_block_end in zip(exon_blocks_loc[:-1], exon_blocks_loc[1:]):
            frame_counts = ':'.join(str(int(frame[exon_block_start:exon_block_end].sum())) for frame in frames)
            exon_block_frame_counts.append(frame_counts)

        exon_coverage.at[row.name, 'orfBlockCoverage'] = ';'.join(frame_counts for
                                                                  frame_counts in exon_block_frame_counts)
        exon_coverage['source'] = condition

    return exon_coverage


def _map_orfs_to_locus(loci, orfs):
    """Merge all predictions at a given locus and
    redefine coordinates relative to locus. These results
    in a BED12+ format with additional columns for the
    predictions."""

    locus = loci['locus'].item()
    m_loc = orfs['locus'] == locus
    orfs_loc = orfs[m_loc]

    orf_ids = []
    orf_types = []
    orf_sources = []
    orf_block_sizes = []
    orf_block_counts = []
    orf_block_starts = []
    orf_bf_means = []
    orf_bf_vars = []
    orf_block_coverage = []

    for row in orfs_loc.itertuples():
        row_starts, row_ends = _get_intervals(row)
        row_block_starts = ','.join(str(start) for start in row_starts - loci['chromStart'].item())
        orf_ids.append(row.name)
        orf_types.append(row.orf_type)
        orf_sources.append(row.source)
        orf_block_sizes.append(row.blockSizes)
        orf_block_counts.append(row.blockCount)
        orf_block_starts.append(row_block_starts)
        orf_block_coverage.append(row.orfBlockCoverage)
        orf_bf_means.append(row.bayes_factor_mean)
        orf_bf_vars.append(row.bayes_factor_var)

    # loc_orfs = loci.drop('locus', axis=1)
    loc_orfs = loci[std_names].copy()

    loc_orfs['orfs'] = ','.join(str(i) for i in orf_ids)
    loc_orfs['types'] = ','.join(str(i) for i in orf_types)
    loc_orfs['sources'] = ','.join(str(i) for i in orf_sources)
    loc_orfs['orfBlockCount'] = ','.join(str(i) for i in orf_block_counts)
    loc_orfs['orfBlockSizes'] = ';'.join(str(i) for i in orf_block_sizes)
    loc_orfs['orfBlockStarts'] = ';'.join(str(i) for i in orf_block_starts)
    loc_orfs['orfBlockCoverage'] = ','.join(str(i) for i in orf_block_coverage)
    loc_orfs['bfMeans'] = ','.join(str(i) for i in orf_bf_means)
    loc_orfs['bfVars'] = ','.join(str(i) for i in orf_bf_vars)

    return loc_orfs


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Group ORF predictions by genomic loci and find exon coverage using
        the ORF profiles. At least the ORF profiles and predictions for the conditions and
        the annotations used must be available. """)

    parser.add_argument('config', help="The (yaml) config file.")
    parser.add_argument('outname', help="The name for the final output "
        "(extended ORF predictions) wo extension.")
    parser.add_argument('outloc', help="The output directory.")

    parser.add_argument('--use-de-novo', help="If de novo annotations are available, "
        "these will be added to the standard annotations.", action='store_true')

    parser.add_argument('--input-list', help="""A space-delimited list of sample names, 
        each quoted separately. They must match the keys in the configuration file. 
        Only these will be converted.""", nargs='*', type=str)


    parser.add_argument('-p', '--num-cpus', help="The number of processors to use",
                        type=int, default=2)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))
    is_unique = not ('keep_riboseq_multimappers' in config)
    note_str = config.get('note', None)

    required_keys = [
        'riboseq_data',
        'riboseq_samples',
        'riboseq_biological_replicates',
        'genome_base_path',
        'genome_name'
    ]
    utils.check_keys_exist(config, required_keys)

    riboseq_replicates = ribo_utils.get_riboseq_replicates(config)


    # Part 1: Cluster annotations into genomic loci

    # ** Single exons falling entirely within intronic regions
    # of multiple exon transcripts are assigned their own locus.
    # In such cases, overlapping single exons are not merged; besides
    # we currently have no way to deal with same sense overlapping genes,
    # i.e. loci are assigned only based on the information in the BED12 file.

    msg = "Clustering annotations into genomic loci"
    logger.info(msg)

    is_de_novo = False
    beds = [filenames.get_bed(config['genome_base_path'], config['genome_name'],
                            is_merged=False, is_annotated=True, is_de_novo=is_de_novo)]
    if args.use_de_novo:
        is_de_novo = True
        beds.append(filenames.get_bed(config['genome_base_path'], config['genome_name'],
                            is_merged=False, is_annotated=False, is_de_novo=is_de_novo))
    bed = bed_utils.concatenate(beds, sort_bed=True)
    bed.reset_index(inplace=True, drop=True)
    bed.rename(columns=names_to_std_names, inplace=True)

    bed['locus'] = np.nan
    zfill = len(str(len(bed)))
    chroms = bed['chrom'].unique()
    strands = ("+", "-")

    counter = itertools.count(0)
    for strand in strands:
        m_strand = bed['strand'] == strand
        for chrom in chroms:
            m_chrom = bed['chrom'] == chrom
            m_filter = m_strand & m_chrom
            if sum(m_filter) == 0:
                continue

            loc_starts = bed.loc[m_filter, 'chromStart']
            loc_ends = bed.loc[m_filter, 'chromEnd']
            loc_info = bed.loc[m_filter, 'name']

            _, _, merged = _merge_intervals(loc_starts, loc_ends, loc_info)

            loci = []
            for loc_tmp in merged:
                loc = bed[bed['name'].isin(loc_tmp)]
                m_single_block = loc['blockCount'] == 1
                single_loc = []
                if m_single_block.any():
                    single_exon_blocks = loc[m_single_block]
                    multi_exon_blocks = loc[~m_single_block]
                    if multi_exon_blocks.empty:
                        loci.append(loc_tmp)
                        continue

                    multi_starts, multi_ends, multi_loc = _get_all_intervals(multi_exon_blocks)
                    single_starts, single_ends, single_loc = _get_all_intervals(single_exon_blocks)

                    overlap = bed_utils.get_interval_overlaps(multi_starts, multi_ends, multi_loc,
                                                              single_starts, single_ends, single_loc)

                    # split info: if no overlap, all single exons are taken apart, else make sure
                    # we only keep the right cluster
                    overlap_loc = list(set([loc_overlap.b_info for loc_overlap in overlap]))
                    single_loc = [loc_single for loc_single in single_loc if loc_single not in overlap_loc]

                merged_loc = [loc_merged for loc_merged in loc_tmp if loc_merged not in single_loc]
                if merged_loc:
                    loci.append(merged_loc)
                for loc_single in single_loc:
                    loci.append([loc_single])

            id_map = {name: loc_format.format(locid, zfill) for locid, locus in zip(counter, loci)
                      for name in locus}
            locids = bed['name'].map(id_map)
            locids_df = pd.DataFrame({'locus': locids})
            bed.update(locids_df)

    # then actually collapse all loci

    loci = bed.groupby('locus')
    collapsed_loci = parallel.apply_parallel_groups(
            loci,
            args.num_cpus,
            _collapse_all
            )
    collapsed_loci_df = pd.concat(collapsed_loci)

    # adjust the fields, and write to disk, also write collapsed transcripts
    dupinfo = collapsed_loci_df[['name', 'dupInfo']].copy()
    dupinfo = dupinfo[~(dupinfo['dupInfo']=='')]
    dupinfo.to_csv(os.path.join(args.outloc, 'dupinfo.txt'), sep=' ',
                   quoting=csv.QUOTE_NONE, index=None, header=False)

    cols = std_names + ['locus']
    bed_trace = filenames.get_bed(args.outloc, config['genome_name'],
                      is_merged=True, is_annotated=True, is_de_novo=is_de_novo)
    collapsed_loci_df.to_csv(bed_trace, sep='\t', quoting=csv.QUOTE_NONE, compression='gzip',
                             index=None, header=False, columns=cols)

    # Part 2: Merge all predictions (conditions only)

    msg = "Merging the ORF predictions."
    logger.info(msg)

    beds = []
    for condition in riboseq_replicates.keys():
        predicted_orfs = filenames.get_riboseq_predicted_orfs(
            config['riboseq_data'],
            condition,
            is_unique=is_unique,
            is_filtered=True,
            note=note_str)

        bed = bed_utils.get_bed_df(predicted_orfs)
        bed = bed[names].copy()
        bed['source'] = condition
        beds.append(bed)

    predicted_orfs_df = pd.concat(beds)
    predicted_orfs_df = bed_utils.sort(predicted_orfs_df)
    predicted_orfs_df.reset_index(inplace=True, drop=True)


    # Part 3: Now actually match predictions with loci (annotations) and group together

    msg = "Finding exon coverage."
    logger.info(msg)

    # For the annotations we need again the names... not the standard names!
    # For the ORFs, they have not yet been modified.
    collapsed_loci_df.rename(columns=std_names_to_names, inplace=True)

    # convert to bed6 and find overlap
    bed6_collapsed_loci = bed_utils.split_bed12(collapsed_loci_df, num_cpus=args.num_cpus)
    bed6_predicted_orfs = bed_utils.split_bed12(predicted_orfs_df, num_cpus=args.num_cpus)

    overlap = bed_utils.get_bed_overlaps(bed6_collapsed_loci, bed6_predicted_orfs)
    overlaps = defaultdict(list)
    for transcript_overlap in overlap:
        overlaps[transcript_overlap.a_info].append(transcript_overlap.b_info)
    # now reconstruct the loci using the predictions
    id_map = {name: collapsed_loci_df.loc[collapsed_loci_df['id'] == overlap, 'locus'].item() for overlap in
              overlaps.keys() for name in overlaps[overlap]}
    locids = predicted_orfs_df['id'].map(id_map)
    predicted_orfs_df['locus'] = locids

    # then one last time adjust the header...
    collapsed_loci_df.rename(columns=names_to_std_names, inplace=True)

    keep = names + ['locus', 'source']
    predicted_orfs_df = predicted_orfs_df[keep].copy()
    predicted_orfs_df.rename(columns=names_to_std_names, inplace=True)

    # find exon coverage, and add info to predictions
    sources = predicted_orfs_df.groupby('source')
    exon_coverage = parallel.apply_parallel_groups(
        sources,
        args.num_cpus,
        _count_exon_coverage_from_profiles,
        config
        )
    exon_coverage_df = pd.concat(exon_coverage)

    exon_coverage_df.reset_index(inplace=True)
    exon_coverage_df.rename(columns={'index': 'name'}, inplace=True)
    predicted_orfs_df = predicted_orfs_df.merge(exon_coverage_df, how='outer', on=['name', 'source'])


    # then actually collapse all loci

    msg = "Matching predictions with genomic loci and writing to disk."
    logger.info(msg)

    # grouping by loci is actually equivalent to grouping by every single row...
    # we should re-write this...
    loci = collapsed_loci_df.groupby('locus')
    predictions_loc = parallel.apply_parallel_groups(
            loci,
            args.num_cpus,
            _map_orfs_to_locus,
            predicted_orfs_df,
            )
    predictions_loc_df = pd.concat(predictions_loc)

    extended_predictions = os.path.join(args.outloc, '{}.bed.gz'.format(args.outname))
    bed_utils.write_bed(predictions_loc_df, extended_predictions)


if __name__ == '__main__':

    main()

