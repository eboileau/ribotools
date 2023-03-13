#! /usr/bin/env python3

"""Pull Rp-Bp ORF predictions into a unique set, assign gene/transcript ids
based on BED operations by replicating the label assignment (label_orfs.py).
This is not particularly efficient or clever, we just want to match predictions
based on the same heuristic as we do with the labels.

Required:
    APPRIS scores

Functions:
    get_transcript_id
"""

import sys
import os
import argparse
import logging
import pandas as pd
import numpy as np
import yaml


import pbio.utils.bed_utils as bed_utils

import pbio.ribo.ribo_filenames as filenames

import pbio.misc.shell_utils as shell_utils
import pbio.misc.logging_utils as logging_utils

from rpbp.defaults import default_num_cpus

logger = logging.getLogger(__name__)


# APPRIS fields
appris_names = ['gene_id', 'gene_name', 'transcript_id', 'protein_id', 'coding_label', 'transcript_biotype', 
                'start_stop', 'ccds', 'tsl', 'missing', 'firestar', 'matador', 'corsair', 'spade', 'thump', 
                'crash', 'inertia', 'cnio', 'appris_score', 'label']

# default filenames
all_orfs_file = 'filtered.predicted-orfs.bed.gz'
unique_orfs_file = 'unique-orfs.appris.bed.gz'

to_keep = bed_utils.bed12_field_names + ['orf_len', 'orf_type']



# Find a "main" transcript using APPRIS, then merge the remaining fields.
# The choice is arbitrary if more than one transcript satisfy the same conditions, and
# in any case should yield the same gene_id.
            
def get_transcript_id(row):
    
    ids = row.matching_transcript_ids
    transcripts = appris[appris['transcript_id'].str.contains(ids)]
    # it's not in APPRIS
    if transcripts.empty:
        return np.nan
    
    # make sure we start with empty dataframe, in case we have an exception
    principal = pd.DataFrame()
    # pick principal, if there are missing values this will raise an exception
    # if there are no principal, we fall back to ordering only
    try:
        principal = transcripts[transcripts.label.str.startswith('PRINCIPAL')]
        principal = principal.sort_values(by='appris_score', ascending=False)
    except:
        pass
    finally:
        # either no principal, or there were some missing values...
        # in any case, we still are "allowed" to sort the scores
        if principal.empty: 
            principal = transcripts.sort_values(by='appris_score', ascending=False)
    
    return principal.iloc[0].transcript_id
    
    
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''Assign gene/transcript ids and APRRIS scores
                                     to Rp-Bp ORF predictions, following label assignment. This script
                                     does NOT apply any filtering (merged replicates, unassigned 
                                     ids, etc.).''')

    parser.add_argument('config', help="The (yaml) configuration file.")

    parser.add_argument('dirloc', help="Path to the working directory for I/O files.")

    parser.add_argument('appris', help="The APPRIS scores file. Note: Fields are hard coded!")
    
    parser.add_argument('--orfs-file', help="""If given, will use this file as input.""", type=str)
    
    parser.add_argument('-n', '--name', help="""If given, basename for output, 'project_name'
                        from the config has precedence.""", type=str, default='rpbp')
    
    parser.add_argument('-p', '--num-cpus', help='''The number of CPUs to use to perform
            BED operations.''', type=int, default=default_num_cpus)
    

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[match-appris-scores]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    required_keys = [
        'genome_base_path',
        'genome_name',
    ]
    utils.check_keys_exist(config, required_keys)
    
    
    annotated_transcripts = filenames.get_bed(config['genome_base_path'],
                                              config['genome_name'],
                                              is_annotated=True)
    
    appris = pd.read_csv(args.appris), 
                         sep='\t', 
                         header=None, 
                         names=appris_names)
    
    basename = config.get("project_name", args.name)
    
    # output
    unique_orfs = '{}.{}'.format(basename, unique_orfs_file)
    unique_orfs = os.path.join(args.dirloc, unique_orfs)
        
    # find all ORF predictions, else create
    if not args.orfs_file:
        orfs_file = '{}.{}'.format(basename, all_orfs_file)
        args.orfs_file = os.path.join(args.dirloc, orfs_file)

    try:
        predicted_orfs = bed_utils.read_bed(args.orfs_file)
    except FileNotFoundError:
        msg = "File {} not found! Calling [get-all-predictions].".format(args.orfs_file)
        logger.info(msg)
        
        logging_str = logging_utils.get_logging_options_string(args)
        cmd = "get-all-predictions {} {} {}".format(args.config,
                                                    args.orfs_file,
                                                    logging_str)
        out_files = [args.orfs_file]
        shell_utils.call_if_not_exists(cmd, 
                                       out_files)
    
        predicted_orfs = bed_utils.read_bed(args.orfs_file)
        
    predicted_orfs = predicted_orfs[to_keep]
    # keep the unique ORFs only based on BED12+orf_len+orf_type
    predicted_orfs.drop_duplicates(inplace=True)
 
    msg = "Starting with {} (unique) predicted ORFs".format(len(predicted_orfs))
    logger.info(msg)
 
    # first extract CDSs from annotated transcripts...
    msg = "Removing the annotated UTRs from the transcripts"
    logger.info(msg)
    canonical_orfs = bed_utils.retain_all_thick_only(annotated_transcripts, 
                                                     num_cpus=args.num_cpus)

    msg = "Splitting the canonical ORFs into exons"
    logger.info(msg)
    canonical_orf_exons = bed_utils.split_bed12(canonical_orfs, 
                                                num_cpus=args.num_cpus, 
                                                progress_bar=False)

    # ... then extract 5'UTRs...
    msg = "Extracting annotated 5' leader regions from the transcripts"
    logger.info(msg)
    five_prime_regions = bed_utils.retain_all_five_prime_of_thick(annotated_transcripts, 
                                                                  num_cpus=args.num_cpus)
    if len(five_prime_regions) == 0:
        msg = "No annotated 5' leader regions were found"
        logger.warning(msg)

    msg = "Splitting the 5' leaders into exons"
    logger.info(msg)
    five_prime_exons = bed_utils.split_bed12(five_prime_regions, 
                                             num_cpus=args.num_cpus, 
                                             progress_bar=False)

    # ... 3'UTRs...
    msg = "Extracting annotated 3' trailer regions"
    logger.info(msg)
    three_prime_regions = bed_utils.retain_all_three_prime_of_thick(annotated_transcripts, 
                                                                    num_cpus=args.num_cpus)
    if len(three_prime_regions) == 0:
        msg = "No annotated 3' trailer regions were found"
        logger.warning(msg)

    msg = "Splitting the 3' trailers into exons"
    logger.info(msg)
    three_prime_exons = bed_utils.split_bed12(three_prime_regions,
                                              num_cpus=args.num_cpus,
                                              progress_bar=False)

    # ... and annotated noncoding 
    msg = "Splitting non-coding transcripts into exons"
    logger.info(msg)

    m_no_thick_start = annotated_transcripts['thick_start'] == -1
    m_no_thick_end = annotated_transcripts['thick_end'] == -1
    m_no_thick = m_no_thick_start & m_no_thick_end
    noncoding_transcripts = annotated_transcripts[m_no_thick]
    noncoding_exons = bed_utils.split_bed12(noncoding_transcripts,
                                            num_cpus=args.num_cpus,
                                            progress_bar=False)
    
    # Split all ORFs into exons, we only do this to check that we have no ORFs remaining
    all_orf_exons = bed_utils.split_bed12(predicted_orfs[bed_utils.bed12_field_names], 
                                          num_cpus=args.num_cpus, 
                                          progress_bar=False)
    
    
    # We know that the labels are supposedly assigned respecting the transcript structure, 
    # as we've done in "label_orfs.py", so we just want to find the "matching transcripts".
    # We use the same "ordering of assignment" as in "label_orfs.py".

    # canonical
    predicted_cds = predicted_orfs[predicted_orfs['orf_type']=='canonical'][bed_utils.bed12_field_names]
    msg = "Found {} predicted canonical".format(len(predicted_cds))
    logger.info(msg)

    msg = "Splitting these into exons"
    logger.info(msg)

    predicted_cds_exons = bed_utils.split_bed12(predicted_cds, 
                                                num_cpus=args.num_cpus, 
                                                progress_bar=False)

    msg = "Finding exact matches between predicted canonical ORFs and transcripts"
    logger.info(msg)

    exact_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                               predicted_cds_exons, 
                                               min_a_overlap=1, 
                                               min_b_overlap=1)

    msg = "Found {} matches".format(len(exact_matches))
    logger.info(msg)

    # Pull that into a dataframe and add APPRIS scores.
    exact_match_pairs = [(m.b_info, m.a_info) for m in exact_matches]
    MATCH_CDS = pd.DataFrame(exact_match_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_CDS = MATCH_CDS.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    match_ids = {m.b_info for m in exact_matches}
    m_match = all_orf_exons['id'].isin(match_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # canonical_variant
    predicted_variants = predicted_orfs[predicted_orfs['orf_type']=='canonical_variant'][bed_utils.bed12_field_names]
    msg = "Found {} predicted canonical_variant".format(len(predicted_variants))
    logger.info(msg)

    msg = "Splitting these into exons"
    logger.info(msg)

    predicted_variants_exons = bed_utils.split_bed12(predicted_variants, 
                                                     num_cpus=args.num_cpus, 
                                                     progress_bar=False)

    msg = "Finding matches between predicted variants (truncated or extended) and transcripts"
    logger.info(msg)

    truncated_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                   predicted_variants_exons,
                                                   min_b_overlap=1)

    truncated_match_ids = {m.b_info for m in truncated_matches}
    m_truncated_matches = predicted_variants_exons['id'].isin(truncated_match_ids)
    predicted_variants_exons = predicted_variants_exons[~m_truncated_matches]
        
    extended_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                  predicted_variants_exons,
                                                  min_a_overlap=1)

    extended_match_ids = {m.b_info for m in extended_matches}
    m_extended_matches = predicted_variants_exons['id'].isin(extended_match_ids)
    predicted_variants_exons = predicted_variants_exons[~m_extended_matches]

    msg = "Remaining {} canonical_variant exons".format(len(predicted_variants_exons))
    logger.info(msg)       

    truncated_matches.extend(extended_matches)

    msg = "Found {} matches in total".format(len(truncated_matches))
    logger.info(msg)

    # Pull that into a dataframe and add APPRIS scores
    variant_match_pairs = [(m.b_info, m.a_info) for m in truncated_matches]
    MATCH_VAR = pd.DataFrame(variant_match_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_VAR = MATCH_VAR.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    match_ids = {m.b_info for m in truncated_matches}
    m_match = all_orf_exons['id'].isin(match_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # within
    predicted_within = predicted_orfs[predicted_orfs['orf_type']=='within'][bed_utils.bed12_field_names]
    msg = "Found {} predicted within".format(len(predicted_within))
    logger.info(msg)

    msg = "Splitting these into exons"
    logger.info(msg)

    predicted_within_exons = bed_utils.split_bed12(predicted_within, 
                                                   num_cpus=args.num_cpus, 
                                                   progress_bar=False)

    msg = "Finding matches between predicted within and transcripts"
    logger.info(msg)

    within_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                predicted_within_exons,
                                                min_b_overlap=1)

    within_ids = {m.b_info for m in within_matches if m.b_info not in truncated_match_ids}
    m_within_matches = predicted_within_exons['id'].isin(within_ids)
    predicted_within_exons = predicted_within_exons[~m_within_matches]
    msg = "Remaining {} within exons".format(len(predicted_within_exons))
    logger.info(msg)  

    # Pull that into a dataframe and add APPRIS scores.
    within_pairs = {(m.b_info, m.a_info) for m in within_matches if m.b_info not in truncated_match_ids}
    MATCH_WITHIN = pd.DataFrame(within_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_WITHIN = MATCH_WITHIN.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    m_match = all_orf_exons['id'].isin(within_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # Find the overlap with predicted uORFs. We know that these are valid uORFs, i.e. they will be contained in 
    # the transcript structure, but we still need to make sure that they are associated with the right transcript!

    predicted_uorfs = predicted_orfs[predicted_orfs['orf_type'].isin(['five_prime', 'five_prime_overlap'])][bed_utils.bed12_field_names]
    msg = "Found {} predicted five_prime or five_prime_overlap".format(len(predicted_uorfs))
    logger.info(msg)
                                    
    msg = "Splitting these into exons"
    logger.info(msg)
                                    
    predicted_uorfs_exons = bed_utils.split_bed12(predicted_uorfs, 
                                                  num_cpus=args.num_cpus, 
                                                  progress_bar=False)

    out_of_frame_matches = bed_utils.get_bed_overlaps(canonical_orf_exons, 
                                                      predicted_uorfs_exons)
    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons, 
                                                predicted_uorfs_exons)
                                    
    leader_match_pairs = {(m.a_info, m.b_info) for m in leader_matches}
    leader_overlaps = {m for m in out_of_frame_matches if (m.a_info, m.b_info) in leader_match_pairs}
    overlap_ids = {m.b_info for m in leader_overlaps}
    m_overlap_matches = predicted_uorfs_exons['id'].isin(overlap_ids)
    predicted_uorfs_exons = predicted_uorfs_exons[~m_overlap_matches]
    msg = "Found {} overlap matches between five_prime_overlap and transcripts".format(len(overlap_ids))
    logger.info(msg)

    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons, 
                                                predicted_uorfs_exons, 
                                                min_b_overlap=1)
    msg = "Found {} independent matches between five_prime and transcripts".format(len(leader_matches))
    logger.info(msg)
            
    leader_ids = {m.b_info for m in leader_matches}
    m_leader_matches = predicted_uorfs_exons['id'].isin(leader_ids)
    predicted_uorfs_exons = predicted_uorfs_exons[~m_leader_matches]
                                    
    msg = "Remaining {} five_prime or five_prime_overlap exons".format(len(predicted_uorfs_exons))
    logger.info(msg)                       

    leader_matches.extend(leader_overlaps)

    msg = "Found {} matches in total".format(len(leader_matches))
    logger.info(msg)

    # Pull that into a dataframe and add APPRIS scores.
    leader_match_pairs = [(m.b_info, m.a_info) for m in leader_matches]
    MATCH_uORF = pd.DataFrame(leader_match_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_uORF = MATCH_uORF.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    match_ids = {m.b_info for m in leader_matches}
    m_match = all_orf_exons['id'].isin(match_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # Do the same with dORFS.
    predicted_dorfs = predicted_orfs[predicted_orfs['orf_type'].isin(['three_prime', 'three_prime_overlap'])][bed_utils.bed12_field_names]
    msg = "Found {} predicted three_prime or three_prime_overlap".format(len(predicted_dorfs))
    logger.info(msg)
                                    
    msg = "Splitting these into exons"
    logger.info(msg)
                                    
    predicted_dorfs_exons = bed_utils.split_bed12(predicted_dorfs, 
                                                  num_cpus=args.num_cpus, 
                                                  progress_bar=False)

    out_of_frame_matches = bed_utils.get_bed_overlaps(canonical_orf_exons, 
                                                      predicted_dorfs_exons)
    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons, 
                                                 predicted_dorfs_exons)
    trailer_match_pairs = {(m.a_info, m.b_info) for m in trailer_matches}
    trailer_overlaps = {m for m in out_of_frame_matches if (m.a_info, m.b_info) in trailer_match_pairs}
    overlap_ids = {m.b_info for m in trailer_overlaps}
    m_overlap_matches = predicted_dorfs_exons['id'].isin(overlap_ids)
    predicted_dorfs_exons = predicted_dorfs_exons[~m_overlap_matches]
    msg = "Found {} overlap matches between three_prime_overlap and transcripts".format(len(overlap_ids))
    logger.info(msg)

    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons, 
                                                 predicted_dorfs_exons, 
                                                 min_b_overlap=1)
    msg = "Found {} independent matches between three_prime and transcripts".format(len(trailer_matches))
    logger.info(msg)
            
    trailer_ids = {m.b_info for m in trailer_matches}
    m_trailer_matches = predicted_dorfs_exons['id'].isin(trailer_ids)
    predicted_dorfs_exons = predicted_dorfs_exons[~m_trailer_matches]

    msg = "Remaining {} three_prime or three_prime_overlap exons".format(len(predicted_dorfs_exons))
    logger.info(msg)  

    trailer_matches.extend(trailer_overlaps)

    msg = "Found {} matches in total".format(len(trailer_matches))
    logger.info(msg)

    # Pull that into a dataframe and add APPRIS scores.
    trailer_match_pairs = [(m.b_info, m.a_info) for m in trailer_matches]
    MATCH_dORF = pd.DataFrame(trailer_match_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_dORF = MATCH_dORF.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    match_ids = {m.b_info for m in trailer_matches}
    m_match = all_orf_exons['id'].isin(match_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # overlap
    predicted_overlap = predicted_orfs[predicted_orfs['orf_type']=='overlap'][bed_utils.bed12_field_names]
    msg = "Found {} predicted overlap".format(len(predicted_overlap))
    logger.info(msg)

    msg = "Splitting these into exons"
    logger.info(msg)

    predicted_overlap_exons = bed_utils.split_bed12(predicted_overlap, 
                                                    num_cpus=args.num_cpus, 
                                                    progress_bar=False)

    msg = "Finding matches between predicted overlap and transcripts"
    logger.info(msg)

    # anything that overlaps a canonical and that has no specific assigned label...
    overlap_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                 predicted_overlap_exons)

    overlap_ids = {m.b_info for m in overlap_matches}
    m_overlap_matches = predicted_overlap_exons['id'].isin(overlap_ids)
    predicted_overlap_exons = predicted_overlap_exons[~m_overlap_matches]
    msg = "Remaining {} overlap exons".format(len(predicted_overlap_exons))
    logger.info(msg)  

    # Pull that into a dataframe and add APPRIS scores
    overlap_pairs = {(m.b_info, m.a_info) for m in overlap_matches}
    MATCH_OVER = pd.DataFrame(overlap_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_OVER = MATCH_OVER.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    m_match = all_orf_exons['id'].isin(overlap_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # noncoding
    predicted_nc = predicted_orfs[predicted_orfs['orf_type']=='noncoding'][bed_utils.bed12_field_names]
    msg = "Found {} predicted noncoding".format(len(predicted_nc))
    logger.info(msg)

    msg = "Splitting these into exons"
    logger.info(msg)

    predicted_nc_exons = bed_utils.split_bed12(predicted_nc, 
                                               num_cpus=args.num_cpus, 
                                               progress_bar=False)

    msg = "Finding matches between predicted noncoding and transcripts"
    logger.info(msg)

    nc_matches = bed_utils.get_bed_overlaps(noncoding_exons,
                                            predicted_nc_exons,
                                            min_b_overlap=1)

    nc_ids = {m.b_info for m in nc_matches}
    m_nc_matches = predicted_nc_exons['id'].isin(nc_ids)
    predicted_nc_exons = predicted_nc_exons[~m_nc_matches]
    msg = "Remaining {} noncoding exons".format(len(predicted_nc_exons))
    logger.info(msg)  


    # Pull that into a dataframe and add APPRIS scores
    nc_pairs = {(m.b_info, m.a_info) for m in nc_matches}
    MATCH_NC = pd.DataFrame(nc_pairs, columns=['orf_id', 'matching_transcript_ids'])
    MATCH_NC = MATCH_NC.groupby('orf_id', as_index=False).agg({'matching_transcript_ids': lambda x: '|'.join(x)})

    # check
    m_match = all_orf_exons['id'].isin(nc_ids)
    all_orf_exons = all_orf_exons[~m_match]

    # we should have no exons left...
    msg = "Remaining {} ORF exons".format(len(all_orf_exons))
    logger.info(msg)  
    
    # Now we sort everything...

    ALL_MATCH = [MATCH_CDS, MATCH_VAR, MATCH_WITHIN, MATCH_uORF, MATCH_dORF, MATCH_OVER, MATCH_NC]
    ALL_DONE = []

    for df in ALL_MATCH:
        DF = df.copy()
        
        DUPLICATES = DF[DF['matching_transcript_ids'].str.contains('\|')]
        DUPLICATES['transcript_id'] = DUPLICATES.apply(get_transcript_id, axis=1)
        
        DF = DF.merge(DUPLICATES[['orf_id', 'transcript_id']], on='orf_id', how='left')
        DF['transcript_id'] = DF.transcript_id.combine_first(DF.matching_transcript_ids)
        DF = DF.merge(appris, on='transcript_id', how='left')
        ALL_DONE.append(DF)
        
    # If there are still unassigned gene_id, we do NOT remove these from the predictions, as long
    # as we have some transcript id in APPRIS. We will filter later on.

    # last check...
    msg = "Started with {} (unique) predicted ORFs...".format(len(predicted_orfs))
    logger.info(msg)
    msg = "... and matched {} ORFs.".format(sum([len(df) for df in ALL_DONE]))
    logger.info(msg)


    # Write to disk in BED12+ format.
    predicted_orfs.rename(columns={'id': 'orf_id'}, inplace=True)
    all_matches = pd.concat(ALL_DONE)
    predicted_orfs = predicted_orfs.merge(all_matches, on='orf_id', how='left')
    predicted_orfs.to_csv(unique_orfs,
                          compression='gzip', 
                          sep='\t', 
                          index=False)
    

if __name__ == '__main__':
    main()
 
