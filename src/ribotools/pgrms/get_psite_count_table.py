#! /usr/bin/env python3

"""Get p-site count table from Rp-Bp predictions.
The script is not very efficient/clever, we initially
created different tables to separate CDS and matched uORFs, but
this is not present anymore. Currently, all counts are aggregated
at gene level (TODO: transcript level), and we either use CDS
or selected ORF types from the Rp-Bp labels.

Required:
    long table ORF predictions
    unique table ORF+gene ids (APPRIS)
"""


import sys
import os
import argparse
import logging
import pandas as pd
import numpy as np
import yaml


import pbio.utils.bed_utils as bed_utils

import pbio.ribo.ribo_utils as ribo_utils
import pbio.ribo.ribo_filenames as filenames

import pbio.misc.shell_utils as shell_utils
import pbio.misc.logging_utils as logging_utils

from rpbp.defaults import default_num_cpus

logger = logging.getLogger(__name__)


# default filenames
all_orfs_file = 'filtered.predicted-orfs.bed.gz'
unique_orfs_file = 'unique-orfs.appris.bed.gz'
final_orfs_file = 'orfs.bed.gz'
   
   
def get_counts_cds(groupby_condition, args):
    
    # groupby_condition after groupby_gene (ie. replicate)
    gene_id = groupby_condition['gene_id'].unique()[0]
    
    # hard coded here: CDS only, one CDS per gene for any given condition
    #if len(groupby_condition) > 1:
    try:
        principal = groupby_condition[groupby_condition.label.str.startswith('PRINCIPAL')]
        principal = principal.sort_values(by=['appris_score', 'frame_count'], ascending=[False, False])
    except:
        pass
    finally:
        # either no principal, or there were some missing values...
        # in any case, we still are "allowed" to sort the scores
        if principal.empty: 
            principal = groupby_condition.sort_values(by=['appris_score', 'frame_count'], ascending=[False, False])

    principal = principal.iloc[0]
    counts = principal.frame_count
    if args.normalise:
        counts /= principal.orf_len

    ret = pd.Series([gene_id, counts])
    
    return ret


def get_tables_cds(groupby_gene, args):
    
    ret = groupby_gene.groupby('condition').apply(get_counts_cds, args)
    pivoted = ret.pivot_table(ret, columns='condition', index=0)

    return pivoted 
   
   
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''Create in-frame p-site count table
                                     using the counts from the rpbp predictions.''')

    parser.add_argument('dirloc', help="Path to the working directory for I/O files.")
    parser.add_argument('output', help="Output name (csv.gz)")
    parser.add_argument('config', help="The (yaml) configuration file for filtering merged replicates.")

    parser.add_argument('--orfs-long', help="""If given, will use this file as input (all ORFs).""", type=str)
    parser.add_argument('--orfs-unique', help="""If given, will use this file as input (ORFs+APPRIS).""", type=str)
    
    parser.add_argument('-minp', '--min_psite', help="""Min. number of in-frame p-sites""", type=int, default=10)
    parser.add_argument('-minrep', '--min-rep', help="""Min. number of replicates for evidence""", type=int, default=1)
    
    parser.add_argument('-n', '--normalise', help="""If this flag is present then number of
                        p-sites is normalised by ORF length, before aggregation at gene level.""", 
                        action='store_true')
    
    parser.add_argument('--types', help="""A space-delimited list of ORF types to consider. If 'canonical' 
                        only, we use the principal transcript, and the ORF with largest number of p-sites
                        before normalisation, otherwise we use everything (TODO:deal with ORF types).""",
                        nargs='*', type=str, default=['canonical']) # ONLY CANONICAL NOW
    
    parser.add_argument('-p', '--num-cpus', help='''The number of CPUs to use.''', 
                        type=int, default=default_num_cpus)
    

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[get-psite-table]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    basename = config.get("project_name", args.name)
    
    msg = 'Removing merged replicates.'
    logger.info(msg)
    condition_name_map = ribo_utils.get_condition_name_map(config)
    conditions = []
    for key in ribo_utils.get_riboseq_replicates(config).keys():
        conditions.append(condition_name_map[key])
    
    # output
    output_file = os.path.join(args.dirloc, args.output)
        
    # find all ORF predictions
    if not args.orfs_long:
        orfs_long_file = '{}.{}'.format(basename, all_orfs_file)
        args.orfs_long = os.path.join(args.dirloc, orfs_long_file)
    if not args.orfs_unique:
        orfs_unique_file = '{}.{}'.format(basename, unique_orfs_file)
        args.orfs_unique = os.path.join(args.dirloc, orfs_unique_file)
    
    appris = pd.read_csv(args.orfs_unique, sep='\t')
    predicted_orfs = bed_utils.read_bed(args.orfs_long)
    
    # Remove "merged replicates"...
    predicted_orfs = predicted_orfs[~predicted_orfs['condition'].isin(conditions)]
    
    # This will make things easier below...
    predicted_orfs.rename(columns={'id': 'orf_id'}, inplace=True)
    
    # filter based on types
    m_types = predicted_orfs['orf_type'].isin(args.types)
    predicted_orfs = predicted_orfs[m_types]
    
    msg = 'Replicates: using {}'.format(', '.join(predicted_orfs.condition.unique()))
    logger.info(msg)
    msg = 'ORF types: using {}'.format(', '.join(predicted_orfs.orf_type.unique()))
    logger.info(msg)
    
    
    # We add APPRIS and egne info to the replicate predictions (P-sites per condition).
    len_before = len(predicted_orfs)
    appris_to_keep = ['orf_id', 'matching_transcript_ids', 'transcript_id', 'gene_id', 'appris_score', 'label']
    predicted_orfs = predicted_orfs.merge(appris[appris_to_keep], on='orf_id', how='left', validate='m:1')

    len_after = len(predicted_orfs)
    if len_before != len_after:
        msg = '{} entries before, and {} after merge...'.format(len_before, len_after)
        logger.critical(msg)

    # We remove unassigned gene ids.
    predicted_orfs = predicted_orfs[~predicted_orfs['gene_id'].isna()]
    msg = "Found {} unique ORFs before filtering (with valid gene_id)...".format(len(predicted_orfs['orf_id'].unique()))
    logger.info(msg)

    # Filter using criteria defined above based on minimum count and minimum number of replicates.
    predicted_orfs = predicted_orfs[predicted_orfs['frame_count'] >= args.min_psite]
    predicted_orfs = predicted_orfs.groupby('orf_id').filter(lambda x: len(x) >= args.min_rep)

    msg = "Found {} unique ORFs after filtering (with valid gene_id)...".format(len(predicted_orfs['orf_id'].unique()))
    logger.info(msg)

    msg = "Found {} unique genes.".format(len(predicted_orfs['gene_id'].unique()))
    logger.info(msg)
    
    # write the final set of ORFs
    to_subset = bed_utils.bed12_field_names + ['orf_len', 'orf_type']
    to_keep = to_subset + ['matching_transcript_ids', 'transcript_id', 'gene_id', 'appris_score', 'label']
    final_list = predicted_orfs.copy()
    final_list.rename(columns={'orf_id': 'id'}, inplace=True)
    final_list = final_list[to_keep]
    # keep the unique ORFs only based on BED12+orf_len+orf_type
    final_list.drop_duplicates(subset=to_subset, inplace=True)
    
    final_orfs_file = '{}-{}'.format(basename, final_orfs_file)
    final_list.to_csv(os.path.join(args.dirloc, final_orfs_file),
                      compression='gzip', 
                      sep='\t', 
                      index=False)
    
    # get the counts (CDS only for now)
    one_table = predicted_orfs.groupby('gene_id').apply(get_tables_cds, args)
    
    # We now split the one table

    table = one_table.xs(1, axis=1).fillna(0)
    table = table.droplevel(level=0)
    table.index.name = None
    table.columns.name = None

    table.to_csv(os.path.join(output_file, 
                              sep=',', 
                              index=True, 
                              compression='gzip')
    

if __name__ == '__main__':
    main()
 
    
