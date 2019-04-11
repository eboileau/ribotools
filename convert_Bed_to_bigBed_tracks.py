#! /usr/bin/env python3

""" Convert Bed files from the Rp-Bp ORF predictions to
bigBed files to be used for visualisation using custom fields.
** single use to rename the predictions using updated labels
and split into used and novel_suspect for manual curation
"""

import sys
import os

import argparse
import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

import yaml
import json
import csv

import bio_utils.bed_utils as bed_utils
import misc.utils as utils
import misc.shell_utils as shell_utils
import misc.pandas_utils as pandas_utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as filenames

# adjust group names otherwise this won't work with 'bedToBigBed'
orf_type_labels_mapping = {
    'Canonical': ['canonical'],
    'Can.(variant)': ['canonical_extended',
                      'novel_canonical_extended',
                      'canonical_truncated'],
    'Can.(oof)': ['within'],
    'lncORF': ['noncoding', 'novel_noncoding'],
    'uORF': ['five_prime', 'five_prime_overlap'],
    'dORF': ['three_prime', 'three_prime_overlap'],
    'Novel': ['novel'],
    'Other': ['novel_within',
              'suspect',
              'novel_suspect',
              'novel_canonical_truncated',
              'novel_three_prime',
              'novel_three_prime_overlap',  # these overlaps are relegated using the relabeled data to suspect
              'novel_five_prime',
              'novel_five_prime_overlap']  # these overlaps are relegated using the relabeled data to suspect
}


color_mapping = {
    'Canonical': '31,119,181',  # #1f77b4
    'Can. (variant)': '174,199,232',  # #aec7e8
    'Can. (oof)': '158,218,229',  # #9edae5
    'lncORF': '197,176,213',  # #c5b0d5
    'uORF': '152,223,138',  # #98df8a
    'dORF': '219,219,141',  # #dbdb8d
    'Novel': '255,187,120'}  # #ffbb78


def _get_bed(input_filename, fields_to_keep, filtered_ids, args):
    """Get Bed12 file and adjust features. The fields
    should correspond to those defined in bed_utils.bed12_field_names
    plus those added by Rp-Bp, and correspond to those passed in via
    [--configure-fields], however there is no check. Some of these fields
    are hard coded below."""

    bed_df = bed_utils.get_bed_df(input_filename)

    # Adjust chrom field.
    if args.add_chr:
        bed_df['seqname'] = 'chr' + bed_df['seqname'].astype(str)
    if args.chr_dict:
        for chrom_old, chrom_new in args.chr_dict.items():
            seqname_m = bed_df['seqname'] == str(chrom_old)
            bed_df.loc[seqname_m, 'seqname'] = str(chrom_new)

    msg = "No. of unique ORFs: {}".format(len(bed_df['id'].unique()))
    logger.info(msg)

    # filter ids
    m_to_discard = bed_df['id'].isin(filtered_ids)
    bed_df = bed_df[~m_to_discard]
    # add group label
    for orf_type in orf_type_labels_mapping.keys():
        bed_df.loc[bed_df['orf_type'].isin(orf_type_labels_mapping[orf_type]), 'orf_group'] = orf_type
    # remove Other
    if not args.keep_other:
        other_m = bed_df['orf_group'] == 'Other'
        bed_df = bed_df[~other_m]

    msg = "No. of unique ORFs (after filtering): {}".format(len(bed_df['id'].unique()))
    logger.info(msg)

    # Sort on the chrom field, and then on the chromStart field.
    bed_df.sort_values(['seqname', 'start'], ascending=[True, True], inplace=True)

    # convert counts to int
    bed_df = bed_df.astype({"x_1_sum": int, "x_2_sum": int, "x_3_sum": int})

    if args.use_color:
        for label, color in color_mapping.items():
            label_m = bed_df['orf_group'] == label
            bed_df.loc[label_m, 'color'] = color

    # remove unused fields, and get order
    bed_df = bed_df[fields_to_keep]

    # Writes bed file to output directory.
    _, name = os.path.split(input_filename)
    filen, ext = os.path.splitext(name)
    if ext == '.gz':
        filen, ext = os.path.splitext(filen)
    output_filename = args.out + filen
    pandas_utils.write_df(bed_df, str(output_filename + '.tmp.bed'), index=False, sep='\t',
                          header=False, do_not_compress=True, quoting=csv.QUOTE_NONE)

    return output_filename


def _convert(bed, bb, use_config_fields, args):
    in_files = [bed, args.chrSizes]
    out_files = [bb]
    if use_config_fields:
        cmd = "bedToBigBed -as={} -type={} -extraIndex={} {} {} {}".format(use_config_fields['as_file'],
                                                                           use_config_fields['bed_type'],
                                                                           "name", bed, args.chrSizes, bb)
        in_files.append(use_config_fields['as_file'])
    else:
        cmd = "bedToBigBed {} {} {}".format(bed, args.chrSizes, bb)
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=True)
    try:
        os.remove(bed)
        msg = "Removing: {}".format(bed)
        logger.info(msg)
    except OSError:
        msg = "Could not remove: {}".format(bed)
        logger.info(msg)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Convert one or more Bed files to bigBed files by calling the executable
        program 'bedToBigBed'. The executable must be available on the user's path, and can be
        downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/'. Unless specified, all 
        samples and conditions in the configuration file are used. ORF types are grouped and
        ORF labeled as "Other" are discarded, unless specified.''')

    parser.add_argument('config', help="The (yaml) config file.")
    parser.add_argument('out', help='''The output directory. All bed files are temporarily re-written
        to this location.''')
    parser.add_argument('chrSizes', help="The 'chrom.sizes' file for the UCSC database.")

    parser.add_argument('--add-chr', help='''If this flag is present then 'chr' will be pre-pended
        to sequence names. This is done before any other changes to sequence names, so this
        must be taken into account if giving a dictionary mapping''', action='store_true')

    parser.add_argument('-d', '--chr-dict', help='''A dictionary mapping of sequence names found
        in the data to the sequence names, as in "chrom.sizes". The format is as follows:
        '{"key":"value"}' ''', type=json.loads)

    parser.add_argument('--configure-fields', help='''A file with comma-separated items (one per line)
        corresponding to fields that will be included in the bigBed file. The field names must
        correspond to the ones used the bed file. Each field name must be followed by
        'type','standard field name', 'description', as needed to generate the AutoSql format (.as)
        file describing these fields. Standard fields must be separated from any extra fields
        by an empty line. See e.g.3 here: https://genome.ucsc.edu/goldenpath/help/bigBed.html.
        One extra index will be created on the name field by default. If multiple bed files are
        passed in, these will be used for all input files.''', required='--use-color' in sys.argv)

    parser.add_argument('-k', '--keep-other', help='''If this flag is present then ORFs labeled
        as "Other" will not be discarded. Note that ORF types/categories are mapped to ORF groups
        according to "local mapping" (to be incorporated in ribo_utils)''', action='store_true')

    parser.add_argument('--use-color', help='''If this flag is present then color (field 9) fields 
        are added. These are currently not configurable, and presumably used only with the ORF 
        predictions from the Rp-Bp pipeline. The [--custom-field] option is required if using 
        colour, even if no extra fields are given.''', action='store_true')

    parser.add_argument('--input-list', help='''A space-delimited list of input files, sample names
        or conditions, each quoted separately. They must either be all files, or sample/condition
        names, and this must be specified with the [--input-type]. Only these will be converted.''',
        nargs='*', required='--input-type' in sys.argv, type=str)

    parser.add_argument('--input-type', help='''The "type" of [--input-list], either f (files),
        s (samples) or c (conditions).''', required='--input-list' in sys.argv, type=str,
        choices=['f', 's', 'c'])

    parser.add_argument('--overwrite', help='''If this flag is present, then existing files
        will be overwritten.''', action='store_true')

    parser.add_argument('--filtered-orfs', help='''A list of ORF ids (one per line) to remove
        from the final prediction set.''', type=str)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    required_keys = [
        'riboseq_data',
        'riboseq_samples',
        'riboseq_biological_replicates'
    ]
    utils.check_keys_exist(config, required_keys)

    filtered_ids = []
    if args.filtered_orfs:
        with open(args.filtered_orfs) as f:
            filtered_ids = f.read().splitlines()

    if args.input_list:
        if args.input_type == 's':
            sample_names = {name: [name] for name in args.input_list}
            condition_names = {}
        elif args.input_type == 'c':
            condition_names = {name: [name] for name in args.input_list}
            sample_names = {}
        else:
            files_only = True
    else:
        files_only = False
        sample_names = config['riboseq_samples']
        condition_names = ribo_utils.get_riboseq_replicates(config)

    # check output path
    if os.path.exists(args.out):
        args.out = os.path.join(args.out, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Quitting.".format(args.out)
        raise OSError(msg)

    use_config_fields = {}
    # generate an AutoSql format (.as) file describing the fields
    if args.configure_fields:
        fields_to_keep = []
        extra_fields = []
        f = open(args.configure_fields, 'r')
        lines = f.readlines()
        f.close()
        as_file = args.out + 'SelectedFields.as'
        f = open(str(as_file), 'w')
        f.write("{} {}\n".format("table", "bedSourceSelectedFields"))
        f.write('''"{}"\n'''.format("Browser extensible data selected fields."))
        f.write("{}\n".format("("))
        n_fields = 0
        for line_no, line in enumerate(lines):
            l = line.strip()
            if not l:
                n_fields = line_no
                break
            fields = l.split(',')
            fields_to_keep.append(fields[0])
            f.write("{}\t{};\t{}\n".format(fields[1], fields[2], fields[3]))
        bed_type = "bed" + str(len(fields_to_keep))
        if n_fields:
            for line in lines[n_fields+1:]:
                l = line.strip()
                fields = l.split(',')
                extra_fields.append(fields[0])
                f.write("{}\t{};\t{}\n".format(fields[1], fields[2], fields[3]))
            bed_type += "+" + str(len(extra_fields))
            fields_to_keep += extra_fields
        f.write("{}\n".format(")"))
        f.close()
        use_config_fields['as_file'] = as_file
        use_config_fields['bed_type'] = bed_type
    else:
        # not sure what happens here...
        fields_to_keep = bed_utils.bed12_field_names

    if files_only:
        for bed_file in args.input_list:
            if not os.path.exists(bed_file):
                msg = "Could not find the bed file: {}. Quitting.".format(bed_file)
                raise FileNotFoundError(msg)
            bed_file_name = _get_bed(bed_file, fields_to_keep, filtered_ids, args)
            bed = bed_file_name + '.tmp.bed'
            bb = bed_file_name + '.bb'
            _convert(bed, bb, use_config_fields, args)

        return

    note_str = config.get('note', None)
    is_unique = not ('keep_riboseq_multimappers' in config)
    # and the smoothing parameters if present
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    # first the samples
    for name in sorted(sample_names.keys()):

        msg = "Processing sample: {}".format(name)
        logger.info(msg)

        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, name, is_unique=is_unique)

        predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], name,
                                                              length=lengths, offset=offsets,
                                                              is_unique=is_unique, note=note_str,
                                                              fraction=fraction,
                                                              reweighting_iterations=reweighting_iterations,
                                                              is_filtered=True)

        if not os.path.exists(predicted_orfs):
            msg = ("Could not find the predictions bed file for {}. ({}). Quitting.".
                   format(name, predicted_orfs))
            raise FileNotFoundError(msg)

        bed_file_name = _get_bed(predicted_orfs, fields_to_keep, filtered_ids, args)
        bed = bed_file_name + '.tmp.bed'
        bb = bed_file_name + '.bb'
        _convert(bed, bb, use_config_fields, args)

    # then conditions
    lengths = None
    offsets = None
    for name in sorted(condition_names.keys()):

        msg = "Processing condition: {}".format(name)
        logger.info(msg)

        predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], name,
                                                              length=lengths, offset=offsets,
                                                              is_unique=is_unique, note=note_str,
                                                              fraction=fraction,
                                                              reweighting_iterations=reweighting_iterations,
                                                              is_filtered=True)

        if not os.path.exists(predicted_orfs):
            msg = ("Could not find the predictions bed file for {}. ({}). Quitting.".
                   format(name, predicted_orfs))
            raise FileNotFoundError(msg)

        bed_file_name = _get_bed(predicted_orfs, fields_to_keep, filtered_ids, args)
        bed = bed_file_name + '.tmp.bed'
        bb = bed_file_name + '.bb'
        _convert(bed, bb, use_config_fields, args)


if __name__ == '__main__':
    main()
