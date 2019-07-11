#! /usr/bin/env python3

"""Convert BED files from the Rp-Bp ORF predictions to bigBed
tracks for visualisation.

Note: Requires the executable "bedToBigBed" and the "chrom.sizes"
      file which can be obtained with the "fetchChromSizes" script.
      The chrom names from the predictions must match those
      of chrom.sizes (UCSC), otherwise they must be re-written by
      passing the right options.

Functions:
    _get_bed
    _convert
"""

import sys
import os
import argparse
import logging
import yaml
import json
import csv

import pbio.misc.utils as utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.pandas_utils as pandas_utils
import pbio.misc.logging_utils as logging_utils

import pbio.utils.bed_utils as bed_utils

import pbio.ribo.ribo_utils as ribo_utils
import pbio.ribo.ribo_filenames as filenames

from rpbp.defaults import metagene_options

logger = logging.getLogger(__name__)


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


def _get_bed(input_filename, pretty_name, fields_to_keep, args):

    """Get BED12 file and adjust features. The fields must match those
    defined in bed_utils.bed12_field_names plus fields passed via
    [--configure-fields], however there is currently no check.
    """

    bed_df = bed_utils.get_bed_df(input_filename)

    # Adjust chrom field
    if args.add_chr:
        bed_df['seqname'] = 'chr' + bed_df['seqname'].astype(str)
    if args.chr_dict:
        for chrom_old, chrom_new in args.chr_dict.items():
            seqname_m = bed_df['seqname'] == str(chrom_old)
            bed_df.loc[seqname_m, 'seqname'] = str(chrom_new)

    msg = "No. of unique ORFs: {}".format(len(bed_df['id'].unique()))
    logger.info(msg)

    # add display label
    for orf_type, labels in ribo_utils.orf_type_labels_mapping.items():
        bed_df.loc[bed_df['orf_type'].isin(labels), 'orf_category'] = display_name_map[orf_type]

    # just to make sure...
    remove_m = bed_df['orf_type'].isna()
    if not args.keep_other:
        other_m = bed_df['orf_category'] == 'Other'
        remove_m = remove_m | other_m
    bed_df = bed_df[~remove_m]

    msg = "No. of unique ORFs (after filtering): {}".format(len(bed_df['id'].unique()))
    logger.info(msg)

    # Sort on the chrom field, and then on the chromStart field.
    bed_df.sort_values(['seqname', 'start'], ascending=[True, True], inplace=True)
    # convert counts to int
    bed_df = bed_df.astype({"x_1_sum": int, "x_2_sum": int, "x_3_sum": int})

    if args.use_color:
        for label, color in color_mapping.items():
            label_m = bed_df['orf_category'] == label
            bed_df.loc[label_m, 'color'] = color

    # remove unused fields, and get order
    bed_df = bed_df[fields_to_keep]

    # Writes bed file to output directory
    pretty_name = pretty_name + '.orfs'
    output_filename = os.path.join(args.dirloc, pretty_name)

    pandas_utils.write_df(bed_df,
                          str(output_filename + '.tmp.bed'),
                          index=False,
                          sep='\t',
                          header=False,
                          do_not_compress=True,
                          quoting=csv.QUOTE_NONE)

    return output_filename


def _convert(bed, bb, use_config_fields, args):

    in_files = [bed, args.chrSizes]
    out_files = [bb]
    if use_config_fields:
        cmd = "bedToBigBed -as={} -type={} -extraIndex={} {} {} {}".format(use_config_fields['as_file'],
                                                                           use_config_fields['bed_type'],
                                                                           "name",
                                                                           bed,
                                                                           args.chrSizes,
                                                                           bb)
        in_files.append(use_config_fields['as_file'])
    else:
        cmd = "bedToBigBed {} {} {}".format(bed, args.chrSizes, bb)

    shell_utils.call_if_not_exists(cmd,
                                   out_files,
                                   in_files=in_files,
                                   overwrite=args.overwrite,
                                   call=True)
    try:
        os.remove(bed)
        msg = "Removing: {}".format(bed)
        logger.info(msg)
    except OSError:
        msg = "Could not remove: {}".format(bed)
        logger.info(msg)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Convert one or more BED files to bigBed 
        files by calling the executable program 'bedToBigBed'. The executable must be 
        available on the user's path, and can be downloaded from 
        'http://hgdownload.soe.ucsc.edu/admin/exe/'.""")

    parser.add_argument('config', help="The (yaml) config file.")

    parser.add_argument('dirloc', help="""The output directory. All BED files are 
        temporarily re-written to this location.""")

    parser.add_argument('chrSizes', help="The 'chrom.sizes' file for the UCSC database.")

    parser.add_argument('--add-chr', help="""If this flag is present then 'chr' will be pre-pended
        to sequence names. This is done before any other changes to sequence names, so this
        must be taken into account if giving a dictionary mapping""", action='store_true')

    parser.add_argument('-d', '--chr-dict', help="""A dictionary mapping of sequence names found
        in the data to the sequence names, as in "chrom.sizes". The format is as follows:
        '{"key":"value"}'""", type=json.loads)

    parser.add_argument('--configure-fields', help="""A file with comma-separated items (one per line)
        corresponding to fields that will be included in the bigBed file. The field names must
        correspond to the ones used the BED file. Each field name must be followed by
        'type', 'standard field name', 'description', as needed to generate the AutoSql format (.as)
        file describing these fields. Standard fields must be separated from any extra fields
        by an empty line. See e.g.3 here: https://genome.ucsc.edu/goldenpath/help/bigBed.html.
        One extra index will be created on the name field by default. If multiple BED files are
        passed in, these will be used for all input files.""", required='--use-color' in sys.argv)

    parser.add_argument('--use-color', help="""If this flag is present then color (field 9) fields 
        are added. These are currently not configurable, and presumably used only with the ORF 
        predictions from the Rp-Bp pipeline. The [--configure-field] option is required if using 
        colour, even if no extra fields are given.""", action='store_true')

    parser.add_argument('-k', '--keep-other', help="""If this flag is present then ORFs labeled
        as "Other" will be included. They are discarded by default.""", action='store_true')

    parser.add_argument('-a', '--all-replicates', help="""If this flag is present then bigBed files
        are created for all replicates in the config file, in addition to the merged replicates or
        conditions. By default, only the latter are created.""", action='store_true')

    parser.add_argument('--input-list', help="""A space-delimited list of input files, sample names
        or conditions, each quoted separately. They must either be all files, or sample/condition
        names, and this must be specified with the [--input-type]. Only these will be converted.""",
        nargs='*', required='--input-type' in sys.argv, type=str)

    parser.add_argument('--input-type', help="""The 'type' of [--input-list], either f (files),
        s (samples) or c (conditions).""", required='--input-list' in sys.argv, type=str,
        choices=['f', 's', 'c'])

    parser.add_argument('--overwrite', help='''If this flag is present, then existing files
        will be overwritten.''', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[create-bigBed-tracks]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    required_keys = [
        'riboseq_data',
        'riboseq_samples',
        'riboseq_biological_replicates'
    ]
    utils.check_keys_exist(config, required_keys)

    sample_name_map = ribo_utils.get_sample_name_map(config)
    condition_name_map = ribo_utils.get_condition_name_map(config)

    files_only = False
    sample_names = {}
    condition_names = {}
    if args.input_list:
        if args.input_type == 's':
            args.all_replicates = True
            sample_names = {name: [name] for name in args.input_list}
        elif args.input_type == 'c':
            condition_names = {name: [name] for name in args.input_list}
        else:
            files_only = True
    else:
        sample_names = config['riboseq_samples']
        condition_names = ribo_utils.get_riboseq_replicates(config)

    # check output path
    if os.path.exists(args.dirloc):
        args.dirloc = os.path.join(args.dirloc, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.dirloc)
        raise OSError(msg)

    use_config_fields = {}
    # generate an AutoSql format (.as) file describing the fields
    if args.configure_fields:
        fields_to_keep = []
        extra_fields = []
        f = open(args.configure_fields, 'r')
        lines = f.readlines()
        f.close()
        as_file = args.dirloc + 'SelectedFields.as'
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
        fields_to_keep = bed_utils.bed12_field_names # not sure what happens here ...
        msg = """Currently no default fields, the [--configure-fields] option 
              must be used. Terminating."""
        logger.critical(msg)

    if files_only:
        for bed_file in args.input_list:
            if not os.path.exists(bed_file):
                msg = "Could not find the bed file: {}. Terminating.".format(bed_file)
                raise FileNotFoundError(msg)
            bed_file_name = _get_bed(bed_file, fields_to_keep, args)
            bed = bed_file_name + '.tmp.bed'
            bb = bed_file_name + '.bb'
            _convert(bed, bb, use_config_fields, args)

        return

    note_str = config.get('note', None)
    is_unique = not ('keep_riboseq_multimappers' in config)
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    if args.all_replicates:

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
            bed_file_name = _get_bed(predicted_orfs, pretty_name, fields_to_keep, args)
            bed = bed_file_name + '.tmp.bed'
            bb = bed_file_name + '.bb'
            _convert(bed, bb, use_config_fields, args)

    # the merged replicates or conditions are always created
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
        bed_file_name = _get_bed(predicted_orfs, pretty_name, fields_to_keep, args)
        bed = bed_file_name + '.tmp.bed'
        bb = bed_file_name + '.bb'
        _convert(bed, bb, use_config_fields, args)


if __name__ == '__main__':
    main()
