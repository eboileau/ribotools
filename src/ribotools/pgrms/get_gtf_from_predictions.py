#! /usr/bin/env python3

"""Convert unique Rp-Bp ORF predictions to GTF"""

import csv
import argparse
import logging

import pandas as pd

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.utils.gtf_utils as gtf_utils
import pbiotools.utils.bed_utils as bed_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.pandas_utils as pandas_utils

logger = logging.getLogger(__name__)

ATTRS = ["gene_id", "gene_name", "gene_biotype", "transcripts", "biotype", "orf_type"]
FIELDS = bed_utils.bed12_field_names + ATTRS


def _clean_attrs(attrs):
    attr_list = attrs.split(";")
    return ";".join([a for a in attr_list if a and "nan" not in a])


def _get_gtf_entries(bed_entry, source: str, id_attribute: str = "transcript_id"):
    """Same as gtf_utils,, but pass the id_attribute!"""

    gtf_exons = gtf_utils._get_gtf_entries(bed_entry, "exon", source, id_attribute)

    if bed_entry["thick_start"] > -1:
        bed_entry_cds = bed_utils.retain_thick_only(bed_entry)
        gtf_cds = gtf_utils._get_gtf_entries(bed_entry_cds, "CDS", source, id_attribute)

        gtf_exons = pd.concat([gtf_exons, gtf_cds])

    gtf_exons = gtf_exons.sort_values("start")
    gtf_exons = gtf_exons.reset_index(drop=True)
    return gtf_exons[gtf_utils.gtf_field_names]


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Convert Rp-Bp ORF predictions (BED12+) into GTF. The
        "thick_start" and "thick_end" fields are used to determine the CDS
        GTF entries. Only creates "exon" and "CDS" entries. Additional
        columns are included as attributes.""",
    )

    parser.add_argument(
        "bed",
        help="The bed12 file. It must conform to the "
        "style expected by utils.bed_utils.",
    )
    parser.add_argument(
        "out",
        help="The (output) gtf file.",
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="The number of CPUs to use",
        type=int,
        default=1,
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading bed file"
    logger.info(msg)
    bed = bed_utils.read_bed(args.bed)

    # rename/wrangle attributes
    try:
        bed = bed[FIELDS]
        bed.rename(
            columns={"transcripts": "transcript_id", "biotype": "transcript_biotype"},
            inplace=True,
        )
        bed["orf_id"] = bed["id"]
        id_attribute = None
        source = "Rp-Bp"
    except:
        msg = (
            f"Input BED file {args.bed} has unexpected or missing fields! "
            f"Additional required fields are {ATTRS}. Truncating to BED12..."
        )
        logger.warning(msg)
        bed = bed.iloc[:, :12]
        id_attribute = "id"
        source = None

    try:
        bed.drop_duplicates(subset=bed_utils.bed12_field_names, inplace=True)
    except KeyError:
        msg = f"Input BED file {args.bed} cannot be processed, missing fields: {bed_utils.bed12_field_names}"
        logger.error(msg)
        raise KeyError(msg)

    msg = "Expanding BED entries to GTF entries"
    logger.info(msg)

    gtf_entries = parallel.apply_parallel(
        bed,
        args.num_cpus,
        _get_gtf_entries,
        source,
        id_attribute,
        progress_bar=True,
    )

    msg = "Joining GTF entries into large data frame"
    logger.info(msg)

    gtf_entries = pd.concat(gtf_entries)
    gtf_entries["attributes"] = gtf_entries["attributes"].apply(_clean_attrs, 1)

    msg = "Sorting GTF entries"
    logger.info(msg)

    gtf_entries = gtf_entries.sort_values(["seqname", "start", "end"])
    gtf_entries = gtf_entries.reset_index(drop=True)

    msg = "Writing GTF to disk"
    logger.info(msg)

    start_field = gtf_entries.columns[3]
    end_field = gtf_entries.columns[4]

    gtf_entries[start_field] = gtf_entries[start_field].astype(int)
    gtf_entries[end_field] = gtf_entries[end_field].astype(int)

    pandas_utils.write_df(
        gtf_entries,
        args.out,
        index=False,
        sep="\t",
        header=None,
        do_not_compress=True,
        quoting=csv.QUOTE_NONE,
    )


if __name__ == "__main__":
    main()
