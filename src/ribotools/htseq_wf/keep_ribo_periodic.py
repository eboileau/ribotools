#! /usr/bin/env python3

"""Filter read lengths out of a BAM file

Functions:
    filter_non_periodic_reads
"""

import argparse
import logging
import sys

import tqdm

import pysam

import pbiotools.utils.bam_utils as bam_utils
import pbiotools.misc.logging_utils as logging_utils

import ribotools.utils.cl_utils as clu

logger = logging.getLogger(__name__)


def filter_non_periodic_reads(alignments, lengths):

    for a in alignments:
        if a.qlen in lengths:
            yield a
        else:
            yield None


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Filter non-periodic read lengths out of
        a BAM file, but does not shift the start positions.""",
    )

    parser.add_argument("bam", help="The (BAM) file to filter.")

    parser.add_argument("out", help="The filtered (BAM) file.")

    parser.add_argument(
        "-l",
        "--lengths",
        help="If any values are given, then only reads "
        "which have those lengths will be included in the final BAM file.",
        type=int,
        default=[],
        nargs="+",
    )

    parser.add_argument(
        "--do-not-call",
        help="Do not execute the program (dry run).",
        action="store_true",
    )

    clu.add_file_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[keep-ribo-periodic]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    msg = "Reading the alignments"
    logger.info(msg)

    bam = pysam.AlignmentFile(args.bam)
    alignments = bam.fetch()
    num_alignments = bam.count()

    # create the output file
    out_bam = pysam.AlignmentFile(args.out, "wb", template=bam)

    msg = "Filtering the alignments"
    logger.info(msg)

    for a in tqdm.tqdm(
        filter_non_periodic_reads(alignments, args.lengths),
        leave=True,
        file=sys.stdout,
        total=num_alignments,
    ):
        if a is not None:
            out_bam.write(a)

    out_bam.close()

    # create the bamtools index if it does not already exists
    bam_utils.index_bam_file(args.out, args)


if __name__ == "__main__":
    main()
