#! /usr/bin/env python3

"""Call htseq-count on a list of samples
"""

import sys
import argparse
import logging
import yaml

import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.utils as utils
import pbio.misc.slurm as slurm

import pbio.ribo.ribo_filenames as filenames
import pbio.ribo.ribo_utils as ribo_utils


logger = logging.getLogger(__name__)

stringtie_executable = 'stringtie'



# now if calling run-all then how we can pass the strandedeness???
# to collect all files, we need to know if unique, then if filtered by lenth
# for both Ribo and RNA independently
# then we need to add all files into the command, and also list the sample
# names (with dict nice) in a separate file

# BETTER, WORK WITH DOWNSREAM DESEQ2, WHICH CAN USE OUTPUT FROM HTSEQ COUNT
# SO IN FACT IN THE PIPLEINE, WE CALL FOR EACH SAMPLE, WHEN READY
# NEED TO PRIVUDE OUTPUT DIR OR ADD TO FILKENAMES
# PASS NAME, CONFIG, HTSEQ ARGS,

# PASS skip_periodicity_estimation OR SIMILAR
# AND trim_rna_to_max_fragment_size

# MAYBE ADD STRANDEDNESS AS OPTION, AND SORT HOW TO IF RUN-ALL
# THEN ADD TO DEFAULTS BEFORE CONSTRUCTING STRING AND PASS,...
# SO HERE JUST HTSEQ-OPTIONS, AND GET FINAL OPTIONS PGRMS


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Call Stringtie and extract read counts [-e], require [-G]. With 
        this option, reads with no reference are skipped. Also use [-C].""")

    parser.add_argument('seq', choices=['rna', 'ribo'])
    parser.add_argument('config', help="The yaml config file.")
    parser.add_argument('name', help="The name of the dataset.")
    parser.add_argument('gtf', help="Stringtie output GTF file")

    parser.add_argument('--strandedness', help="""Library strandedness. If unstranded,
            the "XS" BAM tag must be included to indicate the genomic strand. No check
            will be performed.""")

    parser.add_argument('--not-periodic', help="""Flag: non-periodic read 
            lengths are NOT filtered out; this is to get the right file names 
            (length-). For ribo only.""", action='store_true')

    parser.add_argument('--ribo-config', help="""Optional argument: the Ribo config file
            when seq is rna and rna reads have been trimmed to max ribo fragment lengths.
            If reads are trimmed then this needs to be given, otherwise the program will
            not find the alignment files. In addition, the rna config file must include 
            "matching_samples".""", type=str)

    clu.add_file_options(parser)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    pgrm_utils.add_flexbar_options(parser)
    clu.add_htseq_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)


    args = parser.parse_args()
    logging_utils.update_logging(args)

    # if using slurm, submit the script
    cmd = ' '.join(sys.argv)
    if args.use_slurm:
        slurm.check_sbatch(cmd, args=args)
        return
    logger.info(cmd)

    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    config = yaml.load(open(args.config))
    note = config.get('note', None)
    keep_key = 'keep_' + str(args.seq) + 'seq_multimappers'
    data_key = str(args.seq) + 'seq_data'
    is_unique = not (keep_key in config)

    strandedness = ''
    if args.strandedness != 'un':
        strandedness = '--{}'.format(args.strandedness)

    if args.seq == 'ribo' and not args.not_periodic:
        # get the lengths, we do not use the offsets
        lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            args.name,
            isoform_strategy=args.isoform_strategy,
            is_unique=is_unique
        )

        if len(lengths) == 0:
            msg = "No periodic read lengths and offsets were found!"
            logger.critical(msg)
            return

    elif args.seq == 'rna' and args.ribo_config:
        config_keys = ['matching_samples']
        utils.check_keys_exist(config, config_keys)
        matching_ribo_sample = config['matching_samples'][args.name]

        ribo_config = yaml.load(open(args.ribo_config))
        is_unique_ribo = not ('keep_riboseq_multimappers' in ribo_config)

        # get the lengths, we don't need the offsets
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

    else:
        lengths = None

    bam = filenames.get_seq_bam(
        args.seq,
        config[data_key],
        args.name,
        is_unique=is_unique,
        length=lengths,
        isoform_strategy=args.isoform_strategy,
        stranded=args.strandedness,
        note=note
    )

    # default to stringtie output location of gtf
    cov_ref = args.gtf.split('gtf')[0] + 'cov_refs.gtf'

    cmd = "{} -e {} -p {} -G {} -C {} -o {} {}".format(
        stringtie_executable,
        strandedness,
        args.num_cpus,
        config['gtf'],
        cov_ref,
        args.gtf,
        bam)

    in_files = [config['gtf'], bam]
    out_files = [args.gtf, cov_ref]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call,
                                   keep_delete_files=keep_delete_files)


if __name__ == '__main__':
    main()
