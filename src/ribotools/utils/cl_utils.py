#! /usr/bin/env python3

"""Provide functionalities for the command line parser.
"""


def add_file_options(parser):

    file_options = parser.add_argument_group("file options")

    file_options.add_argument(
        "-t",
        "--tmp",
        help="""Optional argument: where to write
        temporary files. If not specified, programs-specific tmp will be used.""",
        default=None,
    )

    file_options.add_argument(
        "--overwrite",
        help="""Flag: overwrite existing files.
        Note that certain files may be re-created regardless of the presence or absence
        of this flag. """,
        action="store_true",
    )

    file_options.add_argument(
        "-k",
        "--keep-intermediate-files",
        help="""Flag: unless this flag is
        given, all intermediate files (such as discarded reads) will be deleted, unless the
        [--do-not-call] option is also given.""",
        action="store_true",
    )


def get_file_options_string(args):

    """Create a string containing the relevant flags and options based on the
    file options from the arguments.

    Parameters
    ----------
    args: argparse.Namespace
        The parsed arguments
    """

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    keep_intermediate_files_str = ""
    if args.keep_intermediate_files:
        keep_intermediate_files_str = "--keep-intermediate-files"

    options_str = "{} {}".format(overwrite_str, keep_intermediate_files_str)

    return options_str


# Add functions to interact with and handle options
# for HTSeq (htseq-count).


def add_htseq_options(parser):
    """Add options to a cmd parser to call htseq-count.

    N.B. This is primarily intended for use with the htseq workflow.

    Parameters
    ----------
    parser: argparse.ArgumentParser
        The parser to which the options will be added
    """

    htseq_options = parser.add_argument_group("HTSeq options")

    htseq_options.add_argument(
        "--htseq-options",
        help="""Optional argument: a space-delimited
        list of options to pass to htseq-count. Each option must be quoted separately as in
        "--htseqOption value", using soft quotes, where '--htseqOption'
        is the long parameter name from htseq-count and 'value' is the value given to this parameter.
        If specified, htseq-count options will override default settings.""",
        nargs="*",
        type=str,
    )


def get_htseq_options_string(args):
    """Extract the flags and options specified for htseq-count added with
    add_htseq_options.

    Parameters
    ---------
    args: argparse.Namespace
        The parsed arguments

    Returns
    -------
    htseq_options: string
        a string containing htseq options suitable to pass to another command
    """
    import shlex

    args_dict = vars(args)

    s = ""
    if args_dict["htseq_options"]:
        htseq_option_str = "--htseq-options {}".format(
            " ".join(shlex.quote(op) for op in args_dict["htseq_options"])
        )
        s = "{}".format(" ".join([s, htseq_option_str]))

    return s
