"""Define all default parameters and options.

** NOTE: Modifying the default parameter values set in this file
will have NO effect. To override a given parameter value, you must
either provide it via command argument when calling the Rp-Bp
pipeline, or else set the desired value in the configuration file.
Please refer to the documentation: https://rp-bp.readthedocs.io/en/latest/usage-instructions.html#


Default options for external programs (Flexbar, STAR) are
overridden via command line. Currently, call to Bowtie2 is not customisable.

Default parameters for Rp-Bp (default options for Rp-Bp scripts, including
shared MCMC options, but not processing options) are overridden by
providing the option key-value pair in the configuration file.

Processing options (parallel processing, logging options) are given via
command line.
"""

# htseq-count options must be passed via the [--htseq-options]

# In addition, when passing [--run-all], use [--stranded] to override
# default htseq_options

# stranded=yes: read has to be mapped to the same strand as the feature (single-end) or fr or second strand
# stranded=reverse: read has to be mapped to the opposite strand as the feature (single-end) or rf or first strand

# Because the [--additional-attr] must be repeated for each attribute, we cannot currently use a list
# of default attributes. Use e.g. --htseq-options "--additional-attr gene_name" "--additional-attr transcript_id"

htseq_executable = "htseq-count"
htseq_options = {
    "format": "bam",
    "stranded": "no",
    "type": "CDS",  # feature must be 3rd field of GTF file
    "idattr": "gene_id",
    "additional-attr": "gene_name",
    "mode": "intersection-nonempty",
    "secondary-alignments": "ignore",
    "supplementary-alignments": "ignore",
}
