
import sys

"""Define all default parameters and options.

** NOTE: Modifying the default parameter values set in this file
will have NO effect. To override a given parameter value, you must
either provide it via command argument, or else set the desired value 
in the configuration file.

Default options for external programs (Flexbar, STAR, HTSeq) are
overridden via command line. Currently, call to Bowtie2 is not customisable.

Default parameters for periodicity estimation (default for scripts, including
shared MCMC options, but not processing options) are overridden by
providing the option key-value pair in the configuration file. Default values 
for these parameters are the same as those set in the rpbp package.

Processing options (parallel processing, logging options) are given via
command line.
"""


# default: processing
# overridden via command line arguments

default_num_cpus = 1
default_mem = '2G'

default_num_groups = 100  # currently cannot be overridden

# default: periodicity estimation (metagene)
# overridden via command line arguments

default_read_files_command = "zcat"
if sys.platform.startswith("darwin"):
    default_read_files_command = "gzcat"

# NOTES:
#   --limitBAMsortRAM 0 (with --genomeLoad NoSharedMemory)
#   limitBAMsortRAM is set to args.mem or default_mem, if not given, at run time

#   --outTmpDir is set to args.tmp at run time, else default STAR tmp is used (current directory)
#   if given, STAR args.tmp will be created

star_executable = 'STAR'
star_options = {
    'readFilesCommand': default_read_files_command,
    'limitBAMsortRAM': 0,
    'alignIntronMin': 20,
    'alignIntronMax': 100000,
    'outFilterMismatchNmax': 1,
    'outFilterMismatchNoverLmax': 0.04,  # would require a minimum *mapped* length of 25 given max mismatch of 1
    'outFilterType': 'BySJout',
    'outFilterIntronMotifs': 'RemoveNoncanonicalUnannotated',
    'outSAMattributes': ['AS', 'NH', 'HI', 'nM', 'MD'],
    'outSAMtype': 'BAM SortedByCoordinate',
    'sjdbOverhang': 33,  # roughly 90 percentile of a large dataset of varying Ribo-seq fragment lengths
    'seedSearchStartLmaxOverLread': 0.5,  # default seedSearchStartLmax normalised to read length
    'winAnchorMultimapNmax': 100  # increase number of loci anchors are allowed to map to
}

# leave outFilterMultimapNmax to default 20, we filter the multimappers afterwards if desired

flexbar_options = {
    'max-uncalled': 1,
    'pre-trim-left': 0,
    'qtrim-format': 'sanger',
    'qtrim': 'TAIL',
    'qtrim-threshold': 10,
    'zip-output': 'GZ'
}

# htseq-count options must be passed via the [htseq-options]

# [--strandedness] whether the RNA-seq data is from a strand-specific assay
# and [--run-all] is passed, must

# strandedness option must be passed as argument when calling the pipeline,
# otherwise default is unstranded.
# stranded=yes: read has to be mapped to the same strand as the feature (single-end) or fr or second strand
# stranded=reverse: read has to be mapped to the opposite strand as the feature (single-end) or rf or first strand

# all other htseq-count options must be passed via the [htseq-options]

htseq_executable = 'htseq-count'
htseq_options = {
    'format': 'bam',
    'stranded': 'no',
    'type': 'CDS',  # feature must be 3rd field of GTF file
    'idattr': 'gene_id',
    'additional-attr': ['gene_name', 'transcript_id'],
    'mode': 'intersection-nonempty',
    'secondary-alignments': 'ignore',
    'supplementary-alignments': 'ignore'
}

# default: general (for STAR and HTSeq)
# overridden via config file using the key: value pair "gff_spec" : "gff2" or "gff3"
# where gff2 corresponds to the standard gtf format.

# By default, htseq-count expects as gtf file (gff2 specs). If using a gff file,
# then appropriate options must be passed to htseq-count (e.g. modifying type and idattr).
# See also https://htseq.readthedocs.io/en/release_0.11.1/count.html

# The appropriate flags are passed to inform STAR if using a gff file.
# Note that gff files do not generally include the parent gene for each exon,
# but only a parent transcript. This is likely to be required if using --quantMode GeneCounts.
default_gff_spec = 'gff2'

# default: periodicity estimation (metagene)
# overridden via config file

metagene_options = {
    'metagene_start_upstream': 300,
    'metagene_start_downstream': 300,
    'metagene_end_upstream': 300,
    'metagene_end_downstream': 300,
    'periodic_offset_start': -20,
    'periodic_offset_end': 0,
    'metagene_profile_length': 21,
    'seed': 8675309,
    'chains': 2,
    'metagene_iterations': 500,
    'min_metagene_profile_count': 1000,
    'min_metagene_bf_mean': 5,
    'max_metagene_bf_var': None,
    'min_metagene_bf_likelihood': 0.5
}
