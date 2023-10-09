import os
import rpbp.ribo_utils.filenames as filenames


class _return_key_dict(dict):
    def __missing__(self, key):
        return key


def get_stranded_library_string(stranded):
    s = ""
    if (stranded is not None) and (stranded in ["fr", "rf"]):
        s = ".stranded-{}".format(stranded)

    return s


def get_seq_bam(seq, base, name, **kwargs):
    s = ""
    if seq == "rna":
        s = get_rnaseq_bam(base, name, **kwargs)
    if seq == "ribo":
        s = filenames.get_riboseq_bam(base, name, **kwargs)
    return s


def get_rnaseq_bam_base(
    rnaseq_base,
    name,
    length=None,
    is_unique=False,
    stranded=None,
    note=None,
):

    unique = filenames.get_unique_string(is_unique)
    l = filenames.get_length_string(length)
    sl = get_stranded_library_string(stranded)
    n = filenames.get_note_string(note)

    bam_base = "{}{}{}{}{}".format(name, n, unique, sl, l)

    rnaseq_bam_path = get_rnaseq_bam_path(rnaseq_base)
    bam_base = os.path.join(rnaseq_bam_path, bam_base)
    return bam_base


def get_rnaseq_bam(rnaseq_base, name, **kwargs):

    s = get_rnaseq_bam_base(rnaseq_base, name, **kwargs)
    s = s + ".bam"
    return s


def get_rnaseq_bam_path(base_path):
    return os.path.join(
        base_path,
        "without-rrna-mapping",
    )


def get_riboseq_bam_base(riboseq_base, name, **kwargs):
    return filenames.get_riboseq_bam_base(riboseq_base, name, **kwargs)


def get_without_adapters_base(base_path, name, note=None):
    return filenames.get_without_adapters_base(base_path, name, note=note)


def get_without_adapters_fastq(base_path, name, note=None):
    return filenames.get_without_adapters_fastq(base_path, name, note=note)


def get_with_rrna_fastq(base_path, name, note=None):
    return filenames.get_with_rrna_fastq(base_path, name, note=note)


def get_without_rrna_fastq(base_path, name, note=None):
    return filenames.get_without_rrna_fastq(base_path, name, note=note)


def get_gtf(config):
    return filenames.get_gtf(config)


def get_count_table(seq_base, name, length=None, is_unique=False, note=None):

    unique_str = filenames.get_unique_string(is_unique)
    length_str = filenames.get_length_string(length)
    note_str = filenames.get_note_string(note)

    fn = "".join([name, note_str, unique_str, length_str, ".tsv"])

    return os.path.join(seq_base, "count-tables", fn)


def get_sample_name_map(config, seq):

    sample_name_map = _return_key_dict()
    key = "riboseq_sample_name_map"
    if seq == "rna":
        key = "rnaseq_sample_name_map"
    if key in config:
        sample_name_map.update(config[key])

    return sample_name_map
