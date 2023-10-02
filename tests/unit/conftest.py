"""
    confest.py
"""

import pytest


@pytest.fixture()
def get_bed():
    from io import StringIO

    string = """#seqname\t#start\t#end\t#id\t#score\t#strand\t#thick_start\t#thick_end\t#color\t#num_exons\t#exon_lengths\t#exon_genomic_relative_starts\t#condition\t#bayes_factor_mean\t#bayes_factor_var\t#x_1_sum\t#x_2_sum\t#x_3_sum\t#orf_num\t#orf_len\t#orf_type\t#biotype\t#transcripts\t#gene_id\t#gene_name\t#gene_biotype
    1\t634018\t634051\tENST1_1:634018-634051:+\t0\t+\t634018\t634051\t0\t1\t33\t0\tctrl\t5.0\t1.0\t10.0\t2.0\t1.0\t1\t30\tncORF\tunprocessed_pseudogene\tENST1\tENSG1\tGENE1\tunprocessed_pseudogene
    2\t841343\t852064\tENST2_2:841343-852064:+\t0\t-\t841343\t852064\t0\t2\t30,138\t0,10583\tcond\t50.0\t2.0\t20.0\t5.0\t2.0\t2\t168\tuORF\tprotein_coding\tENST2,ENST3\tENSG2\tGENE2\tprotein_coding
    2\t2189711\t2198988\tENST4_2:2189711-2198988:+\t0\t-\t2189711\t2198988\t0\t5\t70,272,136,136,235\t0,3927,4286,8295,9042\tcond\t100.0\t3.0\t45.0\t10.0\t5.0\t3\t849\tCDS\tprotein_coding\tENST4\tENSG3\tGENE3\tprotein_coding
    3\t2228766\t2306762\tENST5_3:2228766-2306762:+\t0\t+\t2228766\t2306762\t0\t7\t969,126,116,269,293,231,186 \t0,74211,74518,75073,75526,77253,77810\tctrl\t2500.0\t2.5\t120.0\t10.0\t1.0\t4\t2190\tNovel altCDS\t\tENST4\tENSG4\t\t"""
    return StringIO(string)


@pytest.fixture()
def get_gtf(tmp_path_factory, get_bed):
    """\
    Run `get-gtf-from-predictions`.
    """
    from pathlib import Path
    import pbiotools.utils.bed_utils as bed_utils
    import pbiotools.misc.shell_utils as shell_utils

    loc = tmp_path_factory.mktemp("data")
    outloc = Path(loc, "test.gtf").as_posix()
    bed = bed_utils.read_bed(get_bed)
    inloc = Path(loc, "test.bed").as_posix()
    bed_utils.write_bed(bed, inloc, compress=False)

    num_cpus = 6
    cmd = f"get-gtf-from-predictions {inloc} {outloc} --num-cpus {num_cpus}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return inloc, outloc
