"""
    Unit tests using the hIPSC-CMs reference dataset.
"""

import logging
import pytest
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def _get_expected_gtf():
    from io import StringIO

    string = """1\tRp-Bp\tCDS\t634019\t634051\t0\t+\t.\tgene_id ENSG1; gene_name GENE1; gene_biotype unprocessed_pseudogene; transcript_id ENST1; transcript_biotype unprocessed_pseudogene; orf_type ncORF; orf_id ENST1_1:634018-634051:+
    1\tRp-Bp\texon\t634019\t634051\t0\t+\t.\tgene_id ENSG1; gene_name GENE1; gene_biotype unprocessed_pseudogene; transcript_id ENST1; transcript_biotype unprocessed_pseudogene; orf_type ncORF; orf_id ENST1_1:634018-634051:+
    2\tRp-Bp\texon\t841344\t841373\t0\t-\t.\tgene_id ENSG2; gene_name GENE2; gene_biotype protein_coding; transcript_id ENST2,ENST3; transcript_biotype protein_coding; orf_type uORF; orf_id ENST2_2:841343-852064:+
    2\tRp-Bp\tCDS\t841344\t841373\t0\t-\t.\tgene_id ENSG2; gene_name GENE2; gene_biotype protein_coding; transcript_id ENST2,ENST3; transcript_biotype protein_coding; orf_type uORF; orf_id ENST2_2:841343-852064:+
    2\tRp-Bp\texon\t851927\t852064\t0\t-\t.\tgene_id ENSG2; gene_name GENE2; gene_biotype protein_coding; transcript_id ENST2,ENST3; transcript_biotype protein_coding; orf_type uORF; orf_id ENST2_2:841343-852064:+
    2\tRp-Bp\tCDS\t851927\t852064\t0\t-\t.\tgene_id ENSG2; gene_name GENE2; gene_biotype protein_coding; transcript_id ENST2,ENST3; transcript_biotype protein_coding; orf_type uORF; orf_id ENST2_2:841343-852064:+
    2\tRp-Bp\texon\t2189712\t2189781\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\tCDS\t2189712\t2189781\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\texon\t2193639\t2193910\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\tCDS\t2193639\t2193910\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\texon\t2193998\t2194133\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\tCDS\t2193998\t2194133\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\texon\t2198007\t2198142\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\tCDS\t2198007\t2198142\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\texon\t2198754\t2198988\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    2\tRp-Bp\tCDS\t2198754\t2198988\t0\t-\t.\tgene_id ENSG3; gene_name GENE3; gene_biotype protein_coding; transcript_id ENST4; transcript_biotype protein_coding; orf_type CDS; orf_id ENST4_2:2189711-2198988:+
    3\tRp-Bp\texon\t2228767\t2229735\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2228767\t2229735\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\texon\t2302978\t2303103\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2302978\t2303103\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\texon\t2303285\t2303400\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2303285\t2303400\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\texon\t2303840\t2304108\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2303840\t2304108\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\texon\t2304293\t2304585\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2304293\t2304585\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\texon\t2306020\t2306250\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2306020\t2306250\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\texon\t2306577\t2306762\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+
    3\tRp-Bp\tCDS\t2306577\t2306762\t0\t+\t.\tgene_id ENSG4; transcript_id ENST4; orf_type 'Novel altCDS'; orf_id ENST5_3:2228766-2306762:+"""
    return StringIO(string)


def test_get_gtf(get_gtf):

    import pbiotools.utils.gtf_utils as gtf_utils

    _, gtf = get_gtf
    # type conversion - gtf_utils
    gtf = gtf_utils.read_gtf(gtf)
    expected_gtf = gtf_utils.read_gtf(_get_expected_gtf())

    # compare only on exon (or CDS), otherwise order might be wrong
    for feature in ["exon", "CDS"]:
        expected = expected_gtf[expected_gtf.feature == feature].copy()
        expected.reset_index(drop=True, inplace=True)
        returned = gtf[gtf.feature == feature].copy()
        returned.reset_index(drop=True, inplace=True)
        pd.testing.assert_frame_equal(returned, expected, check_exact=True)
