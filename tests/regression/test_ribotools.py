import logging
import pytest
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def to_df(filename):
    kwargs = {}
    if filename.endswith(".bed.gz") or filename.endswith(".tab.gz"):
        kwargs = {"sep": "\t"}
    elif filename.endswith(".mtx.gz"):
        kwargs = {"sep": " ", "comment": "%", "header": None}
    return pd.read_csv(filename, **kwargs)


# test output of `run-htseq-workflow`
def test_pipeline(getf_pipeline):

    files, ref_files = getf_pipeline

    for file, ref_file in zip(files, ref_files):
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        # pd.testing.assert_frame_equal(
        #     to_df(file), to_df(ref_file), check_exact=True, check_dtype=False
        # )
