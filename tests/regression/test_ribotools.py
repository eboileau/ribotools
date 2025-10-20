import logging
import pandas as pd

logger = logging.getLogger(__name__)


def to_df(filename):
    kwargs = {}
    names = ["id", "type", "gene", "count"]
    if filename.endswith(".tsv"):
        kwargs = {"sep": "\t", "header": None, "names": names}
    elif filename.endswith(".txt"):
        kwargs = {"sep": "\t", "header": None, "names": names[0:4:2]}
    return pd.read_csv(filename, **kwargs)


# test output of `run-htseq-workflow` and `get-sample-table`
def test_pipeline(getf_pipeline):

    files, ref_files = getf_pipeline

    for file, ref_file in zip(files, ref_files):
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        pd.testing.assert_frame_equal(
            to_df(file), to_df(ref_file), check_exact=True, check_dtype=False
        )
