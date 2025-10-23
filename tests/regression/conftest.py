import pytest
import sys

REF_DATA_URL = "https://data.dieterichlab.org/s/Rqi67d2RHs2zBJM/download"
REF_LOC = "hiPSC-CM-chr1-example"
REF_CONFIG = "hiPSC-CM-test.yaml"
REF_SAMPLE_TBL = "sample-table.csv"
REG_CLASSES = ["buffered.txt", "exclusive.txt", "forwarded.txt", "intensified.txt"]


@pytest.fixture(scope="session")
def data_loc(tmp_path_factory):
    """Download reference dataset for regression testing.

    Parameters
    ----------
    tmp_path_factory
        tmp_path_factory fixture

    Returns
    -------
    :pathlib.Path: base temporary directory
    """
    import pbiotools.misc.shell_utils as shell_utils

    loc = tmp_path_factory.mktemp("data")

    cmd = f"wget --no-verbose {REF_DATA_URL} -O {loc}/data.zip"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    cmd = f"unzip {loc}/data.zip -d {loc}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return loc


@pytest.fixture(scope="session")
def getf_config(data_loc):
    """Set configuration file.

    Parameters
    ----------
    data_loc
        data_loc fixture

    Returns
    -------
    :obj:`tuple`: configuration files
    """
    import yaml
    from pathlib import Path

    loc = Path(data_loc, REF_LOC)
    config = Path(loc, REF_CONFIG)
    config.write_text(
        config.read_text().replace("/path/to/your/hiPSC-CM-example", loc.as_posix())
    )
    sample_tbl = Path(loc, "reference", "tea-results", REF_SAMPLE_TBL)
    sample_tbl.write_text(
        sample_tbl.read_text().replace("/path/to/your/hiPSC-CM-example", loc.as_posix())
    )

    # reference (known) paths from example dataset, keep default options
    ref_config = yaml.load(open(config), Loader=yaml.FullLoader).copy()
    ref_config["riboseq_data"] = Path(loc, "reference", "riboseq-results").as_posix()
    ref_config["rnaseq_data"] = Path(loc, "reference", "rnaseq-results").as_posix()
    ref_config["tea_data"] = Path(loc, "reference", "tea-results").as_posix()
    ref_config["ribo_gtf"] = Path(loc, "input", "ribo-ORFs.chr1.gtf").as_posix()

    return (config, ref_config)


@pytest.fixture(scope="session")
def get_pipeline(getf_config):
    """Run `run-htseq-workflow`.

    Also run `get-sample-table`.

    Parameters
    ----------
    getf_config
        getf_config fixture

    Returns
    -------
    :obj:`tuple`: configuration files
    """
    import pbiotools.misc.shell_utils as shell_utils

    config, ref_config = getf_config

    if sys.platform == "darwin":
        num_cpus = 1  # avoid parallel processing issue on macos: https://github.com/dieterich-lab/rp-bp/issues/140
    else:
        num_cpus = 6  # multiprocessing.cpu_count() see https://github.com/dieterich-lab/rp-bp/issues/144
    opts = (
        "--run-all "
        "--trim-rna-to-max-fragment-size  "
        '--star-options "--quantMode GeneCounts" '
        '--htseq-options "--idattr orf_id" "--additional-attr orf_type" "--additional-attr gene_name" "--stranded yes" '
        "--rna-stranded reverse "
        f'--gtf {ref_config["ribo_gtf"]} '
        "--keep-intermediate-files "
    )
    cmd = (
        f"run-htseq-workflow ribo {config.as_posix()} "
        f"--ribo-config {config.as_posix()} "
        f"--rna-config {config.as_posix()} "
        f"--num-cpus {num_cpus} {opts}"
    )
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    cmd = f"get-sample-table {config.as_posix()}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    opts = "--orfCol 2 --symbolCol 3 --lfcThreshold 0 --alpha .99"
    cmd = f"run-tea {opts} {config.as_posix()}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return config, ref_config


@pytest.fixture(scope="session")
def getf_pipeline(get_pipeline):
    """Get the output file names.

    Count tables (with periodic lengths in file names)
    for the current output and the reference dataset,
    the sample table, and the txt files containing the
    TE features.

    Parameters
    ----------
    get_pipeline
        Fixture calling the pipeline

    Returns
    -------
    :obj:`list`: tuples of output files
    """
    import yaml
    from pathlib import Path
    import rpbp.ribo_utils.utils as ribo_utils
    import ribotools.utils.filenames as filenames

    from rpbp.defaults import metagene_options

    config, ref_config = get_pipeline
    config = yaml.load(open(config), Loader=yaml.FullLoader)

    # identical for config and ref_config
    sample_names = sorted(config["riboseq_samples"].keys())
    sample_names.extend(sorted(config["rnaseq_samples"].keys()))
    keys = [key for key in ["ribo", "rna"] for _ in range(len(sample_names) // 2)]

    lfiles = [[], []]
    lconfigs = [config, ref_config]

    def populate(name, lf, lc, key=None):
        note = lc.get("note", None)
        is_unique = f"keep_{key}seq_multimappers" not in lc
        sample_name_map = filenames.get_sample_name_map(lc, key)

        if key == "ribo":
            lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
                lc, name, is_unique=is_unique, default_params=metagene_options
            )
        else:
            is_unique_ribo = "keep_riboseq_multimappers" not in lc
            matching_ribo_sample = config["matching_samples"][name]
            lengths, _ = ribo_utils.get_periodic_lengths_and_offsets(
                lc,
                matching_ribo_sample,
                is_unique=is_unique_ribo,
                default_params=metagene_options,
            )
            lengths = str(max([int(pl) for pl in lengths]))

        count_file = filenames.get_count_table(
            lc[f"{key}seq_data"],
            sample_name_map[name],
            is_unique=is_unique,
            length=lengths,
            note=note,
        )
        lf.append(count_file)

    for name, key in zip(sample_names, keys):
        for lf, lc in zip(lfiles, lconfigs):
            populate(name, lf, lc, key=key)

    files = [file for file in lfiles[0]]
    ref_files = [file for file in lfiles[1]]

    files.append(Path(config["tea_data"], REF_SAMPLE_TBL).as_posix())
    ref_files.append(Path(ref_config["tea_data"], REF_SAMPLE_TBL).as_posix())

    contrast = [*config["contrasts"]][0]
    for reg_class in REG_CLASSES:
        files.append(Path(config["tea_data"], "LRT", contrast, reg_class).as_posix())
        ref_files.append(
            Path(ref_config["tea_data"], "LRT", "d5_vs_d1", reg_class).as_posix()
        )

    return (files, ref_files)
