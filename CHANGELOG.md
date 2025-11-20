# Change Log

All notable changes to the `ribotools` package will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/),
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] - started 2025-11

## [2.0.0] 2025-11-20

### Changed

- Requires rpbp >= 4.0.1
- Environment (R=4.4.3, argparser)
- Option `--stranded` to `--rna-stranded`
- Scripts
- Documentation

### Added

- Support for Python 3.11, 3.12, and 3.13
- Regression testing and CI
- Option to create ORF profiles for QC
- Slurm option for `get_gtf_from_predictions.py`

### Fixed

- 'coef' should specify same coefficient as in results 'res'
- Allow both 'samples' and 'sample_name_map' in scripts

### Removed

- Support for Python 3.7, 3.8 (EOL), 3.9 and 3.10

## [1.0.3]

### Added

- TEA/DEA scripts
- PyPI deploy
- Documentation (RTD)

### Fixed

- Discouraged "script-files" (wrapped as module for packaging)
- Import

## [1.0.0]

### Changed

- API to conform to Rp-Bp v3

### Removed

- Track-related scripts (moved to `trackhub-utils`)
- `pep_bed6_to_bed12.py`
- `match_appris_scores.py`, `get_psite_count_table.py`, `get_all_seq_filtering_counts.py`, _etc._ (draft scripts)
- `create_rna_profiles.py`, `extract_rna_profiles.py`, _etc._
- Notebooks

## [0.1.0]

Restructure package, add to GitHub. Previous changes not documented.
