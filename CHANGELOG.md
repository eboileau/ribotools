# Change Log

All notable changes to the `ribotools` package will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/),
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Changed

- Requires rpbp >= 4.0.0
- Using 'ashr' instead of 'apeglm'

### Added

- Support for Python 3.11, 3.12, and 3.13

### Fixed

- 'coef' should specify same coefficient as in results 'res'

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
