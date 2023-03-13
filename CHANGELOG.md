# Change Log
All notable changes to the `ribotools` package will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

Change API to conform to Rp-Bp v3.

### Changed
- `pep_bed6_to_bed12.py` added option
- `create_rna_profiles.py`, `extract_rna_profiles.py` removed reference to btea, added defaults, using pbio 1.0.0
- script `create_bigBed_tracks.py` removed, ORF-specific functionalities added to `prep_orf_beds.py` (bed to bigBed format
conversion moved to `trackhub-utils`)

### Removed
- Isoform strategy (Fix reference to `btea`).

### Added
- options to `create_bigBed_tracks.py`
- notebooks
- match_appris_scores.py, get_psite_count_table.py, get_all_seq_filtering_counts.py (draft scripts only)
- genome to dge script
- `get_all_bam_periodic.py` to pgrms
- `create_bigBed_tracks.py` to pgrms
- MANIFEST.in, scripts (`create_trackDb`) and data (`fields.txt`)
- scripts dge/ and R script for DGEA from htseq count workflow
- `pep_bed6_to_bed12.py` to pgrms

## [0.1.0]

Restructure package, add to GitHub. Previous changes not documented.
