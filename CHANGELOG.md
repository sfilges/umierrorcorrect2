# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.32.4] - 2026-02-21

### Added

- **Added downsampling** to the pipeline by default.
- **Updated CI workflows** for Docker to only run on tags.
- **Updated documeentation** to reflect current implementation.
- **Fixed numpy cast error** in `downsampling.py`.

## [0.32.3] - 2026-02-21

### Fixed

- **Docker: fastp installation**: Replaced unverified GitHub binary download of fastp with `apt-get install fastp` (Ubuntu 24.04 provides v0.23.4), ensuring package integrity via GPG-signed apt infrastructure and adding multi-architecture support.
- **CI: Docker workflow branch**: Corrected workflow trigger branch from `main` to `master` so the Docker publish workflow fires on pushes to the correct default branch.
- **CI: Docker workflow tags/labels**: Fixed incorrect `id.meta.outputs` expression (not valid GitHub Actions syntax) to `steps.meta.outputs`, ensuring the published image is correctly tagged and labelled.

## [0.32.2] - 2026-02-21

### Fixed

- **Unmapped reads crash**: Fixed a `TypeError` ("argument of type 'NoneType' is not a container or iterable") that occurred during consensus generation when unmapped mate-pair reads were fetched from the BAM file. Unmapped reads stored at their mapped mate's coordinates are now skipped in `get_cons_dict` (root cause) and `get_cons_info` guards against a `None` cigar string defensively. This crash was reproducible when running without fastp pre-filtering.

## [0.32.1] - 2026-02-02

### Added

- **Console Status**: Added console status updates for FastQC and MultiQC to provide better feedback during batch processing.

### Changed

- **Batch Processing**: Improved batch processing to show progress and status updates for each sample.
- **Mutation BED Template**: Added mutation BED template file for ctDNA tracking to reflect the test data.
- **VCF Output**: VCF output now contains all variants even if not bed file is provided to provide name annotations.

## [0.32.0] - 2026-01-30

### Added

- **Integration Testing**: Added end-to-end integration test `tests/integration/test_pipeline.py` covering the full pipeline execution (FASTQ -> BAM -> VCF).
- **Test Data**: Included slimmed-down real datasets in `tests/test_data/` for rapid integration testing.

### Changed

- **CLI Renamed**: The command-line tool has been renamed from `umierrorcorrect` to `umierrorcorrect2` to avoid conflicts with legacy installations.
- **Type Safety**: Re-enabled and fixed strict `mypy` type checking for core modules: `core.consensus`, `core.group`, and `core.umi_cluster`.
- **Renamed batch to run**: The `batch` command has been renamed to `run` to better reflect its purpose.

### Fixed

- **Testing Infrastructure**: Fixed the `test_data_dir` fixture path in `tests/conftest.py` to correctly locate test resources.

## [0.31.0] - 2026-01-31

### Changed

- **Made numba default** for calculating hamming distances; also included numba in the default dependencies.
- **Removed Python 3.14** from supported versions.
- **Added pre-commit hooks** for code quality checks.
- Updated README.md.

## [0.30.4] - 2026-01-27

### Changed

- **Target Coverage**: `calculate_target_coverage()` now reports both on-target and off-target read counts and fractions. Output columns: family_size, on_target, off_target, total, on_target_fraction, off_target_fraction.
- **Target Coverage Conditional**: The `_target_coverage.txt` file is now only generated when a BED file is provided, as the on/off-target distinction requires region annotations.
- **Moved get_percent_mapped_reads** to utils.py. The function is not used anywhere else in the codebase right now.

## [0.30.3] - 2026-01-26

### Changed

- **Stats Extraction Refactored**: Statistics are now derived entirely from the consensus BAM file, eliminating the need for separate `.stats` files during processing.
  - `get_stat()` now takes a BAM file and optional BED file (instead of BAM + stats file)
  - `run_get_consensus_statistics()` updated to use BAM-only approach
  - `run_downsampling()` updated to use BAM-only approach
  - CLI `stats` and `downsampling` commands: replaced `--stats` parameter with `--bed-file/-bed`
- **Stats Output**: The `_consensus_stats.tsv` file now reports statistics per unique read start position rather than per clustered region. This provides more granular output.
- Updated dependency versions in `pyproject.toml`

### Removed

- Removed intermediate stats file generation from `cluster_consensus_worker()`
- Removed `merge_tmp_stats_files()` function from `umi_error_correct.py`
- Removed `merge_duplicate_stat()` function from `umi_error_correct.py`

### Added

- New `RegionStats` dataclass for BAM-derived statistics
- New `parse_consensus_read_name()` function for parsing consensus read names
- New `get_stats_from_bam()` function for extracting stats directly from BAM
- New `write_stats_file()` function for backward-compatible stats file output
- New `RegionConsensusStats.from_region_stats()` factory method
- Added unit tests for new statistics extraction functions

## [0.30.2] - 2026-01-24

### Fixed

- **Grouping**: R1 and R2 were counted as seprate reads leading to inflated consensus coverage.
- `count_umis_in_region`: Now counts unique Query Names (QNAMEs) per UMI.
- `group_by_position`: Updated the default automatic grouping mode to deduplicate mates within each 20bp region.
- `read_bam_from_tag`: Updated the tag-based grouping mode to deduplicate mates sharing the same UG tag.
- Bounding Box Robustness: In read_bam_from_tag, I also ensured that the starts dictionary tracks the minimum start position of any read in the group, ensuring the consensus covers the entire fragment range.
- Added CI testing for Python 3.13

## [0.30.1] - 2026-01-23

### Added

- **Shell Completion**: Re-enabled `typer` shell completion support via `--install-completion`.

### Fixed

- **Numba Stability**: Added robust error handling for broken Numba/LLVM installations to ensure graceful fallback.
- **Recursive Discovery**: Fixed recursive FASTQ discovery in subdirectories for the `batch` command.
- **Type Safety**: Resolved 80+ `mypy` errors and improved type annotations across the codebase.
- **CI Stability**: Adjusted `mypy` and `ruff` configurations for reliable production CI/CD.
- **Installation**: Removed unnecessary `[all]` extra from `typer` dependency to avoid installation warnings.

## [0.30.0] - 2026-01-23

Initial release of umierrorcorrect2 (a fork of umierrorcorrect v0.29).

### Added

- **Modern CLI**: Completely redesigned CLI using `typer` for better user experience and documentation.
- **GitHub Actions CI**: Automated linting, type checking, and testing for Python 3.9-3.12.
- **Pydantic Models**: Added robust configuration and data validation using `pydantic`.
- **Fastp Integration**: Support for `fastp` for high-performance preprocessing and UMI extraction.
- **Improved Logging**: Integrated `loguru` for structured and informative logging.
- **Comprehensive Testing**: Added extensive unit and integration tests with `pytest`.
- **Containerization**: Added Docker support for easy deployment.
- **Documentation**: New `USER_GUIDE.md`, `IMPLEMENTATION.md`, and updated `README.md`.
- **Modern Tooling**: Switched to `hatch` for build management and `uv` for dependency resolution.

### Changed

- **Codebase Refactoring**: Major refactor of core modules (consensus, alignment, clustering) with PEP 8 compliance and type hints.
- **Dependency Management**: Migrated to `pyproject.toml` (PEP 621) and `uv.lock`.
- **Performance**: Optimized UMI clustering and consensus generation for parallel processing.
- **License**: Updated author to Stefan Filges and extended copyright year.

### Fixed

- Fixed various bugs in UMI extraction and merging logic.
- Resolved issues with single-end read processing in `fastp`.
- Fixed path handling issues using `pathlib`.

---
[0.32.1]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.32.1
[0.32.0]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.32.0
[0.31.0]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.31.0
[0.30.4]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.4
[0.30.3]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.3
[0.30.2]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.2
[0.30.1]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.1
[0.30.0]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.0
