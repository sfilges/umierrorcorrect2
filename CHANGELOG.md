# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.30.2] - 2026-01-24

### Fixed

- **Grouping**: R1 and R2 were counted as seprate reads leading to inflated consensus coverage.
- `count_umis_in_region`: Now counts unique Query Names (QNAMEs) per UMI.
- `group_by_position`: Updated the default automatic grouping mode to deduplicate mates within each 20bp region.
- `read_bam_from_tag`: Updated the tag-based grouping mode to deduplicate mates sharing the same UG tag.
- Bounding Box Robustness: In read_bam_from_tag, I also ensured that the starts dictionary tracks the minimum start position of any read in the group, ensuring the consensus covers the entire fragment range.

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
[0.30.1]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.1
[0.30.0]: https://github.com/sfilges/umierrorcorrect2/releases/tag/v0.30.0
