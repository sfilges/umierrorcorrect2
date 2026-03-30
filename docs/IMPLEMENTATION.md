# Implementation Status & Architecture: UMIErrorCorrect

This document outlines the current architecture, utilized technologies, and implementation details of the UMIErrorCorrect package. It serves as a guide for developers to understand the codebase structure and design decisions.

## Code Guidelines

- **Parallel Processing**: Uses Python's `multiprocessing.Pool` for parallel region processing.
- **BAM/SAM Handling**: All alignment operations use the `pysam` library.
- **Region Detection**: Regions can be auto-detected from BAM or defined via BED file.
- **Package Management**: Orchestration using `uv`.
- **Build System**: PEP 621 compliant `pyproject.toml` with `hatchling` build backend.
- **CLI**: Unified CLI using `typer`.
- **Logging**: Structured logging with `loguru` and `rich`.
- **Linting & Formatting**: Enforced by `ruff`.
- **Validation**: Data validation using `pydantic`.
- **Type Hints**: Extensive usage of type annotations and `mypy` for static analysis.

---

## 1. Architecture Overview

The package has been modernized to use standard Python packaging and tooling.

### Directory Structure

```text
umierrorcorrect/
├── umierrorcorrect/        # Main package source
│   ├── core/               # Core libraries and utilities
│   │   ├── consensus.py    # Consensus generation logic
│   │   ├── constants.py    # Global constants and type aliases
│   │   ├── logging_config.py # Loguru configuration
│   │   ├── umi_cluster.py  # UMI clustering algorithms
│   │   └── ...
│   ├── models/             # Pydantic data models
│   ├── cli.py              # Main CLI entry point (Typer)
│   ├── batch.py            # Batch processing logic
│   ├── glue logic files... # (align.py, preprocess.py, etc.)
│   └── version.py
├── tests/                  # Test suite
│   ├── unit/
│   ├── integration/
│   └── conftest.py
├── pyproject.toml          # Project configuration
└── uv.lock                 # Dependency lock file
```

### Build & Dependency Management

- **Build Backend**: `hatchling`
- **Dependency Manager**: `uv`
- **Configuration**: All metadata and tool configuration (Ruff, Pytest, Coverage, Mypy) is consolidated in `pyproject.toml`.

### Command Line Interface

The package provides a single entry point `umierrorcorrect` which exposes multiple subcommands via `Typer`.

**Entry Point**: `umierrorcorrect.cli:main_cli`

**Subcommands**:

- `batch`: Run the complete pipeline (preprocessing -> mapping -> consensus -> stats -> variants).
- `preprocess`: UMI extraction and adapter trimming.
- `mapping`: BWA alignment.
- `consensus`: Generate consensus sequences (error correction).
- `stats`: Generate consensus statistics.
- `variants`: Call variants from consensus reads.
- `filter-bam`: Filter BAM files by consensus depth.
- `downsampling`: Analysis of downsampling.

---

## 2. Pipeline Implementation

The pipeline processes raw sequencing data into error-corrected consensus reads and variants.

### Pipeline Flow Diagrams

The preprocessing pipeline handles adapter trimming and UMI extraction.

#### Overview: Complete Pipeline

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                         UMI Error Correct Pipeline                          │
└─────────────────────────────────────────────────────────────────────────────┘

Input: FASTQ R1 (+ R2 for paired-end)
                    │
                    ▼
        ┌───────────────────────┐
        │   PREPROCESSING       │ ◄── Handles adapter trimming & UMI extraction
        │   (fastp or cutadapt) │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   BWA MAPPING         │
        │   (align_bwa)         │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   UMI ERROR CORRECT   │
        │   - UMI clustering    │
        │   - Consensus gen     │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   CONSENSUS STATS     │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   VARIANT CALLING     │
        └───────────────────────┘
                    │
                    ▼
Output: VCF, consensus BAM, statistics
```

#### Preprocessing Logic

Configuration flags determine whether `fastp` (default, faster) or `cutadapt` is used.

##### **Path A: With fastp (Default)**

- Quality filtering (Q20 default)
- Adapter trimming
- Read merging (optional)
- UMI extraction

##### **Path B: Without fastp (`--no-fastp`)**

- Adapter trimming via `cutadapt` (if enabled)
- No quality filtering
- UMI extraction

---

## 3. Core Components

### 3.1 UMI Clustering & Error Correction

**Location**: `umierrorcorrect/core/umi_cluster.py` & `umierrorcorrect/umi_error_correct.py`

- **Algorithm**: Adjacency-based clustering using Hamming or Levenshtein distance.
- **Optimization**: Accesses `pysam` objects directly and uses `multiprocessing` to handle genomic regions in parallel.
- **Security**: Subprocess calls (e.g., for sorting/uniq) use `shutil.which` and avoid `shell=True` to prevent command injection. Secure temporary files are used.

### 3.2 Consensus Generation

**Location**: `umierrorcorrect/core/consensus.py`

- Generates consensus sequences from clustered reads.
- Supports position-based grouping.
- Filters consensus bases based on frequency thresholds (default 60%).

### 3.3 Logging System

**Location**: `umierrorcorrect/core/logging_config.py`

- Implemented using `loguru`.
- **Console**: Uses `rich` for pretty formatting.
- **File**: Rotated, compressed log files saved to the output directory.
- **Context**: Captures function names and line numbers for debugging.

---

## 4. Development & Testing Infrastructure

### 4.1 Development Setup

Using `uv` for fast, reproducible environments.

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Sync dependencies
uv sync --dev

# Run tests
uv run pytest
```

### 4.2 Test Suite

**Location**: `tests/`

- **Framework**: `pytest`
- **Unit Tests**: `tests/unit/` - Test individual functions and classes (formatting, clustering logic, args).
- **Integration Tests**: `tests/integration/` - Test full pipeline steps (preprocessing, alignment).
- **Fixtures**: Shared test data and artifacts defined in `tests/conftest.py`.

### 4.3 Code Quality Tools

All configured in `pyproject.toml`.

- **Linting**: `ruff check .`
- **Formatting**: `ruff format .`
- **Type Checking**: `mypy .`
- **Coverage**: `pytest --cov=umierrorcorrect`

---

## 5. Completed Modernization Tasks

The following improvements have been fully implemented:

1. **Refactored to `src` layout (implicit)**: Code moved to `umierrorcorrect/` package structure.
2. **Modern Packaging**: Replaced `setup.py` with `pyproject.toml`.
3. **Typer CLI**: Replaced `argparse` with `Typer` for a robust, self-documenting CLI.
4. **Security Fixes**: Removed vulnerable `shell=True` usage in subprocess calls.
5. **Resource Management**: Improved temporary file handling with `tempfile` and `contextlib`.
6. **Type Safety**: Added type hints to core modules (`constants.py`, `umi_cluster.py`, etc.).
7. **Constants**: Externalized magic numbers to `core/constants.py`.
