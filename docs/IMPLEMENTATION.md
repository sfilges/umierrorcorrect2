# Implementation Plan: UMIErrorCorrect Modernization

This document outlines a comprehensive plan to modernize the UMIErrorCorrect codebase, improve code quality, and establish proper testing infrastructure.

## Code guidelines

- Uses Python's `multiprocessing.Pool` for parallel region processing
- All BAM/SAM handling via `pysam` library
- Regions can be auto-detected from BAM or defined via BED file
- Version stored in `umierrorcorrect/version.py`
- Orchestration using `uv`
- Linting formatting with `ruff`
- Use `ruff check . --fix` to automatically fix linting errors
- Use `pydantic` for data validation
- Use `Path` instead of `os` for file handling
- Use type hints

## Phase 1: Modern Build System and Tooling (Foundation)

### 1.1 Create pyproject.toml with Hatchling

Replace `setup.py` and `setup.cfg` with a modern `pyproject.toml`:

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "umierrorcorrect"
dynamic = ["version"]
description = "Pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI)"
readme = "README.md"
license = "MIT"
requires-python = ">=3.9"
authors = [
    { name = "Tobias Osterlund", email = "tobias.osterlund@gu.se" }
]
classifiers = [
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
]
dependencies = [
    "pysam>=0.8.4",
    "scipy",
    "matplotlib",
]

[project.optional-dependencies]
dev = [
    "ruff",
    "pytest",
    "pytest-cov",
    "mypy",
]

[project.scripts]
run_umierrorcorrect = "umierrorcorrect.run_umierrorcorrect:main_cli"
preprocess = "umierrorcorrect.preprocess:main_cli"
run_mapping = "umierrorcorrect.run_mapping:main_cli"
umi_error_correct = "umierrorcorrect.umi_error_correct:main_cli"
get_consensus_statistics = "umierrorcorrect.get_consensus_statistics:main_cli"
call_variants = "umierrorcorrect.call_variants:main_cli"
filter_bam = "umierrorcorrect.filter_bam:main_cli"
filter_cons = "umierrorcorrect.filter_cons:main_cli"
downsampling_plots = "umierrorcorrect.downsampling_plots:main_cli"
fit_background_model = "umierrorcorrect.fit_background_model:main_cli"

[project.urls]
Homepage = "https://github.com/stahlberggroup/umierrorcorrect"
Documentation = "https://github.com/stahlberggroup/umierrorcorrect/wiki"
Repository = "https://github.com/stahlberggroup/umierrorcorrect"

[tool.hatch.version]
path = "umierrorcorrect/version.py"

[tool.hatch.build.targets.sdist]
include = [
    "/umierrorcorrect",
    "/test_data",
    "/doc",
]

[tool.hatch.build.targets.wheel]
packages = ["umierrorcorrect"]
```

### 1.2 Configure Ruff for Linting and Formatting

Add to `pyproject.toml`:

```toml
[tool.ruff]
target-version = "py39"
line-length = 120
src = ["umierrorcorrect"]

[tool.ruff.lint]
select = [
    "E",      # pycodestyle errors
    "W",      # pycodestyle warnings
    "F",      # Pyflakes
    "I",      # isort
    "B",      # flake8-bugbear
    "C4",     # flake8-comprehensions
    "UP",     # pyupgrade
    "ARG",    # flake8-unused-arguments
    "SIM",    # flake8-simplify
    "S",      # flake8-bandit (security)
    "PTH",    # flake8-use-pathlib
]
ignore = [
    "E501",   # line too long (handled by formatter)
    "S101",   # assert usage (ok in tests)
]
# Always auto-fix when possible
fixable = ["ALL"]

[tool.ruff.lint.per-file-ignores]
"tests/*" = ["S101", "ARG"]

[tool.ruff.format]
quote-style = "double"
indent-style = "space"
```

### 1.3 Set Up uv for Package Management

Create a `uv.lock` file and development workflow:

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Initialize project with uv
uv sync

# Development installation
uv sync --dev

# Run tools via uv
uv run ruff check --fix .
uv run ruff format .
uv run pytest
```

Add to documentation:

```markdown
## Development Setup

1. Install [uv](https://docs.astral.sh/uv/):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. Clone and install:

   ```bash
   git clone https://github.com/stahlberggroup/umierrorcorrect.git
   cd umierrorcorrect
   uv sync --dev
   ```

3. Run tests:

   ```bash
   uv run pytest
   ```
```

### 1.4 Migration Tasks

| Task | Description | Breaking Change |
|------|-------------|-----------------|
| Refactor entry points | Change from `scripts` to `project.scripts` with proper `main_cli()` functions | No (same CLI names) |
| Remove setup.py | Delete after pyproject.toml is validated | No |
| Remove setup.cfg | Delete after pyproject.toml is validated | No |
| Update version.py | Ensure compatible with hatch version reading | No |

**Entry point refactoring example** (for each script):

```python
# Current (run_umierrorcorrect.py):
if __name__ == '__main__':
    args = parseArgs()
    main(args)

# New:
def main_cli():
    """CLI entry point."""
    args = parseArgs()
    main(args)

if __name__ == '__main__':
    main_cli()
```

---

## Phase 2: Test Suite Infrastructure

### 2.1 Test Directory Structure

```text
tests/
├── __init__.py
├── conftest.py              # Shared fixtures
├── test_data/               # Small test fixtures (symlink or copy)
│   ├── small_R1.fastq.gz    # Minimal FASTQ for unit tests
│   ├── small.bam            # Minimal BAM for unit tests
│   └── regions.bed
├── unit/
│   ├── __init__.py
│   ├── test_umi_cluster.py
│   ├── test_get_consensus.py
│   ├── test_get_cons_info.py
│   ├── test_group.py
│   ├── test_get_regions_from_bed.py
│   ├── test_check_args.py
│   └── test_read_fastq_records.py
└── integration/
    ├── __init__.py
    ├── test_preprocess.py
    ├── test_run_mapping.py
    ├── test_umi_error_correct.py
    ├── test_pipeline.py      # End-to-end pipeline test
    └── test_call_variants.py
```

### 2.2 Pytest Configuration

Add to `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_functions = ["test_*"]
addopts = [
    "-v",
    "--tb=short",
    "--strict-markers",
]
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    "integration: marks tests as integration tests",
    "requires_bwa: marks tests that require bwa installed",
]
filterwarnings = [
    "ignore::DeprecationWarning:pysam.*",
]
```

### 2.3 Test Fixtures (conftest.py)

```python
import pytest
from pathlib import Path
import tempfile
import shutil

@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "test_data"

@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)

@pytest.fixture
def sample_umi_groups():
    """Sample UMI barcode groups for clustering tests."""
    return {
        "ACGTACGTACGT": ["read1", "read2", "read3"],
        "ACGTACGTACGA": ["read4"],  # 1 edit distance from first
        "GGGGGGGGGGGG": ["read5", "read6"],
    }

@pytest.fixture
def sample_reads():
    """Sample aligned reads for consensus testing."""
    # Return mock pysam AlignedSegment objects or dict representations
    pass
```

### 2.4 Priority Test Cases

**Unit tests for core algorithms (highest priority):**

| Module | Function | Test Cases |
|--------|----------|------------|
| `umi_cluster.py` | `hamming_distance()` | Identical strings, single diff, all diff, different lengths |
| `umi_cluster.py` | `edit_distance()` | Same as above plus insertions/deletions |
| `umi_cluster.py` | `cluster_barcodes()` | Various distance thresholds, edge cases |
| `umi_cluster.py` | `get_connected_components()` | Graph connectivity scenarios |
| `consensus.py` | `getConsensus3()` | Uniform reads, mixed reads, low quality, indels |
| `get_cons_info.py` | `get_cons_info()` | Various coverage depths, variant positions |
| `get_regions_from_bed.py` | `read_bed()` | Valid BED, empty, malformed |
| `get_regions_from_bed.py` | `merge_regions()` | Overlapping, adjacent, disjoint regions |

**Integration tests (medium priority):**

| Test | Description |
|------|-------------|
| `test_preprocess.py` | FASTQ preprocessing with various UMI configurations |
| `test_pipeline.py` | End-to-end with test_data files |

### 2.5 Coverage Configuration

Add to `pyproject.toml`:

```toml
[tool.coverage.run]
source = ["umierrorcorrect"]
branch = true
omit = [
    "*/test_*.py",
    "*/__init__.py",
]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "if __name__ == .__main__.:",
    "raise NotImplementedError",
]
show_missing = true
fail_under = 60  # Start low, increase over time
```

---

## Phase 3: Critical Bug Fixes and Security

### 3.1 CRITICAL: Fix Command Injection Vulnerability

**Location:** `umierrorcorrect/umi_error_correct.py` lines 227-229

**Current (vulnerable):**
```python
command1 = ['sort tmp.txt | uniq -d']
p1 = subprocess.Popen(command1, shell=True, stdout=g)
```

**Fixed:**
```python
import shutil

# Use proper subprocess chaining without shell=True
sort_cmd = [shutil.which("sort") or "sort", tmp_file.name]
uniq_cmd = [shutil.which("uniq") or "uniq", "-d"]

p1 = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE)
p2 = subprocess.Popen(uniq_cmd, stdin=p1.stdout, stdout=g)
p1.stdout.close()
p2.communicate()
```

### 3.2 Fix Resource Leak

**Location:** `umierrorcorrect/src/check_args.py` line 37-38

**Current:**

```python
devnull = open(os.devnull)
subprocess.Popen([name,'--version'], stdout=devnull, stderr=devnull)
```

**Fixed:**

```python
with open(os.devnull, 'w') as devnull:
    subprocess.run([name, '--version'], stdout=devnull, stderr=devnull, check=False)
```

### 3.3 Fix Unreachable Code

**Location:** `umierrorcorrect/src/check_args.py` lines 77-83

**Current:**

```python
except ValueError as e:
    raise(e + " Barcode length needs to be an integer")
    sys.exit(1)  # Unreachable
```

**Fixed:**

```python
except ValueError as e:
    raise ValueError(f"Barcode length must be an integer: {e}") from e
```

### 3.4 Use Temporary Files Properly

**Location:** `umierrorcorrect/umi_error_correct.py` lines 219, 228

**Current:**

```python
f = open('tmp.txt', 'w')
```

**Fixed:**

```python
import tempfile

with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
    tmp_file = f.name
    # ... use f ...
# Clean up later
os.unlink(tmp_file)
```

---

## Phase 4: Code Quality Improvements

### 4.1 Add Type Hints

**Priority files** (core algorithms):

```python
# umierrorcorrect/src/umi_cluster.py

from typing import Dict, List, Set, Tuple

def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings of equal length."""
    ...

def edit_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein edit distance between two strings."""
    ...

def cluster_barcodes(
    barcode_counts: Dict[str, int],
    distance_threshold: int = 1
) -> Dict[str, List[str]]:
    """Cluster barcodes by edit distance."""
    ...
```

Add mypy configuration to `pyproject.toml`:

```toml
[tool.mypy]
python_version = "3.9"
warn_return_any = true
warn_unused_configs = true
ignore_missing_imports = true
exclude = ["tests/", "build/"]

[[tool.mypy.overrides]]
module = "pysam.*"
ignore_missing_imports = true
```

### 4.2 Extract Constants

Create `umierrorcorrect/constants.py`:

```python
"""Constants used throughout the UMIErrorCorrect pipeline."""

# Default family sizes for consensus statistics
DEFAULT_FAMILY_SIZES: list[int] = [0, 1, 2, 3, 4, 5, 7, 10, 20, 30]

# Maximum reads per chunk for parallel processing
MAX_READS_PER_CHUNK: int = 100_000

# Default edit distance threshold for UMI clustering
DEFAULT_EDIT_DISTANCE_THRESHOLD: int = 1

# Default position threshold for grouping
DEFAULT_POSITION_THRESHOLD: int = 10

# Default consensus frequency threshold (percent)
DEFAULT_CONSENSUS_FREQUENCY: float = 60.0

# Default indel frequency threshold (percent)
DEFAULT_INDEL_FREQUENCY: float = 60.0
```

### 4.3 Refactor Long Functions

**Target: `consensus.py:getConsensus3()` (170 lines)**

Split into:
- `_initialize_consensus_read()` - Set up initial state
- `_process_alignment_column()` - Handle single column
- `_apply_quality_thresholds()` - Filter by quality
- `_build_cigar_string()` - Construct CIGAR
- `getConsensus3()` - Orchestrate the above

**Target: `umi_error_correct.py:cluster_consensus_worker()` (70 lines)**

Split into:
- `_cluster_umis()` - Just the clustering logic
- `_generate_consensus_for_cluster()` - Consensus generation
- `_write_cluster_output()` - File I/O
- `cluster_consensus_worker()` - Orchestrate

### 4.4 Modernize Path Handling

Replace string concatenation with pathlib:

```python
# Before
output_file = output_path + '/' + sample_name + '_consensus_reads.bam'
basename = filename.split('/')[-1]

# After
from pathlib import Path

output_file = Path(output_path) / f"{sample_name}_consensus_reads.bam"
basename = Path(filename).name
```

### 4.5 Modernize String Formatting

Replace `.format()` with f-strings:

```python
# Before
parser = argparse.ArgumentParser(description="UmiErrorCorrect v. {}. Pipeline...".format(__version__))

# After
parser = argparse.ArgumentParser(description=f"UmiErrorCorrect v. {__version__}. Pipeline...")
```

### 4.6 Remove Python 2 Compatibility Code

Remove from all files:
```python
from __future__ import division  # Not needed in Python 3
```

### 4.7 Fix Typos in User-Facing Text

| File | Line | Current | Fixed |
|------|------|---------|-------|
| `run_umierrorcorrect.py` | 44 | "emove the original" | "Remove the original" |
| `umi_error_correct.py` | 54 | "emove the original" | "Remove the original" |
| `run_mapping.py` | 26 | "emove the original" | "Remove the original" |
| `umi_error_correct.py` | 491 | "0cluster umis" | "cluster UMIs" |

### 4.8 Refactor CLI with Typer + Rich

Consolidate all CLI argument parsing into a single `cli.py` module using [Typer](https://typer.tiangolo.com/) and [Rich](https://rich.readthedocs.io/) for a modern, user-friendly CLI experience.

**Current state:** Each of the 10 CLI scripts uses argparse with mixed argument parsing and business logic.

**Target structure:**

```
umierrorcorrect/
├── cli.py              # Typer app with subcommands, Rich console output
├── preprocess.py       # preprocess() function (no CLI code)
├── run_mapping.py      # run_mapping() function
├── umi_error_correct.py # umi_error_correct() function
├── filter_bam.py       # filter_bam() function
└── ...
```

**Benefits:**
- **Testability:** Core functions can be unit tested without mocking CLI
- **Programmatic use:** Functions importable from Jupyter notebooks or custom pipelines
- **Single entry point:** `umierrorcorrect <command>` instead of 10 separate scripts
- **Beautiful output:** Rich progress bars, colored output, formatted tables
- **Type-driven:** Arguments derived from type hints, automatic validation
- **Auto-generated help:** Beautiful help text with examples

**Add dependencies to pyproject.toml:**

```toml
dependencies = [
    "pysam>=0.8.4",
    "scipy",
    "matplotlib",
    "typer[all]>=0.9.0",  # Includes Rich
]
```

**Example refactor for `filter_bam.py`:**

```python
# filter_bam.py (after refactor - core function only)
from pathlib import Path
import pysam

def filter_bam(infile: Path, outfile: Path, consensus_cutoff: int = 3) -> int:
    """Filter BAM file by consensus depth cutoff.

    Args:
        infile: Path to input BAM file.
        outfile: Path to output BAM file.
        consensus_cutoff: Minimum consensus depth to retain read.

    Returns:
        Number of reads written.
    """
    count = 0
    with pysam.AlignmentFile(infile, "rb") as f, \
         pysam.AlignmentFile(outfile, "wb", template=f) as g:
        for read in f.fetch():
            size = int(read.qname.rsplit("=", 1)[-1])
            if size >= consensus_cutoff:
                g.write(read)
                count += 1
    return count
```

```python
# cli.py (new file - Typer CLI)
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from umierrorcorrect import __version__
from umierrorcorrect.filter_bam import filter_bam
from umierrorcorrect.preprocess import preprocess

app = typer.Typer(
    name="umierrorcorrect",
    help="Pipeline for analyzing barcoded amplicon sequencing data with UMIs.",
    add_completion=False,
)
console = Console()


def version_callback(value: bool):
    if value:
        console.print(f"[bold blue]UMIErrorCorrect[/] v{__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        Optional[bool],
        typer.Option("--version", "-v", callback=version_callback, is_eager=True),
    ] = None,
):
    """UMIErrorCorrect: UMI-based error correction for amplicon sequencing."""
    pass


@app.command()
def filter(
    infile: Annotated[Path, typer.Option("--infile", "-i", help="Input BAM file")],
    outfile: Annotated[Path, typer.Option("--outfile", "-o", help="Output BAM file")],
    consensus_cutoff: Annotated[
        int, typer.Option("--consensus-cutoff", "-c", help="Minimum consensus depth")
    ] = 3,
):
    """Filter BAM file by consensus depth."""
    if not infile.exists():
        console.print(f"[red]Error:[/] Input file not found: {infile}")
        raise typer.Exit(1)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        progress.add_task("Filtering BAM...", total=None)
        count = filter_bam(infile, outfile, consensus_cutoff)

    console.print(f"[green]✓[/] Wrote {count:,} reads to {outfile}")


@app.command()
def preprocess(
    infile1: Annotated[Path, typer.Option("--read1", "-r1", help="R1 FASTQ file")],
    infile2: Annotated[Path, typer.Option("--read2", "-r2", help="R2 FASTQ file")],
    outdir: Annotated[Path, typer.Option("--output-dir", "-o", help="Output directory")],
    umi_length: Annotated[int, typer.Option("--umi-length", "-ul", help="UMI barcode length")] = 12,
    spacer_length: Annotated[int, typer.Option("--spacer-length", "-sl", help="Spacer length")] = 0,
):
    """Preprocess FASTQ files: extract UMIs and merge reads."""
    # ... implementation
    pass


@app.command()
def run(
    infile1: Annotated[Path, typer.Option("--read1", "-r1", help="R1 FASTQ file")],
    infile2: Annotated[Path, typer.Option("--read2", "-r2", help="R2 FASTQ file")],
    outdir: Annotated[Path, typer.Option("--output-dir", "-o", help="Output directory")],
    reference: Annotated[Path, typer.Option("--reference", "-r", help="Reference genome")],
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads")] = 1,
):
    """Run the complete UMI error correction pipeline."""
    console.print("[bold]UMIErrorCorrect Pipeline[/]")
    console.print()

    with Progress(console=console) as progress:
        task = progress.add_task("Running pipeline...", total=5)

        # Step 1: Preprocess
        progress.update(task, description="[cyan]Preprocessing reads...")
        # preprocess(...)
        progress.advance(task)

        # Step 2: Mapping
        progress.update(task, description="[cyan]Mapping to reference...")
        # run_mapping(...)
        progress.advance(task)

        # ... etc

    console.print("[green]✓ Pipeline complete![/]")


if __name__ == "__main__":
    app()
```

**pyproject.toml update:**

```toml
[project.scripts]
umierrorcorrect = "umierrorcorrect.cli:app"
# Keep legacy entry points for backwards compatibility (optional)
filter_bam = "umierrorcorrect.cli:filter"
preprocess = "umierrorcorrect.cli:preprocess"
```

**CLI output examples:**

```
$ umierrorcorrect --help

 Usage: umierrorcorrect [OPTIONS] COMMAND [ARGS]...

 UMIErrorCorrect: UMI-based error correction for amplicon sequencing.

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --version  -v        Show version and exit.                                  │
│ --help               Show this message and exit.                             │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ───────────────────────────────────────────────────────────────────╮
│ filter       Filter BAM file by consensus depth.                             │
│ preprocess   Preprocess FASTQ files: extract UMIs and merge reads.           │
│ run          Run the complete UMI error correction pipeline.                 │
│ ...                                                                          │
╰──────────────────────────────────────────────────────────────────────────────╯

$ umierrorcorrect filter -i input.bam -o output.bam -c 5
⠋ Filtering BAM...
✓ Wrote 1,234,567 reads to output.bam
```

**Migration steps:**

1. Add `typer[all]` to dependencies in pyproject.toml
2. Create `cli.py` with Typer app and subcommands
3. Move argument definitions from each script to `cli.py` as Typer options
4. Ensure core functions have clean signatures (Path objects, no argparse)
5. Add Rich progress bars and formatted output where appropriate
6. Update `pyproject.toml` entry points
7. Add deprecation warnings to old entry points (optional)

### 4.9 Comprehensive Logging with Loguru

Replace ad-hoc print statements and basic logging with [Loguru](https://github.com/Delgan/loguru) for structured, configurable logging throughout the pipeline.

**Why Loguru over stdlib logging:**
- Zero boilerplate - no handlers, formatters, or config files
- Automatic exception formatting with full traceback
- Built-in rotation, retention, and compression
- Lazy evaluation of log messages (no f-string overhead when disabled)
- Easy integration with Rich for colored console output
- Structured logging (JSON) for production environments

**Add dependency to pyproject.toml:**

```toml
dependencies = [
    "pysam>=0.8.4",
    "scipy",
    "matplotlib",
    "typer[all]>=0.9.0",
    "loguru>=0.7.0",
]
```

**Create logging configuration module:**

```python
# umierrorcorrect/logging_config.py
import sys
from pathlib import Path
from typing import Optional

from loguru import logger

# Remove default handler
logger.remove()

def setup_logging(
    verbosity: int = 0,
    log_file: Optional[Path] = None,
    json_logs: bool = False,
) -> None:
    """Configure logging for the pipeline.

    Args:
        verbosity: 0=WARNING, 1=INFO, 2=DEBUG
        log_file: Optional path to log file
        json_logs: If True, output structured JSON logs
    """
    # Map verbosity to log level
    level_map = {0: "WARNING", 1: "INFO", 2: "DEBUG"}
    level = level_map.get(verbosity, "DEBUG")

    # Console handler with Rich-compatible formatting
    console_format = (
        "<level>{level: <8}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
        "<level>{message}</level>"
    )

    # Simpler format for normal users (verbosity=0)
    simple_format = "<level>{level: <8}</level> | <level>{message}</level>"

    logger.add(
        sys.stderr,
        format=simple_format if verbosity == 0 else console_format,
        level=level,
        colorize=True,
    )

    # File handler (if specified)
    if log_file:
        file_format = (
            "{time:YYYY-MM-DD HH:mm:ss.SSS} | {level: <8} | "
            "{name}:{function}:{line} - {message}"
        )

        if json_logs:
            logger.add(
                log_file,
                format="{message}",
                level="DEBUG",
                serialize=True,  # JSON output
                rotation="10 MB",
                retention="1 week",
                compression="gz",
            )
        else:
            logger.add(
                log_file,
                format=file_format,
                level="DEBUG",
                rotation="10 MB",
                retention="1 week",
                compression="gz",
            )


def get_logger(name: str):
    """Get a logger instance bound with module name.

    Args:
        name: Module name (typically __name__)

    Returns:
        Bound logger instance
    """
    return logger.bind(name=name)
```

**Integrate with CLI:**

```python
# cli.py
from pathlib import Path
from typing import Annotated, Optional

import typer
from loguru import logger

from umierrorcorrect.logging_config import setup_logging

app = typer.Typer()


@app.callback()
def main(
    verbose: Annotated[
        int,
        typer.Option(
            "--verbose", "-v",
            count=True,
            help="Increase verbosity (-v=INFO, -vv=DEBUG)"
        ),
    ] = 0,
    log_file: Annotated[
        Optional[Path],
        typer.Option("--log-file", "-l", help="Write logs to file"),
    ] = None,
    json_logs: Annotated[
        bool,
        typer.Option("--json-logs", help="Output structured JSON logs"),
    ] = False,
):
    """UMIErrorCorrect: UMI-based error correction for amplicon sequencing."""
    setup_logging(verbosity=verbose, log_file=log_file, json_logs=json_logs)


@app.command()
def run(
    infile1: Annotated[Path, typer.Option("--read1", "-r1")],
    # ... other args
):
    """Run the complete pipeline."""
    logger.info("Starting UMIErrorCorrect pipeline")
    logger.debug(f"Input file: {infile1}")

    try:
        # ... pipeline steps
        logger.info("Preprocessing reads")
        # preprocess(...)

        logger.info("Mapping to reference")
        # run_mapping(...)

        logger.success("Pipeline completed successfully")

    except Exception as e:
        logger.exception("Pipeline failed")
        raise typer.Exit(1)
```

**Usage in core modules:**

```python
# umierrorcorrect/src/consensus.py
from umierrorcorrect.logging_config import get_logger

logger = get_logger(__name__)


def getConsensus3(reads, contig, regionid, ...):
    """Generate consensus sequence from reads."""
    logger.debug(f"Processing region {contig}:{regionid} with {len(reads)} reads")

    # ... processing logic

    if read_count < min_reads:
        logger.warning(
            f"Low read count for region {regionid}: {read_count} < {min_reads}"
        )

    logger.trace(f"Consensus quality scores: {qual_scores}")  # Very verbose

    return consensus_read
```

```python
# umierrorcorrect/umi_error_correct.py
from umierrorcorrect.logging_config import get_logger

logger = get_logger(__name__)


def cluster_consensus_worker(args):
    """Worker function for parallel consensus generation."""
    contig, pos, umis, bamfile, ... = args

    logger.debug(f"Worker processing {contig}:{pos} with {len(umis)} UMIs")

    with logger.contextualize(region=f"{contig}:{pos}"):
        # All logs in this block include region context
        for umi, reads in umis.items():
            logger.trace(f"Clustering UMI {umi} with {len(reads)} reads")
            # ... clustering logic

    logger.debug(f"Worker completed {contig}:{pos}")
    return results
```

**Log output examples:**

```bash
# Default (WARNING only)
$ umierrorcorrect run -r1 reads_R1.fq.gz -r2 reads_R2.fq.gz -o output/
WARNING  | Low read count for region chr1:12345: 2 < 3

# Verbose (-v = INFO)
$ umierrorcorrect run -v -r1 reads_R1.fq.gz ...
INFO     | Starting UMIErrorCorrect pipeline
INFO     | Preprocessing reads
INFO     | Mapping to reference
INFO     | Running UMI error correction
WARNING  | Low read count for region chr1:12345: 2 < 3
INFO     | Calling variants
SUCCESS  | Pipeline completed successfully

# Very verbose (-vv = DEBUG)
$ umierrorcorrect run -vv -r1 reads_R1.fq.gz ...
DEBUG    | consensus:getConsensus3:45 - Processing region chr1:12345 with 150 reads
DEBUG    | umi_error_correct:cluster_consensus_worker:89 - Worker processing chr1:12345 with 23 UMIs
...

# With log file (captures everything)
$ umierrorcorrect run -v --log-file pipeline.log -r1 reads_R1.fq.gz ...

# JSON logs for parsing/monitoring
$ umierrorcorrect run --json-logs --log-file pipeline.jsonl -r1 reads_R1.fq.gz ...
```

**JSON log output (for monitoring/parsing):**

```json
{"text": "Starting UMIErrorCorrect pipeline", "record": {"elapsed": {"repr": "0:00:00.001", "seconds": 0.001}, "level": {"name": "INFO", "no": 20}, "name": "umierrorcorrect.cli", "time": {"repr": "2024-01-15 10:30:00.123", "timestamp": 1705312200.123}}}
{"text": "Processing region chr1:12345 with 150 reads", "record": {"level": {"name": "DEBUG", "no": 10}, "extra": {"region": "chr1:12345"}, ...}}
```

**Timing and performance logging:**

```python
from loguru import logger
from time import perf_counter
from contextlib import contextmanager


@contextmanager
def log_duration(description: str):
    """Context manager to log operation duration."""
    start = perf_counter()
    logger.info(f"Starting: {description}")
    try:
        yield
    finally:
        elapsed = perf_counter() - start
        logger.info(f"Completed: {description} ({elapsed:.2f}s)")


# Usage:
with log_duration("consensus generation for chr1"):
    consensus = generate_consensus(reads)
```

**Integration with multiprocessing:**

```python
# For multiprocessing, configure logging in each worker
from loguru import logger
import multiprocessing


def worker_init(log_queue):
    """Initialize logging in worker process."""
    logger.remove()
    logger.add(log_queue.put, format="{message}", level="DEBUG")


def run_parallel(regions, num_workers):
    """Run parallel processing with centralized logging."""
    log_queue = multiprocessing.Queue()

    with multiprocessing.Pool(
        num_workers,
        initializer=worker_init,
        initargs=(log_queue,)
    ) as pool:
        # ... parallel work

    # Drain log queue in main process
    while not log_queue.empty():
        logger.info(log_queue.get())
```

**Migration steps:**

1. Add `loguru>=0.7.0` to dependencies
2. Create `umierrorcorrect/logging_config.py`
3. Add `--verbose`, `--log-file`, `--json-logs` options to CLI
4. Replace `print()` statements with appropriate log levels:
   - `print(f"Processing {x}")` → `logger.info(f"Processing {x}")`
   - `print(f"Warning: {x}")` → `logger.warning(x)`
   - Debug/trace output → `logger.debug()` or `logger.trace()`
5. Add `logger.exception()` in exception handlers
6. Use `logger.contextualize()` for region/UMI context in workers
7. Add timing logs for performance monitoring

---

## Phase 5: Pre-commit Hooks

Create `.pre-commit-config.yaml`:

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.8.0
    hooks:
      - id: ruff
        args: [--fix]
      - id: ruff-format

  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v5.0.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files
```

---

## Phase 6: Documentation Updates

### 6.1 Add Docstrings

Use Google-style docstrings:

```python
def cluster_barcodes(
    barcode_counts: Dict[str, int],
    distance_threshold: int = 1
) -> Dict[str, List[str]]:
    """Cluster UMI barcodes by edit distance.

    Groups barcodes that are within the specified edit distance
    threshold, with the most abundant barcode as the cluster center.

    Args:
        barcode_counts: Dictionary mapping barcode sequences to their counts.
        distance_threshold: Maximum edit distance for clustering (default: 1).

    Returns:
        Dictionary mapping cluster center barcodes to lists of member barcodes.

    Example:
        >>> counts = {"ACGT": 10, "ACGA": 2, "GGGG": 5}
        >>> clusters = cluster_barcodes(counts, distance_threshold=1)
        >>> clusters
        {"ACGT": ["ACGT", "ACGA"], "GGGG": ["GGGG"]}
    """
```

### 6.2 Update README.md

Add sections for:

- Development setup with uv
- Running tests
- Contributing guidelines
- Link to this implementation plan

---

## Implementation Priority Matrix

| Phase | Priority | Effort | Impact | Dependencies |
|-------|----------|--------|--------|--------------|
| Phase 1 (Build System) | HIGH | Medium | High | None |
| Phase 2 (Test Suite) | HIGH | High | High | Phase 1 |
| Phase 3 (Security Fixes) | CRITICAL | Low | Critical | None |
| Phase 4 (Code Quality) | MEDIUM | High | Medium | Phase 1, 2 |
| Phase 5 (Pre-commit Hooks) | MEDIUM | Low | Medium | Phase 1 |
| Phase 6 (Documentation) | LOW | Medium | Medium | Phase 4 |

## Recommended Implementation Order

1. **Immediate (Week 1):**
   - Phase 3.1: Fix command injection vulnerability
   - Phase 3.2-3.4: Fix other critical bugs
   - Phase 1.1-1.2: Set up pyproject.toml with hatchling and ruff

2. **Short-term (Weeks 2-3):**
   - Phase 1.3-1.4: Complete build system migration
   - Phase 2.1-2.3: Set up test infrastructure
   - Phase 5: Set up pre-commit hooks

3. **Medium-term (Weeks 4-6):**
   - Phase 2.4-2.5: Write priority unit tests
   - Phase 4.1: Add type hints to core modules
   - Phase 4.2: Extract constants

4. **Long-term (Ongoing):**
   - Phase 4.3-4.6: Refactoring and modernization
   - Phase 6: Documentation improvements
   - Increase test coverage to 80%+

---

## Files to Create

| File | Purpose |
|------|---------|
| `pyproject.toml` | Modern build configuration |
| `.pre-commit-config.yaml` | Pre-commit hooks |
| `tests/conftest.py` | Test fixtures |
| `tests/unit/test_umi_cluster.py` | Unit tests |
| `umierrorcorrect/constants.py` | Extracted constants |

## Files to Delete (after migration)

| File | Reason |
|------|--------|
| `setup.py` | Replaced by pyproject.toml |
| `setup.cfg` | Replaced by pyproject.toml |

## Files to Modify

| File | Changes |
|------|---------|
| All scripts in `umierrorcorrect/` | Add `main_cli()` entry points |
| `umierrorcorrect/umi_error_correct.py` | Security fixes, refactoring |
| `umierrorcorrect/src/check_args.py` | Bug fixes, resource leak |
| `README.md` | Development instructions |
| `CLAUDE.md` | Update build/test commands |

---

## Appendix: Current State Assessment

### Strengths

- Well-structured pipeline with clear separation of concerns
- Modular design allows running steps independently
- Published and validated algorithm (Clinical Chemistry 2022)
- Docker support for reproducibility
- Reasonable documentation for end users

### Weaknesses

- No automated testing (only manual integration test)
- No type hints (harder to maintain, no IDE support)
- Legacy build system (setup.py)
- Security vulnerability (shell=True)
- Code debt (long functions, hardcoded values, old patterns)
- No CI/CD pipeline
- Inconsistent code style

### Risk Assessment

- **High Risk:** Command injection vulnerability could allow arbitrary code execution
- **Medium Risk:** Lack of tests means regressions go undetected
- **Low Risk:** Code style issues affect maintainability but not functionality
