#!/usr/bin/env python3
"""Quality control functions for UMI Error Correct.

Provides wrappers for FastQC and MultiQC quality control tools.
"""

import subprocess
from pathlib import Path

from umierrorcorrect2.core.check_args import is_tool
from umierrorcorrect2.core.logging_config import get_logger, log_subprocess_stderr

logger = get_logger(__name__)


def run_fastqc(files: list[Path], output_dir: Path, threads: int = 4) -> bool:
    """Run FastQC on a list of FASTQ files.

    Args:
        files: List of FASTQ files to analyze.
        output_dir: Output directory for FastQC reports.
        threads: Number of threads.

    Returns:
        True if successful, False otherwise.
    """
    if not is_tool("fastqc"):
        logger.warning("fastqc not found in PATH, skipping QC")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["fastqc", "-o", str(output_dir), "-t", str(threads)]
    cmd.extend([str(f) for f in files if f.exists()])

    logger.info(f"Running FastQC on {len(files)} files")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_subprocess_stderr(result.stderr, "fastqc")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"FastQC failed: {e.stderr}")
        return False


def run_multiqc(input_dir: Path, output_dir: Path) -> bool:
    """Run MultiQC to aggregate QC reports.

    Args:
        input_dir: Directory containing QC reports to aggregate.
        output_dir: Output directory for MultiQC report.

    Returns:
        True if successful, False otherwise.
    """
    if not is_tool("multiqc"):
        logger.warning("multiqc not found in PATH, skipping aggregation")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["multiqc", str(input_dir), "-o", str(output_dir), "-f"]

    logger.info("Running MultiQC")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_subprocess_stderr(result.stderr, "multiqc")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"MultiQC failed: {e.stderr}")
        return False


# TODO: Generate QC summary: on-target fraction, coverage, etc. per sample
# TODO: What is in the "_target_coverage.txt" file?
# TODO: Move the downsampling logic to the QC module
# TODO: Add a function to calculate on-target fraction
# TODO: Estimate background error if mutation bed is provided and known mutations can be ruled out
