#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

import pysam

from umierrorcorrect.core.logging_config import get_logger, log_subprocess_stderr

logger = get_logger(__name__)


def _cleanup_files(*files: Path) -> None:
    """Remove files if they exist."""
    for f in files:
        if f.exists():
            f.unlink()


def check_bwa_index(reference_file: str | Path) -> None:
    """Check if BWA index files exists, otherwise create."""
    ref_path = Path(reference_file)
    if not ref_path.is_file():
        logger.error(f"Reference genome file {reference_file} does not exist, exiting")
        sys.exit(1)

    if not Path(str(reference_file) + ".bwt").is_file():  # check if index exists
        logger.warning(f"BWA index for reference genome file {reference_file} does not exist")
        answer = input("Do you want to create a BWA index now? (y/n) ").lower().strip()
        while answer not in ("y", "yes", "n", "no"):
            logger.warning("Answer yes or no")
            answer = input("Do you want to create a BWA index now? (y/n) ").lower().strip()
        if answer[0] != "y":
            sys.exit(1)
        try:
            logger.info("Creating BWA index...")
            result = subprocess.run(
                ["bwa", "index", reference_file],
                capture_output=True,
                check=True,
            )
            log_subprocess_stderr(result.stderr, "bwa-index")
        except subprocess.CalledProcessError as e:
            logger.error(f"bwa index failed: {e.stderr.decode() if e.stderr else 'Unknown error'}")
            sys.exit(1)


def align_bwa(
    num_threads: int,
    reference_file: str | Path,
    fastq_files: list[str | Path],
    output_path: str | Path,
    sample_name: str,
    remove_large_files: bool,
) -> str | None:
    """Align reads with BWA to create a SAM file, then convert it to BAM, sort and index the file."""
    logger.info("Starting alignment with BWA")

    # Validate inputs
    output_dir = Path(output_path)
    if not output_dir.is_dir():
        logger.error(f"Output directory {output_path} does not exist")
        return None

    if not 1 <= len(fastq_files) <= 2:
        logger.error(f"Expected 1 or 2 FASTQ files, got {len(fastq_files)}")
        return None

    check_bwa_index(reference_file)

    output_base = output_dir / sample_name
    sam_file = output_base.with_suffix(".sam")
    bam_file = output_base.with_suffix(".bam")
    sorted_bam = output_base.parent / f"{sample_name}.sorted.bam"
    logger.info(f"Creating output file: {sorted_bam}")

    bwacommand = ["bwa", "mem", "-t", str(num_threads), str(reference_file), *[str(f) for f in fastq_files]]

    try:
        with sam_file.open("w") as g:
            result = subprocess.run(bwacommand, stdout=g, stderr=subprocess.PIPE, check=True)
            log_subprocess_stderr(result.stderr, "bwa-mem")
    except subprocess.CalledProcessError as e:
        logger.error(f"bwa mem failed: {e.stderr.decode() if e.stderr else 'Unknown error'}")
        _cleanup_files(sam_file)
        return None

    try:
        pysam.view("-Sb", "-@", str(num_threads), str(sam_file), "-o", str(bam_file), catch_stdout=False)
        pysam.sort("-@", str(num_threads), str(bam_file), "-o", str(sorted_bam), catch_stdout=False)
        pysam.index(str(sorted_bam), catch_stdout=False)
    except pysam.SamtoolsError as e:
        logger.error(f"SAM/BAM processing failed: {e}")
        _cleanup_files(sam_file, bam_file)
        return None
    finally:
        _cleanup_files(sam_file, bam_file)

    if remove_large_files:
        for fastq in fastq_files:
            Path(fastq).unlink()

    logger.info("Finished alignment")
    return str(sorted_bam)
