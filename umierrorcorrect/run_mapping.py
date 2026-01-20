#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

import pysam

from umierrorcorrect.core.logging_config import get_logger, log_subprocess_stderr

logger = get_logger(__name__)


def check_bwa_index(reference_file):
    """Check if BWA index files exists, otherwise create."""
    ref_path = Path(reference_file)
    if not ref_path.is_file():
        logger.error(f"Reference genome file {reference_file} does not exist, exiting")
        sys.exit(1)

    if not Path(reference_file + ".bwt").is_file():  # check if index exists
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


def run_mapping(num_threads, reference_file, fastq_files, output_path, sample_name, remove_large_files):
    """Run mapping with bwa to create a SAM file, then convert it to BAM, sort and index the file"""
    logger.info("Starting mapping with BWA")
    check_bwa_index(reference_file)
    output_base = Path(output_path) / sample_name
    sam_file = f"{output_base}.sam"
    bam_file = f"{output_base}.bam"
    sorted_bam = f"{output_base}.sorted.bam"
    logger.info(f"Creating output file: {sorted_bam}")
    if len(fastq_files) == 1:
        bwacommand = ["bwa", "mem", "-t", num_threads, reference_file, fastq_files[0]]
    if len(fastq_files) == 2:
        bwacommand = ["bwa", "mem", "-t", num_threads, reference_file, fastq_files[0], fastq_files[1]]

    try:
        with Path(sam_file).open("w") as g:
            result = subprocess.run(bwacommand, stdout=g, stderr=subprocess.PIPE, check=True)
            log_subprocess_stderr(result.stderr, "bwa-mem")
    except subprocess.CalledProcessError as e:
        logger.error(f"bwa mem failed: {e.stderr.decode() if e.stderr else 'Unknown error'}")
        return None
    pysam.view("-Sb", "-@", num_threads, sam_file, "-o", bam_file, catch_stdout=False)

    pysam.sort("-@", num_threads, bam_file, "-o", sorted_bam, catch_stdout=False)
    pysam.index(sorted_bam, catch_stdout=False)
    Path(sam_file).unlink()
    Path(bam_file).unlink()
    if remove_large_files:
        Path(fastq_files[0]).unlink()
        if len(fastq_files) == 2:
            Path(fastq_files[1]).unlink()
    logger.info("Finished mapping")
    return sorted_bam
