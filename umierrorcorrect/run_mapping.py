#!/usr/bin/env python3
import logging
import subprocess
import sys
from pathlib import Path

import pysam

from umierrorcorrect.core.logging_config import log_subprocess_stderr


def check_bwa_index(reference_file):
    """Check if BWA index files exists, otherwise create"""
    ref_path = Path(reference_file)
    if not ref_path.is_file():
        print(f"Reference genome file {reference_file} does not exist, exiting")
        sys.exit(1)
    else:
        if not Path(reference_file + ".bwt").is_file():  # check if index exists
            print(f"BWA index for reference genome file {reference_file} does not exist")
            answer = input("Do you want to create a BWA index now? (y/n) ".lower().strip())
            while not (answer == "y" or answer == "yes" or answer == "n" or answer == "no"):
                print("Answer yes or no ")
                answer = input("Do you want to create a BWA index now? (y/n) ".lower().strip())
            if answer[0] == "y":
                create_index = True
            else:
                create_index = False
                sys.exit(1)
            if create_index:
                a = subprocess.Popen(["bwa", "index", reference_file], stderr=subprocess.PIPE)
                _, stderr = a.communicate()
                log_subprocess_stderr(stderr, "bwa-index")


def run_mapping(num_threads, reference_file, fastq_files, output_path, sample_name, remove_large_files):
    """Run mapping with bwa to create a SAM file, then convert it to BAM, sort and index the file"""
    logging.info("Starting mapping with BWA")
    check_bwa_index(reference_file)
    output_base = Path(output_path) / sample_name
    sam_file = f"{output_base}.sam"
    bam_file = f"{output_base}.bam"
    sorted_bam = f"{output_base}.sorted.bam"
    logging.info(f"Creating output file: {sorted_bam}")
    if len(fastq_files) == 1:
        bwacommand = ["bwa", "mem", "-t", num_threads, reference_file, fastq_files[0]]
    if len(fastq_files) == 2:
        bwacommand = ["bwa", "mem", "-t", num_threads, reference_file, fastq_files[0], fastq_files[1]]

    with Path(sam_file).open("w") as g:
        p1 = subprocess.Popen(bwacommand, stdout=g, stderr=subprocess.PIPE)
        _, stderr = p1.communicate()
        log_subprocess_stderr(stderr, "bwa-mem")
    p1.wait()
    pysam.view("-Sb", "-@", num_threads, sam_file, "-o", bam_file, catch_stdout=False)

    pysam.sort("-@", num_threads, bam_file, "-o", sorted_bam, catch_stdout=False)
    pysam.index(sorted_bam, catch_stdout=False)
    Path(sam_file).unlink()
    Path(bam_file).unlink()
    if remove_large_files:
        Path(fastq_files[0]).unlink()
        if len(fastq_files) == 2:
            Path(fastq_files[1]).unlink()
    logging.info("Finished mapping")
    return sorted_bam
