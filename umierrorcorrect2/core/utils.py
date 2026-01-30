#!/usr/bin/env python3
"""Consolidated utility functions for UMIErrorCorrect.

This module provides common utility functions used across multiple modules
to avoid code duplication.
"""

import re
from pathlib import Path
from typing import Literal

import pysam


# TODO: This is currently not used, but should be used in the future
def get_percent_mapped_reads(num_fastq_reads, bamfile):
    """Get the number of mapped reads from the BAM index statistics.

    Args:
        num_fastq_reads (int): Total number of reads in the original FASTQ file.
        bamfile (str): Path to the BAM file.

    Returns:
        tuple: A tuple containing (number of mapped reads, ratio of mapped reads).
    """
    with pysam.AlignmentFile(bamfile, "rb") as f:
        stats = f.get_index_statistics()
        num_mapped = 0
        for s in stats:
            num_mapped += s.mapped
    ratio = (num_mapped / num_fastq_reads) * 1.0
    return (num_mapped, ratio)


def check_output_directory(outdir: str) -> str:
    """Check if outdir exists, otherwise create it.

    Args:
        outdir: Path to the output directory.

    Returns:
        The output directory path as a string.
    """
    outdir_path = Path(outdir)
    if outdir_path.is_dir():
        return outdir
    else:
        outdir_path.mkdir(parents=True, exist_ok=True)
        return outdir


def get_sample_name(filename: str, mode: Literal["single", "paired", "bam"] = "bam") -> str:
    """Get the sample name as the basename of the input files.

    Args:
        filename: Input filename to extract sample name from.
        mode: Processing mode - 'single' for single-end FASTQ, 'paired' for
              paired-end FASTQ, or 'bam' for BAM files.

    Returns:
        Extracted sample name.
    """
    basename = Path(filename).name

    if mode == "single":
        sample_name = basename.removesuffix(".gz").removesuffix(".fastq")
        if "_umis_in_header" in sample_name:
            sample_name = sample_name.replace("_umis_in_header", "")
    elif mode == "paired":
        sample_name = basename.removesuffix(".gz").removesuffix(".fastq")
        if "_umis_in_header" in sample_name:
            sample_name = sample_name.replace("_umis_in_header", "")
        if sample_name.endswith("_001"):
            sample_name = sample_name[:-4]
        if re.match(".*R[1-2]$", sample_name):
            sample_name = sample_name[:-2]
        sample_name = sample_name.rstrip("_")
        if re.search(".*_L00[0-9]$", sample_name):
            sample_name = sample_name[:-5]
    elif mode == "bam":
        sample_name = basename
        if ".sorted" in sample_name:
            sample_name = sample_name.replace(".sorted", "")
        sample_name = sample_name.replace(".bam", "")

    return sample_name


def get_sample_name_from_cons(cons_name: str) -> str:
    """Get the sample name from a consensus file path.

    Args:
        cons_name: Path to the consensus file.

    Returns:
        Extracted sample name.
    """
    sample_name = Path(cons_name).name
    sample_name = sample_name.replace("_cons.tsv", "")
    return sample_name


# TODO: This is only used for the variant calling pipeline. Maybe the variant calling pipeline should be refactored in the same module?
def parse_cons_file(filename: str, fsize: int = 3, include_position: bool = False) -> tuple:
    """
    Parse a consensus file and extract variant information.

    Args:
        filename: Path to the consensus file.
        fsize: Family size to filter on.
        include_position: If True, also return position information.

    Returns:
        If include_position is False: (frequencies, coverages, counts, data_lines)
        If include_position is True: (frequencies, coverages, counts, positions, data_lines)
    """
    n1 = []  # coverage
    f1 = []  # frequency
    c1 = []  # count
    posx = []  # positions
    data = []  # raw lines

    encoding = "iso-8859-1"
    with Path(filename).open(encoding=encoding) as f:
        header = f.readline()
        # Skip header if present
        if header.startswith("Sample Name"):
            pass  # Already read header
        else:
            # No header, process this line
            f.seek(0)
            f.readline()  # re-read to skip

        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            name = parts[3]

            if name not in "":
                famsize = parts[-4]
                if int(famsize) == fsize:
                    frac = float(parts[-2])
                    alt = parts[-1]
                    count = parts[-3]
                    if frac > 0 and alt not in "N":
                        cov = int(parts[-5])
                        f1.append(float(frac))
                        n1.append(int(cov))
                        c1.append(int(count))
                        if include_position:
                            pos = parts[1] + ":" + parts[2]
                            posx.append(pos)
                        data.append(line)

    if include_position:
        return (f1, n1, c1, posx, data)
    return (f1, n1, c1, data)
