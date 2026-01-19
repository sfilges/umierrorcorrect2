#!/usr/bin/env python3

import errno
import os
import re
import subprocess
from argparse import Namespace
from pathlib import Path
from typing import Literal


def check_output_directory(outdir: str) -> str:
    """Check if outdir exists, otherwise create it."""
    if os.path.isdir(outdir):
        return outdir
    else:
        os.mkdir(outdir)
        return outdir


def get_sample_name(filename: str, mode: Literal["single", "paired", "bam"]) -> str:
    """Get the sample name as the basename of the input files."""
    basename = Path(filename).name
    if mode == "single":
        sample_name = basename.removesuffix(".gz").removesuffix(".fastq")
    elif mode == "paired":
        sample_name = basename.removesuffix(".gz").removesuffix(".fastq")
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


def is_tool(name: str) -> bool:
    """Check if a command-line tool is available.

    Args:
        name: Name of the tool to check.

    Returns:
        True if the tool is available, False otherwise.
    """
    try:
        with open(os.devnull, "w") as devnull:
            subprocess.run(
                [name, "--version"],
                stdout=devnull,
                stderr=devnull,
                check=False,
            )
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def check_args_fastq(args: Namespace) -> Namespace:
    """Function for checking arguments."""
    args.output_path = check_output_directory(args.output_path)
    is_pigz = is_tool("pigz")
    is_gzip = is_tool("gzip")
    is_bwa = is_tool("bwa")
    if args.adapter_trimming:
        is_cutadapt = is_tool("cutadapt")
        if not is_cutadapt:
            raise ValueError('Cannot find program "cutadapt". Please install it and add it to the path.')
    if not is_bwa:
        raise ValueError('Cannot find program "bwa". Please install it and add it to the path.')
    if is_pigz:
        args.gziptool = "pigz"
    elif is_gzip:
        args.gziptool = "gzip"
    else:
        raise ValueError('Cannot find program "gzip" or "pigz". Install one of them and add to the path.')

    # determine the mode (single or paired)
    if not args.read2:
        args.mode = "single"
    else:
        args.mode = "paired"
    if not args.sample_name:
        args.sample_name = get_sample_name(args.read1, args.mode)
    if args.dual_index and not args.mode == "paired":
        raise ValueError("Dual index can only be used when both an R1 and R2 file are supplied, exiting.")
    if args.reverse_index and not args.mode == "paired":
        raise ValueError("Reverse index can only be used when both an R1 and R2 file are supplied, exiting")
    try:
        args.umi_length = int(args.umi_length)
    except ValueError as e:
        raise ValueError(f"Barcode length must be an integer: {e}") from e
    try:
        args.spacer_length = int(args.spacer_length)
    except ValueError as e:
        raise ValueError(f"Spacer length must be an integer: {e}") from e
    # check if fastq files exist
    if not os.path.isfile(args.read1):
        raise ValueError(f"The file specified as r1 ({args.read1}) does not exist.")
    if args.mode == "paired":
        if not os.path.isfile(args.read2):
            raise ValueError(f"The file specified as r2 ({args.read2}) does not exist.")
    # check if umis_in_header file exists
    output_path = Path(args.output_path)
    if args.mode == "paired":
        f1file = output_path / f"{args.sample_name}_R1_umis_in_header.fastq.gz"
        f2file = output_path / f"{args.sample_name}_R2_umis_in_header.fastq.gz"
        if f1file.is_file() or f2file.is_file():
            if not args.force:
                raise ValueError(
                    f"The file {f1file} already exists. Overwrite it by including --force in the command line"
                )
            else:
                f1file.unlink(missing_ok=True)
                f2file.unlink(missing_ok=True)
    elif args.mode == "single":
        f1file = output_path / f"{args.sample_name}_umis_in_header.fastq.gz"
        if f1file.is_file():
            if not args.force:
                raise ValueError(
                    f"The file {f1file} already exists. Overwrite it by including --force in the command line"
                )
            else:
                f1file.unlink()
    return args


def check_args_bam(args: Namespace) -> Namespace:
    """Function for checking arguments."""
    args.output_path = check_output_directory(args.output_path)
    basenamefile = args.read1
    if not args.sample_name:
        args.sample_name = get_sample_name(basenamefile, args.mode)
    if args.regions_from_bed and not args.bed_file:
        raise ValueError("To use option regions_from_bed a bedfile needs to be supplied, using -bed option")
    return args


if __name__ == "__main__":
    is_pigz = is_tool("pigz")
    is_gzip = is_tool("gzip")
    is_bwa = is_tool("bwa")
    print(is_pigz)
    print(is_gzip)
    print(is_bwa)
