#!/usr/bin/env python3

import shutil
from argparse import Namespace
from pathlib import Path

from umierrorcorrect2.core.utils import check_output_directory, get_sample_name


def is_tool(name: str) -> bool:
    """Check if a command-line tool is available in PATH.

    Args:
        name: Name of the tool to check.

    Returns:
        True if the tool is available, False otherwise.
    """
    return shutil.which(name) is not None


def check_args_fastq(args: Namespace) -> Namespace:
    """Function for checking arguments."""
    args.output_path = check_output_directory(args.output_path)
    is_pigz = is_tool("pigz")
    is_gzip = is_tool("gzip")
    is_bwa = is_tool("bwa")
    # Only require cutadapt if adapter trimming is needed and fastp won't handle it
    fastp_config = getattr(args, "fastp_config", None)
    fastp_handles_adapters = (
        fastp_config is not None
        and getattr(fastp_config, "enabled", False)
        and getattr(fastp_config, "trim_adapters", False)
    )
    if args.adapter_trimming and not fastp_handles_adapters:
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
    if not Path(args.read1).is_file():
        raise ValueError(f"The file specified as r1 ({args.read1}) does not exist.")
    if args.mode == "paired" and not Path(args.read2).is_file():
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
