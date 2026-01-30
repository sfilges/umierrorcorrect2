#!/usr/bin/env python3
"""
umierrorcorrect, preprocess.py - Preprocess FASTQ sequences.
============================================================

:Authors: Stefan Filges, Tobias Osterlund

Purpose
-------

Preprocess FASTQ files for UMI error correction. This module handles:
1. Quality filtering and adapter trimming (using fastp).
2. Adapter trimming with cutadapt (optional).
3. Extraction of Unique Molecular Identifiers (UMIs) from reads and moving them to the read header.
4. Handling of both single-end and paired-end sequencing data.

"""

import datetime
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from umierrorcorrect2.core.check_args import is_tool
from umierrorcorrect2.core.logging_config import get_logger, log_subprocess_stderr
from umierrorcorrect2.core.read_fastq_records import read_fastq, read_fastq_paired_end
from umierrorcorrect2.models.models import FastpConfig, FastpResult, PreprocessConfig

logger = get_logger(__name__)


def run_fastp(
    read1: Path,
    read2: Optional[Path],
    output_dir: Path,
    sample_name: str,
    config: FastpConfig,
) -> Optional[FastpResult]:
    """Run fastp on FASTQ files for quality filtering.

    Args:
        read1: Path to first FASTQ file (R1).
        read2: Path to second FASTQ file (R2), or None for single-end.
        output_dir: Directory for filtered FASTQ files.
        sample_name: Sample name for output file naming.
        config: FastpConfig with filtering options.

    Returns:
        FastpResult with paths to filtered files, or None if fastp not available.
    """
    if not is_tool("fastp"):
        logger.warning("fastp not found in PATH, skipping pre-filtering")
        return None

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine output file naming based on whether UMI extraction is enabled
    umi_extraction = config.umi_enabled and config.umi_length > 0

    # Initialize output paths
    filtered_r1: Optional[Path] = None
    filtered_r2: Optional[Path] = None
    merged_reads: Optional[Path] = None

    fastp_json = output_dir / f"{sample_name}.fastp.json"
    fastp_html = output_dir / f"{sample_name}.fastp.html"

    cmd = [
        "fastp",
        "--in1",
        str(read1),
        "--json",
        str(fastp_json),
        "--html",
        str(fastp_html),
        "--qualified_quality_phred",
        str(config.phred_score),
        "--thread",
        str(config.threads),
    ]

    # Add adapter trimming for paired-end if enabled, single-end is automatically detected
    if config.trim_adapters and read2:
        cmd.append("--detect_adapter_for_pe")

    # Disable adapter trimming if not enabled (it's auto enabled otherwise)
    if not config.trim_adapters:
        cmd.append("--disable_adapter_trimming")

    # Add UMI extraction if enabled
    if umi_extraction:
        # TODO: What happens if umi_loc is in ["read2", "per_read"], but read2 is not provided?
        cmd.extend(
            [
                "--umi",
                "--umi_loc",
                config.umi_loc,
                "--umi_len",
                str(config.umi_length),
            ]
        )
        if config.umi_skip > 0:
            cmd.extend(["--umi_skip", str(config.umi_skip)])

    # Handle Paired-end vs Single-end and Merging logic
    if read2:
        cmd.extend(["--in2", str(read2)])

        if config.merge_reads:
            # Merge mode
            if umi_extraction:
                merged_reads_name = f"{sample_name}_umis_in_header.fastq.gz"
            else:
                merged_reads_name = f"{sample_name}.merged.fastq.gz"

            merged_reads = output_dir / merged_reads_name
            cmd.extend(["--merge", "--merged_out", str(merged_reads)])

            # Optionally keep unmerged reads
            if config.keep_unmerged:
                if umi_extraction:
                    filtered_r1 = output_dir / f"{sample_name}_R1_unmerged_umis_in_header.fastq.gz"
                    filtered_r2 = output_dir / f"{sample_name}_R2_unmerged_umis_in_header.fastq.gz"
                else:
                    filtered_r1 = output_dir / f"{sample_name}.unmerged.R1.fastq.gz"
                    filtered_r2 = output_dir / f"{sample_name}.unmerged.R2.fastq.gz"

                cmd.extend(["--out1", str(filtered_r1), "--out2", str(filtered_r2)])
        else:
            # Paired-end NO merge
            if umi_extraction:
                filtered_r1 = output_dir / f"{sample_name}_R1_umis_in_header.fastq.gz"
                filtered_r2 = output_dir / f"{sample_name}_R2_umis_in_header.fastq.gz"
            else:
                filtered_r1 = output_dir / f"{sample_name}.filtered.R1.fastq.gz"
                filtered_r2 = output_dir / f"{sample_name}.filtered.R2.fastq.gz"

            cmd.extend(["--out1", str(filtered_r1), "--out2", str(filtered_r2)])

    else:
        # Single-end mode
        if umi_extraction:
            filtered_r1 = output_dir / f"{sample_name}_umis_in_header.fastq.gz"
        else:
            filtered_r1 = output_dir / f"{sample_name}.filtered.R1.fastq.gz"

        cmd.extend(["--out1", str(filtered_r1)])

    logger.info(f"Running fastp on {sample_name}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_subprocess_stderr(result.stderr, "fastp")

        return FastpResult(
            filtered_read1=filtered_r1,
            filtered_read2=filtered_r2,
            merged_reads=merged_reads,
            fastp_json=fastp_json,
            fastp_html=fastp_html,
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"fastp failed for {sample_name}: {e.stderr}")
        return None


def generate_random_dir(tmpdir: str) -> str:
    """Generate a directory for storing temporary files, using a timestamp."""
    Path(tmpdir).mkdir(parents=True, exist_ok=True)
    prefix = f"r{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}_"
    return tempfile.mkdtemp(prefix=prefix, dir=str(tmpdir))


def run_unpigz(filename: str, tmpdir: str, num_threads: int, program: str) -> str:
    """Unzip fastq.gz files using parallel gzip (pigz) or gunzip."""
    input_path = Path(filename)
    outfilename = Path(tmpdir) / input_path.name.removesuffix(".gz")

    if program == "pigz":
        command = ["unpigz", "-p", str(num_threads), "-c", filename]
    else:
        command = ["gunzip", "-c", filename]

    with outfilename.open("w") as g:
        result = subprocess.run(command, stdout=g, stderr=subprocess.PIPE, check=True)
        log_subprocess_stderr(result.stderr, program)

    return str(outfilename)


def run_pigz(filename: str, num_threads: int, program: str) -> None:
    """Zip fastq files using parallel gzip (pigz) or gzip in-place."""
    if program == "pigz":
        command = ["pigz", "-p", str(num_threads), filename]
    else:
        command = ["gzip", filename]

    result = subprocess.run(command, capture_output=True, check=True)
    log_subprocess_stderr(result.stderr, program)


def run_cutadapt(
    r1file: str,
    r2file: str | None,
    output_path: Path,
    sample_name: str,
    adapter_sequence: str,
    mode: str,
) -> tuple[str, str | None]:
    """Run cutadapt to trim adapters."""
    if adapter_sequence.lower() == "illumina":
        adapter = "AGATCGGAAGAGC"
    elif adapter_sequence.lower() == "nextera":
        adapter = "CTGTCTCTTATA"
    elif adapter_sequence.lower() == "small-rna":
        adapter = "ATGGAATTCTCG"
    else:
        adapter = adapter_sequence.upper()

    if mode == "single":
        outfilename = str(output_path / f"{sample_name}_trimmed.fastq")
        command = ["cutadapt", "-a", adapter, "-o", outfilename, "-O", "3", "-m", "20", r1file]
        outfile1 = outfilename
        outfile2 = None
    else:
        outfile1 = str(output_path / f"{sample_name}_R1_trimmed.fastq")
        outfile2 = str(output_path / f"{sample_name}_R2_trimmed.fastq")
        if r2file is None:
            raise ValueError("r2file cannot be None in paired mode")
        command = [
            "cutadapt",
            "-a",
            adapter,
            "-A",
            adapter,
            "-o",
            outfile1,
            "-p",
            outfile2,
            "-O",
            "3",
            "-m",
            "20",
            r1file,
            r2file,
        ]

    logger.info(f"Performing adapter trimming using cutadapt with adapter sequence {adapter}")
    result = subprocess.run(command, capture_output=True, check=True)
    log_subprocess_stderr(result.stderr, "cutadapt")

    # Clean up input files
    Path(r1file).unlink()
    if mode != "single" and r2file:
        Path(r2file).unlink()

    return outfile1, outfile2


def preprocess_se(infilename: str, outfilename: str, barcode_length: int, spacer_length: int) -> int:
    """Run the preprocessing for single end data (one fastq file)."""
    with Path(infilename).open() as f, Path(outfilename).open("w") as g:
        read_start = barcode_length + spacer_length
        nseqs = 0
        for name, seq, qual in read_fastq(f):
            nseqs += 1
            barcode = seq[:barcode_length]
            # g.write(name+':'+barcode+'\n'+rest+'\n'+qualname+'\n'+qual[12+11:]+'\n')
            parts = name.split()
            newname = ":".join([parts[0], barcode]) + " " + parts[-1]
            g.write("\n".join([newname, seq[read_start:], "+", qual[read_start:]]) + "\n")
    return nseqs


def preprocess_pe(
    r1file: str,
    r2file: str,
    outfile1: str,
    outfile2: str,
    barcode_length: int,
    spacer_length: int,
    dual_index: bool,
) -> int:
    """Run the preprocessing for paired end data (two fastq files)."""
    read_start = barcode_length + spacer_length
    with (
        Path(r1file).open() as f1,
        Path(r2file).open() as f2,
        Path(outfile1).open("w") as g1,
        Path(outfile2).open("w") as g2,
    ):
        nseqs = 0
        for name1, seq1, qual1, name2, seq2, qual2 in read_fastq_paired_end(f1, f2):
            nseqs += 1
            if dual_index:
                barcode = seq1[:barcode_length] + seq2[:barcode_length]
            else:
                barcode = seq1[:barcode_length]
            parts1 = name1.split()
            parts2 = name2.split()
            newname1 = ":".join([parts1[0], barcode]) + " " + parts1[-1]
            newname2 = ":".join([parts2[0], barcode]) + " " + parts2[-1]
            g1.write("\n".join([newname1, seq1[read_start:], "+", qual1[read_start:]]) + "\n")
            if dual_index:
                g2.write("\n".join([newname2, seq2[read_start:], "+", qual2[read_start:]]) + "\n")
            else:
                g2.write("\n".join([newname2, seq2, "+", qual2]) + "\n")
    return 2 * nseqs


def prepare_input_files(
    input_read1: Path,
    input_read2: Optional[Path],
    effective_mode: str,
    tmpdir: str,
    num_threads: int,
    gziptool: str,
) -> tuple[str, Optional[str], bool]:
    """Unzip input files if necessary.

    Args:
        input_read1: Path to first input file.
        input_read2: Path to second input file (optional).
        effective_mode: "single" or "paired".
        tmpdir: Directory for temporary files.
        num_threads: Number of threads for pigz.
        gziptool: "pigz" or "gzip".

    Returns:
        Tuple of (r1file path, r2file path, removerfiles flag).
    """
    if not str(input_read1).endswith("gz"):
        r1file = str(input_read1)
        removerfiles = False
        r2file = str(input_read2) if input_read2 else None

        if effective_mode == "paired" and not input_read2:
            raise ValueError("Read2 not provided in paired mode")
    else:
        removerfiles = True
        if effective_mode == "paired":
            if not input_read2:
                raise ValueError("Read2 not provided in paired mode")
            r1file = run_unpigz(str(input_read1), tmpdir, num_threads, gziptool)
            r2file = run_unpigz(str(input_read2), tmpdir, num_threads, gziptool)
        else:
            r1file = run_unpigz(str(input_read1), tmpdir, num_threads, gziptool)
            r2file = None

    return r1file, r2file, removerfiles


def process_umi_extraction(
    r1file: str,
    r2file: Optional[str],
    effective_mode: str,
    config: PreprocessConfig,
    removerfiles: bool,
    tmpdir: str,
) -> tuple[list[str], int]:
    """Extract UMIs and cleanup temporary files.

    Args:
        r1file: Path to R1 file (unzipped).
        r2file: Path to R2 file (unzipped, optional).
        effective_mode: "single" or "paired".
        config: Preprocess config.
        removerfiles: Whether to remove input files after processing.
        tmpdir: Temporary directory path.

    Returns:
        Tuple of (list of output FASTQ files, number of sequences).
    """
    output_path = Path(config.output_path)

    if effective_mode == "single":
        outfilename = str(output_path / f"{config.sample_name}_umis_in_header.fastq")
        nseqs = preprocess_se(r1file, outfilename, config.umi_length, config.spacer_length)
        run_pigz(outfilename, config.num_threads, config.gziptool)

        if removerfiles:
            Path(r1file).unlink()
        Path(tmpdir).rmdir()
        fastqfiles = [f"{outfilename}.gz"]

    else:
        if r2file is None:
            raise ValueError("r2file is None in paired mode")

        if config.reverse_index:
            # switch forward and reverse read
            r1filetmp = r1file
            r1file = r2file
            r2file = r1filetmp
            outfile1 = str(output_path / f"{config.sample_name}_R2_umis_in_header.fastq")
            outfile2 = str(output_path / f"{config.sample_name}_R1_umis_in_header.fastq")
        else:
            outfile1 = str(output_path / f"{config.sample_name}_R1_umis_in_header.fastq")
            outfile2 = str(output_path / f"{config.sample_name}_R2_umis_in_header.fastq")

        nseqs = preprocess_pe(
            r1file,
            r2file,
            outfile1,
            outfile2,
            config.umi_length,
            config.spacer_length,
            config.dual_index,
        )
        run_pigz(outfile1, config.num_threads, config.gziptool)
        run_pigz(outfile2, config.num_threads, config.gziptool)

        r1path = Path(r1file)
        r2path = Path(r2file)
        if removerfiles and r1path.is_file():
            r1path.unlink()
        if removerfiles and r2path.is_file():
            r2path.unlink()
        Path(tmpdir).rmdir()
        fastqfiles = [outfile1 + ".gz", outfile2 + ".gz"]

    return fastqfiles, nseqs


def run_preprocessing(config: PreprocessConfig) -> tuple[list[str], int]:
    """Start preprocessing.

    Handles the complete preprocessing workflow:
    1. Optionally runs fastp for quality filtering and adapter trimming
    2. Performs UMI extraction (moving UMI from read to header)
    3. Conditionally runs cutadapt (only if adapter_trimming=True AND fastp didn't handle adapters)

    Args:
        config: PreprocessConfig with input/output paths, fastp settings, and options.

    Returns:
        Tuple of (list of output FASTQ file paths, number of sequences processed).
    """
    logger.info(f"Start preprocessing of sample {config.sample_name}")

    # Track inputs - may be modified by fastp
    input_read1 = config.read1
    input_read2 = config.read2
    fastp_trimmed_adapters = False

    # Step 1: Optional fastp preprocessing
    if config.fastp_config is not None and config.fastp_config.enabled:
        # Output to main directory if UMI extraction is enabled (final output), otherwise subdirectory
        umi_extraction = config.fastp_config.umi_enabled and config.fastp_config.umi_length > 0
        if umi_extraction:
            fastp_output_dir = config.output_path
        else:
            fastp_output_dir = config.output_path / "fastp_filtered"

        if config.sample_name is None:
            raise ValueError("Sample name must be provided for fastp processing")

        fastp_result = run_fastp(
            read1=config.read1,
            read2=config.read2,
            output_dir=fastp_output_dir,
            sample_name=config.sample_name,
            config=config.fastp_config,
        )

        if fastp_result is not None:
            # If fastp handled UMI extraction, return early with the output files
            if umi_extraction:
                logger.info("UMI extraction handled by fastp")

                if fastp_result.merged_reads and fastp_result.merged_reads.exists():
                    logger.info(f"Using merged reads with UMIs from fastp: {fastp_result.merged_reads}")
                    return [str(fastp_result.merged_reads)], -1
                elif fastp_result.filtered_read2:
                    logger.info("Using paired reads with UMIs from fastp")
                    return [str(fastp_result.filtered_read1), str(fastp_result.filtered_read2)], -1
                else:
                    logger.info("Using single-end reads with UMIs from fastp")
                    return [str(fastp_result.filtered_read1)], -1

            # Handle merged reads case (paired-end with merge enabled)
            if fastp_result.merged_reads and fastp_result.merged_reads.exists():
                input_read1 = fastp_result.merged_reads
                input_read2 = None
                logger.info(f"Using merged reads from fastp: {input_read1}")
            else:
                if fastp_result.filtered_read1 is None:
                    raise ValueError("fastp returned None for filtered_read1")
                input_read1 = fastp_result.filtered_read1
                input_read2 = fastp_result.filtered_read2
                logger.info("Using filtered reads from fastp")

            # Track if fastp handled adapter trimming
            if config.fastp_config.trim_adapters:
                fastp_trimmed_adapters = True
                logger.info("Adapter trimming handled by fastp")
        else:
            logger.warning("fastp failed, continuing with original reads")

    # Determine effective mode (may change from paired to single if fastp merged reads)
    effective_mode = "paired" if input_read2 else "single"

    # Step 2: Setup temp directory
    if config.tmpdir:
        newtmpdir = generate_random_dir(str(config.tmpdir))
    else:
        newtmpdir = generate_random_dir(str(config.output_path))

    # Step 3: Unzip input files (only reached if fastp didn't handle UMI extraction)
    r1file, r2file, removerfiles = prepare_input_files(
        input_read1, input_read2, effective_mode, newtmpdir, config.num_threads, config.gziptool
    )

    logger.info(f"Writing output files to {config.output_path}")

    # Step 4: Optional cutadapt (only if adapter trimming requested and fastp didn't handle it)
    if config.adapter_trimming is True and not fastp_trimmed_adapters:
        if config.sample_name is None:
            raise ValueError("Sample name must be provided for adapter trimming")
        output_path = Path(config.output_path)
        r2_input = r2file if effective_mode == "paired" else None
        r1file, r2file_out = run_cutadapt(
            r1file,
            r2_input,
            output_path,
            config.sample_name,
            config.adapter_sequence,
            effective_mode,
        )
        if effective_mode == "paired":
            if r2file_out is None:
                raise ValueError("Cutadapt returned None for R2 in paired mode")
            r2file = r2file_out

    # Step 5: UMI extraction
    return process_umi_extraction(r1file, r2file, effective_mode, config, removerfiles, newtmpdir)
