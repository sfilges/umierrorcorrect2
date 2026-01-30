#!/usr/bin/env python3
"""Batch processing for UMI Error Correct.

Process multiple samples in parallel with optional fastp pre-filtering
and QC report generation.
"""

import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn

from umierrorcorrect2.core.constants import HISTOGRAM_SUFFIX
from umierrorcorrect2.core.logging_config import (
    add_file_handler,
    get_log_path,
    get_logger,
    setup_logging,
)
from umierrorcorrect2.models.models import FastpConfig, Sample
from umierrorcorrect2.qc import run_fastqc, run_multiqc

logger = get_logger("batch")
console = Console()


@dataclass
class ProcessingResult:
    """Result of processing a single sample."""

    sample_name: str
    success: bool
    output_dir: Path
    error_message: Optional[str] = None
    consensus_bam: Optional[Path] = None
    vcf_file: Optional[Path] = None
    hist_file: Optional[Path] = None


def discover_samples(input_dir: Path) -> list[Sample]:
    """Discover paired-end FASTQ samples in a directory.

    Looks for files matching common patterns:
    - *_R1*.fastq.gz / *_R2*.fastq.gz
    - *_1.fastq.gz / *_2.fastq.gz
    - *R1*.fastq.gz / *R2*.fastq.gz

    Args:
        input_dir: Directory to search for FASTQ files.

    Returns:
        List of Sample objects with discovered read pairs.
    """
    samples: list[Sample] = []
    seen_r1_files: set[Path] = set()

    # Patterns for R1 files: (glob_pattern, r1_marker_for_replace, r2_marker_for_replace, r1_marker_for_name)
    # r1_marker_for_name is the marker to find in the sample name (after stripping extensions)
    patterns = [
        ("*_R1*.fastq.gz", "_R1", "_R2", "_R1"),
        ("*_R1*.fastq", "_R1", "_R2", "_R1"),
        ("*_1.fastq.gz", "_1.", "_2.", "_1"),
        ("*_1.fastq", "_1.", "_2.", "_1"),
        ("*R1*.fastq.gz", "R1", "R2", "R1"),
        ("*R1*.fastq", "R1", "R2", "R1"),
    ]

    for pattern, r1_marker, r2_marker, name_marker in patterns:
        for r1_path in input_dir.rglob(pattern):
            if r1_path in seen_r1_files:
                continue
            seen_r1_files.add(r1_path)

            # Try to find matching R2
            r2_name = r1_path.name.replace(r1_marker, r2_marker)
            r2_path = r1_path.parent / r2_name

            # Extract sample name (everything before the R1 marker)
            sample_name = r1_path.stem
            for suffix in [".fastq", ".fq"]:
                if sample_name.endswith(suffix):
                    sample_name = sample_name[: -len(suffix)]
            # Remove R1 marker from sample name
            if name_marker in sample_name:
                sample_name = sample_name.split(name_marker)[0].rstrip("_")

            sample = Sample(
                name=sample_name,
                read1=r1_path,
                read2=r2_path if r2_path.exists() else None,
            )
            samples.append(sample)
            logger.debug(
                f"Discovered sample: {sample.name} (R1: {r1_path}, R2: {r2_path if r2_path.exists() else 'None'})"
            )

    # Sort by sample name for consistent ordering
    samples.sort(key=lambda s: s.name)

    return samples


def parse_sample_sheet(csv_path: Path) -> list[Sample]:
    """Parse a CSV/TSV sample sheet.

    Expected format:
        sample_name,read1,read2
        sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz

    Also supports TSV format (auto-detected by .tsv extension or tab detection).

    Args:
        csv_path: Path to the sample sheet file.

    Returns:
        List of Sample objects parsed from the file.

    Raises:
        ValueError: If required columns are missing or files don't exist.
    """
    samples: list[Sample] = []

    with csv_path.open() as f:
        # Auto-detect delimiter
        first_line = f.readline()
        f.seek(0)

        delimiter = "\t" if csv_path.suffix.lower() == ".tsv" or "\t" in first_line else ","

        reader = csv.DictReader(f, delimiter=delimiter)

        # Check required columns
        if reader.fieldnames is None:
            raise ValueError("Sample sheet is empty or has no header")

        fieldnames_lower = [f.lower() for f in reader.fieldnames]
        if "sample_name" not in fieldnames_lower and "sample" not in fieldnames_lower:
            raise ValueError("Sample sheet must have 'sample_name' or 'sample' column")
        if "read1" not in fieldnames_lower and "r1" not in fieldnames_lower:
            raise ValueError("Sample sheet must have 'read1' or 'r1' column")

        # Map column names (case-insensitive)
        col_map: dict[str, str] = {}
        for fn in reader.fieldnames:
            fn_lower = fn.lower()
            if fn_lower in ("sample_name", "sample"):
                col_map["sample_name"] = fn
            elif fn_lower in ("read1", "r1"):
                col_map["read1"] = fn
            elif fn_lower in ("read2", "r2"):
                col_map["read2"] = fn

        for row in reader:
            sample_name = row[col_map["sample_name"]].strip()
            read1 = Path(row[col_map["read1"]].strip())

            read2: Optional[Path] = None
            if "read2" in col_map and row.get(col_map["read2"]):
                read2_val = row[col_map["read2"]].strip()
                if read2_val:
                    read2 = Path(read2_val)

            # Sample model validates file existence automatically
            samples.append(Sample(name=sample_name, read1=read1, read2=read2))

    return samples


def process_sample(
    sample: Sample,
    reference: Path,
    output_dir: Path,
    bed_file: Optional[Path] = None,
    umi_length: int = 12,
    spacer_length: int = 0,
    threads: int = 4,
    edit_distance: int = 1,
    dual_index: bool = False,
    reverse_index: bool = False,
    adapter_trimming: bool = False,
    remove_large_files: bool = False,
    fastp_config: Optional[FastpConfig] = None,
) -> ProcessingResult:
    """Process a single sample through the UMI Error Correct pipeline.

    Args:
        sample: Sample to process.
        reference: Reference genome FASTA.
        output_dir: Output directory for this sample.
        bed_file: Optional BED file for targeted regions.
        umi_length: Length of UMI sequence.
        spacer_length: Length of spacer between UMI and read.
        threads: Number of threads.
        edit_distance: Edit distance threshold for UMI clustering.
        dual_index: Use dual indices (UMIs on R1 and R2).
        reverse_index: UMI is on R2 instead of R1.
        adapter_trimming: Perform 3' adapter trimming.
        remove_large_files: Remove original FASTQ and BAM files.
        fastp_config: Optional FastpConfig for quality filtering and adapter trimming.

    Returns:
        ProcessingResult with outcome of processing.
    """
    from argparse import Namespace

    from umierrorcorrect2.pipeline import run_pipeline

    sample_output_dir = output_dir / sample.name
    sample_output_dir.mkdir(parents=True, exist_ok=True)

    # Create args namespace matching the expected format
    # TODO: Change to pydantic validation. Many args not exposed to the user here.
    args = Namespace(
        read1=str(sample.read1),
        read2=str(sample.read2) if sample.read2 else None,
        reference_file=str(reference),
        output_path=str(sample_output_dir),
        bed_file=str(bed_file) if bed_file else None,
        umi_length=str(umi_length),
        spacer_length=str(spacer_length),
        sample_name=sample.name,
        num_threads=str(threads),
        edit_distance_threshold=edit_distance,
        dual_index=dual_index,
        reverse_index=reverse_index,
        adapter_trimming=adapter_trimming,
        adapter_sequence="illumina",
        remove_large_files=remove_large_files,
        force=False,
        regions_from_bed=False,
        include_singletons=True,
        position_threshold=20,
        indel_frequency_threshold=60.0,
        consensus_frequency_threshold=60.0,
        group_method="position",
        consensus_method="position",
        count_cutoff=3,
        qvalue_threshold=0.01,
        fsize=3,
        vc_method="count",
        params_file=None,
        output_json=False,
        tmpdir=None,
        fastp_config=fastp_config,
    )

    try:
        run_pipeline(args)

        # Check for expected output files
        consensus_bam = sample_output_dir / f"{sample.name}_consensus_reads.bam"
        vcf_file = sample_output_dir / f"{sample.name}.vcf"
        hist_file = sample_output_dir / f"{sample.name}{HISTOGRAM_SUFFIX}"

        return ProcessingResult(
            sample_name=sample.name,
            success=True,
            output_dir=sample_output_dir,
            consensus_bam=consensus_bam if consensus_bam.exists() else None,
            vcf_file=vcf_file if vcf_file.exists() else None,
            hist_file=hist_file if hist_file.exists() else None,
        )

    except Exception as e:
        logger.error(f"Failed to process sample {sample.name}: {e}")
        return ProcessingResult(
            sample_name=sample.name,
            success=False,
            output_dir=sample_output_dir,
            error_message=str(e),
        )


def _process_sample_wrapper(args: tuple) -> ProcessingResult:
    """Wrapper for process_sample to work with ProcessPoolExecutor."""
    (
        sample,
        reference,
        output_dir,
        bed_file,
        umi_length,
        spacer_length,
        threads_per_sample,
        edit_distance,
        dual_index,
        reverse_index,
        adapter_trimming,
        remove_large_files,
        fastp_config,
    ) = args

    # Configure logging in subprocess:
    # - Console at INFO to avoid DEBUG spam on terminal
    # - File at DEBUG to capture external tool stderr
    setup_logging(level="INFO")
    sample_output_dir = output_dir / sample.name
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(sample_output_dir)
    add_file_handler(log_path)

    return process_sample(
        sample=sample,
        reference=reference,
        output_dir=output_dir,
        bed_file=bed_file,
        umi_length=umi_length,
        spacer_length=spacer_length,
        threads=threads_per_sample,
        edit_distance=edit_distance,
        dual_index=dual_index,
        reverse_index=reverse_index,
        adapter_trimming=adapter_trimming,
        remove_large_files=remove_large_files,
        fastp_config=fastp_config,
    )


def batch_process(
    samples: list[Sample],
    reference: Path,
    output_dir: Path,
    bed_file: Optional[Path] = None,
    umi_length: int = 19,
    spacer_length: int = 16,
    threads: int = 8,
    samples_parallel: int = 2,
    edit_distance: int = 1,
    dual_index: bool = False,
    reverse_index: bool = False,
    fastp: bool = True,
    merge_reads: bool = True,
    phred_score: int = 20,
    run_qc: bool = False,
    adapter_trimming: bool = True,
    remove_large_files: bool = False,
) -> list[ProcessingResult]:
    """Process samples through the UMI Error Correct pipeline.

    Args:
        samples: List of samples to process.
        reference: Reference genome FASTA.
        output_dir: Base output directory.
        bed_file: Optional BED file for targeted regions.
        umi_length: Length of UMI sequence.
        spacer_length: Length of spacer between UMI and read.
        threads: Total number of threads.
        samples_parallel: Number of samples to process in parallel.
        edit_distance: Edit distance threshold for UMI clustering.
        dual_index: Use dual indices.
        reverse_index: UMI is on R2.
        fastp: Enable fastp quality filtering and adapter trimming.
        merge_reads: Merge overlapping reads with fastp.
        phred_score: Minimum Phred quality score for fastp.
        run_qc: Whether to run FastQC/MultiQC.
        adapter_trimming: Perform 3' adapter trimming.
        remove_large_files: Remove original FASTQ and BAM files.

    Returns:
        List of ProcessingResult objects.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples_dir = output_dir / "samples"
    samples_dir.mkdir(exist_ok=True)

    # Calculate threads per sample
    threads_per_sample = max(1, threads // samples_parallel)

    # Create FastpConfig if fastp is enabled (will be passed to run_preprocessing)
    fastp_config: Optional[FastpConfig] = None
    if fastp:
        # Determine UMI location from dual_index and reverse_index flags
        if dual_index:
            umi_loc = "per_read"
        elif reverse_index:
            umi_loc = "read2"
        else:
            umi_loc = "read1"

        fastp_config = FastpConfig(
            enabled=True,
            phred_score=phred_score,
            merge_reads=merge_reads,
            trim_adapters=adapter_trimming,
            threads=threads_per_sample,
            umi_enabled=True,
            umi_length=umi_length,
            umi_skip=spacer_length,
            umi_loc=umi_loc,  # type: ignore
        )

    # Run FastQC if requested
    if run_qc:
        qc_dir = output_dir / "qc_reports"
        fastqc_dir = qc_dir / "fastqc"

        all_fastq_files = []
        for sample in samples:
            all_fastq_files.append(sample.read1)
            if sample.read2:
                all_fastq_files.append(sample.read2)

        run_fastqc(all_fastq_files, fastqc_dir, threads=threads)

    # Process samples in parallel
    results: list[ProcessingResult] = []

    console.print(f"[bold blue]Processing {len(samples)} samples ({samples_parallel} in parallel)...[/bold blue]")

    # Prepare arguments for each sample
    process_args = [
        (
            sample,
            reference,
            samples_dir,
            bed_file,
            umi_length,
            spacer_length,
            threads_per_sample,
            edit_distance,
            dual_index,
            reverse_index,
            adapter_trimming,
            remove_large_files,
            fastp_config,
        )
        for sample in samples
    ]

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]Processing samples...", total=len(samples))

        with ProcessPoolExecutor(max_workers=samples_parallel) as executor:
            future_to_sample = {executor.submit(_process_sample_wrapper, args): args[0].name for args in process_args}

            for future in as_completed(future_to_sample):
                sample_name = future_to_sample[future]
                try:
                    result = future.result()
                    results.append(result)

                    if result.success:
                        progress.console.print(f"  [green]Completed:[/green] {sample_name}")
                    else:
                        progress.console.print(f"  [red]Failed:[/red] {sample_name} - {result.error_message}")

                except Exception as e:
                    logger.error(f"Error processing {sample_name}: {e}")
                    results.append(
                        ProcessingResult(
                            sample_name=sample_name,
                            success=False,
                            output_dir=samples_dir / sample_name,
                            error_message=str(e),
                        )
                    )
                    progress.console.print(f"  [red]Failed:[/red] {sample_name} - {e}")

                progress.update(task, advance=1)

    # Run MultiQC if QC was requested
    if run_qc:
        qc_dir = output_dir / "qc_reports"
        run_multiqc(output_dir, qc_dir)

    # Generate summary report
    write_batch_summary(results, output_dir)

    return results


def write_batch_summary(results: list[ProcessingResult], output_dir: Path) -> None:
    """Write a summary TSV file with processing results.

    Args:
        results: List of processing results.
        output_dir: Output directory for summary file.
    """
    summary_file = output_dir / "batch_summary.tsv"

    with summary_file.open("w") as f:
        f.write("sample_name\tstatus\toutput_dir\tconsensus_bam\tvcf_file\thist_file\terror_message\n")

        for result in sorted(results, key=lambda r: r.sample_name):
            status = "SUCCESS" if result.success else "FAILED"
            f.write(
                f"{result.sample_name}\t"
                f"{status}\t"
                f"{result.output_dir}\t"
                f"{result.consensus_bam or ''}\t"
                f"{result.vcf_file or ''}\t"
                f"{result.hist_file or ''}\t"
                f"{result.error_message or ''}\n"
            )

    logger.info(f"Batch summary written to {summary_file}")

    # Print summary to console
    successful = sum(1 for r in results if r.success)
    failed = len(results) - successful

    console.print()
    console.print("[bold]Batch Processing Summary:[/bold]")
    console.print(f"  [green]Successful:[/green] {successful}")
    console.print(f"  [red]Failed:[/red] {failed}")
    console.print(f"  [blue]Summary file:[/blue] {summary_file}")
