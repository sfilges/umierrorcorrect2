#!/usr/bin/env python3
"""Unified CLI for UMI Error Correct using Typer."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console

from umierrorcorrect.src.logging_config import get_logger, setup_logging
from umierrorcorrect.version import __version__

# Create main app and subcommand apps
app = typer.Typer(
    name="umierrorcorrect",
    help="Pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI).",
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
)

console = Console()
logger = get_logger("cli")


def version_callback(value: bool) -> None:
    """Show version and exit."""
    if value:
        console.print(f"[bold green]UMI Error Correct[/bold green] version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option("--version", "-v", callback=version_callback, is_eager=True, help="Show version and exit."),
    ] = False,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-V", help="Enable verbose output (DEBUG level)."),
    ] = False,
) -> None:
    """UMI Error Correct - Pipeline for analyzing barcoded amplicon sequencing data."""
    log_level = "DEBUG" if verbose else "INFO"
    setup_logging(level=log_level)


@app.command()
def run(
    read1: Annotated[Path, typer.Option("-r1", "--read1", help="Path to first FASTQ file (R1).")],
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")],
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    read2: Annotated[Optional[Path], typer.Option("-r2", "--read2", help="Path to second FASTQ file (R2).")] = None,
    bed_file: Annotated[
        Optional[Path], typer.Option("-bed", "--bed-file", help="Path to BED file defining targeted regions.")
    ] = None,
    umi_length: Annotated[int, typer.Option("-ul", "--umi-length", help="Length of UMI sequence.")] = 12,
    spacer_length: Annotated[
        int, typer.Option("-sl", "--spacer-length", help="Length of spacer between UMI and read.")
    ] = 0,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name for output files.")] = None,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Number of threads.")] = 1,
    edit_distance: Annotated[int, typer.Option("-d", "--edit-distance", help="Edit distance threshold for UMI clustering.")] = 1,
    dual_index: Annotated[bool, typer.Option("--dual-index", help="Use dual indices (UMIs on R1 and R2).")] = False,
    reverse_index: Annotated[bool, typer.Option("--reverse-index", help="UMI is on R2 instead of R1.")] = False,
    adapter_trimming: Annotated[bool, typer.Option("--trim", help="Perform 3' adapter trimming.")] = False,
    remove_large_files: Annotated[bool, typer.Option("--remove", help="Remove original FASTQ and BAM files.")] = False,
) -> None:
    """Run the complete UMI Error Correct pipeline.

    This command runs all pipeline steps: preprocessing, mapping, UMI error correction,
    consensus statistics, and variant calling.
    """
    from argparse import Namespace

    from umierrorcorrect.run_umierrorcorrect import main as run_pipeline

    # Convert to legacy args format
    args = Namespace(
        read1=str(read1),
        read2=str(read2) if read2 else None,
        reference_file=str(reference),
        output_path=str(output),
        bed_file=str(bed_file) if bed_file else None,
        umi_length=str(umi_length),
        spacer_length=str(spacer_length),
        sample_name=sample_name,
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
        consensus_frequency_threshold=None,
        group_method="position",
        consensus_method="position",
        count_cutoff=3,
        qvalue_threshold=0.01,
        fsize=3,
        vc_method="count",
        params_file=None,
        output_json=False,
    )

    logger.info("Starting UMI Error Correct pipeline")
    run_pipeline(args)
    logger.info("Pipeline complete!")


@app.command()
def preprocess(
    read1: Annotated[Path, typer.Option("-r1", "--read1", help="Path to first FASTQ file (R1).")],
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    umi_length: Annotated[int, typer.Option("-ul", "--umi-length", help="Length of UMI sequence.")],
    read2: Annotated[Optional[Path], typer.Option("-r2", "--read2", help="Path to second FASTQ file (R2).")] = None,
    spacer_length: Annotated[int, typer.Option("-sl", "--spacer-length", help="Length of spacer.")] = 0,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Number of threads.")] = 1,
    dual_index: Annotated[bool, typer.Option("--dual-index", help="Use dual indices.")] = False,
    reverse_index: Annotated[bool, typer.Option("--reverse-index", help="UMI is on R2.")] = False,
    adapter_trimming: Annotated[bool, typer.Option("--trim", help="Perform 3' adapter trimming.")] = False,
    adapter_sequence: Annotated[str, typer.Option("-a", "--adapter", help="Adapter sequence to trim.")] = "illumina",
    force: Annotated[bool, typer.Option("-f", "--force", help="Overwrite existing files.")] = False,
) -> None:
    """Preprocess FASTQ files by extracting UMIs and adding them to read headers."""
    from argparse import Namespace

    from umierrorcorrect.preprocess import main as run_preprocess

    args = Namespace(
        read1=str(read1),
        read2=str(read2) if read2 else None,
        output_path=str(output),
        umi_length=str(umi_length),
        spacer_length=str(spacer_length),
        sample_name=sample_name,
        num_threads=str(threads),
        dual_index=dual_index,
        reverse_index=reverse_index,
        adapter_trimming=adapter_trimming,
        adapter_sequence=adapter_sequence,
        force=force,
        tmpdir=None,
    )

    logger.info("Starting preprocessing")
    run_preprocess(args)
    logger.info("Preprocessing complete!")


@app.command()
def consensus(
    bam: Annotated[Path, typer.Option("-b", "--bam", help="Path to input BAM file.")],
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    bed_file: Annotated[Optional[Path], typer.Option("-bed", "--bed-file", help="Path to BED file.")] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    edit_distance: Annotated[int, typer.Option("-d", "--edit-distance", help="Edit distance threshold.")] = 1,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Number of threads.")] = 1,
    include_singletons: Annotated[bool, typer.Option("--singletons", help="Include singleton reads.")] = True,
    consensus_method: Annotated[str, typer.Option("-c", "--consensus-method", help="Consensus method.")] = "position",
) -> None:
    """Generate consensus sequences from UMI-tagged BAM files."""
    from argparse import Namespace

    from umierrorcorrect.umi_error_correct import main as run_consensus

    args = Namespace(
        bam_file=str(bam),
        output_path=str(output),
        bed_file=str(bed_file) if bed_file else None,
        sample_name=sample_name,
        edit_distance_threshold=edit_distance,
        num_threads=str(threads),
        include_singletons=include_singletons,
        consensus_method=consensus_method,
        group_method="position",
        position_threshold=20,
        indel_frequency_threshold=60.0,
        consensus_frequency_threshold=None,
        regions_from_bed=False,
        regions_from_tag=False,
        remove_large_files=False,
        output_json=False,
    )

    logger.info("Starting consensus generation")
    run_consensus(args)
    logger.info("Consensus generation complete!")


@app.command()
def stats(
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    cons_bam: Annotated[Optional[Path], typer.Option("-c", "--cons-bam", help="Path to consensus BAM file.")] = None,
    hist_file: Annotated[Optional[Path], typer.Option("--hist", help="Path to histogram file.")] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    output_raw: Annotated[bool, typer.Option("--raw", help="Output raw consensus group counts.")] = False,
) -> None:
    """Generate consensus statistics from processed BAM files."""
    from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics

    logger.info("Generating consensus statistics")
    run_get_consensus_statistics(
        str(output),
        str(cons_bam) if cons_bam else None,
        str(hist_file) if hist_file else None,
        output_raw,
        sample_name,
    )
    logger.info("Statistics generation complete!")


@app.command()
def variants(
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")],
    cons_file: Annotated[Optional[Path], typer.Option("--cons", help="Path to cons.tsv file.")] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    fsize: Annotated[int, typer.Option("-f", "--fsize", help="Family size cutoff.")] = 3,
    method: Annotated[str, typer.Option("-m", "--method", help="Variant calling method (count or bbmodel).")] = "count",
    count_cutoff: Annotated[int, typer.Option("--cutoff", help="Minimum count for variant calling.")] = 3,
) -> None:
    """Call variants from consensus sequences."""
    from argparse import Namespace

    from umierrorcorrect.call_variants import run_call_variants

    args = Namespace(
        output_path=str(output),
        reference_file=str(reference),
        cons_file=str(cons_file) if cons_file else None,
        sample_name=sample_name,
        fsize=fsize,
        vc_method=method,
        count_cutoff=count_cutoff,
        qvalue_threshold=0.01,
        params_file=None,
    )

    logger.info("Starting variant calling")
    run_call_variants(args)
    logger.info("Variant calling complete!")


def main_cli() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main_cli()
