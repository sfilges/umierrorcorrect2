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
    sample_name: Annotated[
        Optional[str], typer.Option("-s", "--sample-name", help="Sample name for output files.")
    ] = None,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Number of threads.")] = 1,
    edit_distance: Annotated[
        int, typer.Option("-d", "--edit-distance", help="Edit distance threshold for UMI clustering.")
    ] = 1,
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


@app.command()
def mapping(
    read1: Annotated[Path, typer.Option("-r1", "--read1", help="Path to first FASTQ file (R1).")],
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")],
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    read2: Annotated[Optional[Path], typer.Option("-r2", "--read2", help="Path to second FASTQ file (R2).")] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Number of threads.")] = 1,
    remove_files: Annotated[bool, typer.Option("--remove", help="Remove original FASTQ files after mapping.")] = False,
) -> None:
    """Run BWA mapping to reference genome."""

    from umierrorcorrect.run_mapping import check_output_directory, get_sample_name, run_mapping

    output_path = check_output_directory(str(output))

    if read2 is None:
        fastq_files = [str(read1)]
        mode = "single"
    else:
        fastq_files = [str(read1), str(read2)]
        mode = "paired"

    if not sample_name:
        sample_name = get_sample_name(str(read1), mode)

    logger.info("Starting BWA mapping")
    run_mapping(str(threads), str(reference), fastq_files, output_path, sample_name, remove_files)
    logger.info("Mapping complete!")


@app.command(name="filter-bam")
def filter_bam(
    infile: Annotated[Path, typer.Option("-i", "--infile", help="Path to input BAM file.")],
    outfile: Annotated[Path, typer.Option("-o", "--outfile", help="Path to output BAM file.")],
    cutoff: Annotated[int, typer.Option("-c", "--cutoff", help="Consensus depth cutoff.")] = 3,
) -> None:
    """Filter BAM file by removing reads below consensus depth threshold."""
    from umierrorcorrect.filter_bam import filter_bam as run_filter_bam

    logger.info("Filtering BAM file")
    run_filter_bam(str(infile), str(outfile), cutoff)
    logger.info("BAM filtering complete!")


@app.command(name="filter-cons")
def filter_cons(
    infile: Annotated[Path, typer.Option("-i", "--infile", help="Path to input cons.tsv file.")],
    depth: Annotated[int, typer.Option("-d", "--depth", help="Raw depth cutoff.")] = 150,
    family_sizes: Annotated[
        str, typer.Option("-f", "--family-sizes", help="Family sizes to include, comma-separated.")
    ] = "0,1,2,3,4,5,7,10,20,30",
    write_raw: Annotated[bool, typer.Option("--write-raw", help="Include raw reads in output.")] = False,
) -> None:
    """Filter consensus file by depth and family sizes."""
    from umierrorcorrect.filter_cons import filter_cons as run_filter_cons

    logger.info("Filtering consensus file")
    run_filter_cons(str(infile), depth, family_sizes, write_raw)
    logger.info("Consensus filtering complete!")


@app.command()
def downsampling(
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    cons_bam: Annotated[Optional[Path], typer.Option("-c", "--cons-bam", help="Path to consensus BAM file.")] = None,
    hist_file: Annotated[Optional[Path], typer.Option("--hist", help="Path to histogram file.")] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    fsize: Annotated[int, typer.Option("-f", "--fsize", help="Family size cutoff for downsampling.")] = 3,
) -> None:
    """Generate downsampling analysis plots."""
    from umierrorcorrect.downsampling_plots import run_downsampling

    logger.info("Generating downsampling plots")
    run_downsampling(
        str(output),
        str(cons_bam) if cons_bam else None,
        str(hist_file) if hist_file else None,
        fsize,
        sample_name,
    )
    logger.info("Downsampling analysis complete!")


@app.command(name="fit-model")
def fit_model(
    cons_file: Annotated[Optional[Path], typer.Option("--cons", help="Path to cons.tsv file.")] = None,
    nonbg_file: Annotated[
        Optional[Path], typer.Option("--nonbg", help="Path to file with non-background positions.")
    ] = None,
    out_file: Annotated[Path, typer.Option("-o", "--out", help="Output file for model parameters.")] = Path(
        "bgmodel.params"
    ),
    fsize: Annotated[int, typer.Option("-f", "--fsize", help="Family size cutoff.")] = 3,
) -> None:
    """Fit beta-binomial background model for variant calling."""
    from argparse import Namespace

    from umierrorcorrect.fit_background_model import run_fit_bgmodel

    args = Namespace(
        cons_file=str(cons_file) if cons_file else None,
        nonbgposfile=str(nonbg_file) if nonbg_file else None,
        out_file=str(out_file),
        fsize=fsize,
    )

    logger.info("Fitting background model")
    run_fit_bgmodel(args)
    logger.info("Model fitting complete!")


@app.command()
def batch(
    input_dir: Annotated[
        Optional[Path], typer.Option("-i", "--input-dir", help="Directory containing FASTQ files to process.")
    ] = None,
    sample_sheet: Annotated[
        Optional[Path],
        typer.Option("--sample-sheet", help="CSV/TSV sample sheet with sample_name,read1,read2 columns."),
    ] = None,
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")] = ...,
    output_dir: Annotated[Path, typer.Option("-o", "--output-dir", help="Output directory for all samples.")] = ...,
    bed_file: Annotated[
        Optional[Path], typer.Option("-bed", "--bed-file", help="Path to BED file defining targeted regions.")
    ] = None,
    umi_length: Annotated[int, typer.Option("-ul", "--umi-length", help="Length of UMI sequence.")] = 12,
    spacer_length: Annotated[
        int, typer.Option("-sl", "--spacer-length", help="Length of spacer between UMI and read.")
    ] = 0,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Total number of threads to use.")] = 8,
    samples_parallel: Annotated[
        int, typer.Option("-j", "--jobs", help="Number of samples to process in parallel.")
    ] = 2,
    prefilter: Annotated[Optional[str], typer.Option("--prefilter", help="Pre-filtering tool to use (fastp).")] = None,
    merge_reads: Annotated[
        bool, typer.Option("--merge-reads", help="Merge overlapping reads with fastp (paired-end only).")
    ] = False,
    phred_score: Annotated[
        int, typer.Option("-q", "--phred-score", help="Minimum Phred quality score for fastp filtering.")
    ] = 20,
    qc: Annotated[bool, typer.Option("--qc", help="Generate FastQC and MultiQC reports.")] = False,
    edit_distance: Annotated[
        int, typer.Option("-d", "--edit-distance", help="Edit distance threshold for UMI clustering.")
    ] = 1,
    dual_index: Annotated[bool, typer.Option("--dual-index", help="Use dual indices (UMIs on R1 and R2).")] = False,
    reverse_index: Annotated[bool, typer.Option("--reverse-index", help="UMI is on R2 instead of R1.")] = False,
) -> None:
    """Process multiple samples in batch.

    Discovers FASTQ pairs in a directory or reads from a sample sheet,
    then processes each sample through the UMI Error Correct pipeline in parallel.

    Examples:

        # Process all FASTQ files in a directory
        umierrorcorrect batch -i /path/to/fastqs -r genome.fa -o results/ -ul 12

        # With sample sheet
        umierrorcorrect batch --sample-sheet samples.csv -r genome.fa -o results/

        # With pre-filtering and QC
        umierrorcorrect batch -i /path/to/fastqs -r genome.fa -o results/ --prefilter fastp --qc -t 16
    """
    from umierrorcorrect.batch import batch_process, discover_samples, parse_sample_sheet

    # Validate input options
    if input_dir is None and sample_sheet is None:
        console.print("[red]Error:[/red] Either --input-dir or --sample-sheet must be provided.")
        raise typer.Exit(1)

    if input_dir is not None and sample_sheet is not None:
        console.print("[red]Error:[/red] Cannot use both --input-dir and --sample-sheet. Choose one.")
        raise typer.Exit(1)

    # Discover or parse samples
    if sample_sheet is not None:
        if not sample_sheet.exists():
            console.print(f"[red]Error:[/red] Sample sheet not found: {sample_sheet}")
            raise typer.Exit(1)
        try:
            samples = parse_sample_sheet(sample_sheet)
        except ValueError as e:
            console.print(f"[red]Error parsing sample sheet:[/red] {e}")
            raise typer.Exit(1) from None
    else:
        if not input_dir.exists():
            console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
            raise typer.Exit(1)
        samples = discover_samples(input_dir)

    if not samples:
        console.print("[red]Error:[/red] No samples found.")
        raise typer.Exit(1)

    # Validate reference
    if not reference.exists():
        console.print(f"[red]Error:[/red] Reference file not found: {reference}")
        raise typer.Exit(1)

    # Validate prefilter option
    if prefilter and prefilter.lower() not in ("fastp",):
        console.print(f"[red]Error:[/red] Unknown prefilter tool: {prefilter}. Supported: fastp")
        raise typer.Exit(1)

    console.print("[bold green]UMI Error Correct - Batch Processing[/bold green]")
    console.print(f"  Samples found: {len(samples)}")
    console.print(f"  Reference: {reference}")
    console.print(f"  Output: {output_dir}")
    console.print(f"  Threads: {threads} ({samples_parallel} samples in parallel)")
    if prefilter:
        console.print(f"  Pre-filter: {prefilter}")
    if qc:
        console.print("  QC reports: enabled")
    console.print()

    logger.info(f"Starting batch processing of {len(samples)} samples")

    results = batch_process(
        samples=samples,
        reference=reference,
        output_dir=output_dir,
        bed_file=bed_file,
        umi_length=umi_length,
        spacer_length=spacer_length,
        threads=threads,
        samples_parallel=samples_parallel,
        edit_distance=edit_distance,
        dual_index=dual_index,
        reverse_index=reverse_index,
        prefilter=prefilter,
        merge_reads=merge_reads,
        phred_score=phred_score,
        run_qc=qc,
    )

    # Exit with non-zero code if any sample failed
    failed_count = sum(1 for r in results if not r.success)
    if failed_count > 0:
        logger.warning(f"{failed_count} sample(s) failed to process")
        raise typer.Exit(1)

    logger.info("Batch processing complete!")


def main_cli() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main_cli()
