#!/usr/bin/env python3
"""Unified CLI for UMI Error Correct using Typer."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from umierrorcorrect2.core.constants import ASCII_ART, DEFAULT_FAMILY_SIZES_STR
from umierrorcorrect2.core.logging_config import add_file_handler, get_log_path, get_logger, setup_logging
from umierrorcorrect2.version import __version__

# Create main app and subcommand apps
app = typer.Typer(
    name="umierrorcorrect2",
    help="Pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI).",
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
    setup_logging(level=log_level)  # type: ignore


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
    adapter_trimming: Annotated[bool, typer.Option("--trim/--no-trim", help="Perform 3' adapter trimming.")] = True,
    adapter_sequence: Annotated[
        str,
        typer.Option("-a", "--adapter", help="Adapter sequence to trim (used by cutadapt; fastp uses auto-detection)."),
    ] = "illumina",
    force: Annotated[bool, typer.Option("-f", "--force", help="Overwrite existing files.")] = False,
    fastp: Annotated[
        bool, typer.Option("--fastp/--no-fastp", help="Enable fastp quality filtering and UMI extraction.")
    ] = True,
    fastp_phred: Annotated[int, typer.Option("--fastp-phred", "-q", help="Minimum Phred score for fastp.")] = 20,
    fastp_merge: Annotated[
        bool, typer.Option("--fastp-merge/--no-fastp-merge", help="Merge overlapping reads with fastp.")
    ] = True,
    fastp_extract_umi: Annotated[
        bool, typer.Option("--fastp-extract-umi/--no-fastp-extract-umi", help="Extract UMIs with fastp.")
    ] = True,
) -> None:
    """Preprocess FASTQ files by extracting UMIs and adding them to read headers."""
    from umierrorcorrect2.models.models import FastpConfig, PreprocessConfig
    from umierrorcorrect2.preprocess import run_preprocessing

    # Set up file logging
    output.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(output)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

    # Create FastpConfig if fastp is enabled
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
            phred_score=fastp_phred,
            merge_reads=fastp_merge,
            trim_adapters=adapter_trimming,
            threads=threads,
            umi_enabled=fastp_extract_umi,
            umi_length=umi_length,
            umi_skip=spacer_length,
            umi_loc=umi_loc,  # type: ignore
        )

    # Create unified PreprocessConfig - fastp handling is now internal
    config = PreprocessConfig(
        read1=read1,
        read2=read2,
        output_path=output,
        umi_length=umi_length,
        spacer_length=spacer_length,
        sample_name=sample_name,
        num_threads=threads,
        dual_index=dual_index,
        reverse_index=reverse_index,
        adapter_trimming=adapter_trimming,
        adapter_sequence=adapter_sequence,
        force=force,
        tmpdir=None,
        fastp_config=fastp_config,
    )

    logger.info("Starting preprocessing")
    _, nseqs = run_preprocessing(config)
    logger.info(f"Preprocessing complete! Processed {nseqs} sequences.")


@app.command()
def consensus(
    bam: Annotated[Path, typer.Option("-b", "--bam", help="Path to input BAM file.")],
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")],
    bed_file: Annotated[Optional[Path], typer.Option("-bed", "--bed-file", help="Path to BED file.")] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    edit_distance: Annotated[int, typer.Option("-d", "--edit-distance", help="Edit distance threshold.")] = 1,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Number of threads.")] = 1,
    include_singletons: Annotated[bool, typer.Option("--singletons", help="Include singleton reads.")] = True,
    consensus_method: Annotated[str, typer.Option("-c", "--consensus-method", help="Consensus method.")] = "position",
) -> None:
    """Generate consensus sequences from UMI-tagged BAM files."""
    from umierrorcorrect2.models.models import UMIErrorCorrectConfig
    from umierrorcorrect2.umi_error_correct import run_umi_errorcorrect

    # Set up file logging
    output.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(output)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

    config = UMIErrorCorrectConfig(
        reference_file=reference,
        output_path=output,
        bam_file=bam,
        sample_name=sample_name,
        bed_file=bed_file,
        edit_distance_threshold=edit_distance,
        num_threads=threads,
        include_singletons=include_singletons,
        consensus_method=consensus_method,  # type: ignore
    )

    logger.info("Starting consensus generation")
    run_umi_errorcorrect(config)
    logger.info("Consensus generation complete!")


@app.command()
def stats(
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    cons_bam: Annotated[Optional[Path], typer.Option("-c", "--cons-bam", help="Path to consensus BAM file.")] = None,
    bed_file: Annotated[
        Optional[Path], typer.Option("-bed", "--bed-file", help="Path to BED file for region annotations.")
    ] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    output_raw: Annotated[bool, typer.Option("--raw", help="Output raw consensus group counts.")] = False,
) -> None:
    """Generate consensus statistics from processed BAM files.

    Statistics are derived directly from the consensus BAM file - no separate
    stats file is required.
    """
    from umierrorcorrect2.get_consensus_statistics import run_get_consensus_statistics

    # Set up file logging
    output.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(output)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

    logger.info("Generating consensus statistics")
    run_get_consensus_statistics(
        str(output),
        str(cons_bam) if cons_bam else None,
        str(bed_file) if bed_file else None,
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

    from umierrorcorrect2.call_variants import run_call_variants

    # Set up file logging
    output.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(output)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

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

    from umierrorcorrect2.align import align_bwa
    from umierrorcorrect2.core.utils import check_output_directory, get_sample_name

    output_path = check_output_directory(str(output))

    # Set up file logging
    log_path = get_log_path(output_path)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

    if read2 is None:
        fastq_files = [str(read1)]
        mode = "single"
    else:
        fastq_files = [str(read1), str(read2)]
        mode = "paired"

    if not sample_name:
        sample_name = get_sample_name(str(read1), mode)  # type: ignore

    logger.info("Starting BWA mapping")
    align_bwa(threads, reference, [Path(f) for f in fastq_files], output_path, sample_name, remove_files)
    logger.info("Mapping complete!")


@app.command(name="filter-bam")
def filter_bam(
    infile: Annotated[Path, typer.Option("-i", "--infile", help="Path to input BAM file.")],
    outfile: Annotated[Path, typer.Option("-o", "--outfile", help="Path to output BAM file.")],
    cutoff: Annotated[int, typer.Option("-c", "--cutoff", help="Consensus depth cutoff.")] = 3,
) -> None:
    """Filter BAM file by removing reads below consensus depth threshold."""
    from umierrorcorrect2.core.filter import filter_bam as run_filter_bam

    logger.info("Filtering BAM file")
    run_filter_bam(str(infile), str(outfile), cutoff)
    logger.info("BAM filtering complete!")


@app.command(name="filter-cons")
def filter_cons(
    infile: Annotated[Path, typer.Option("-i", "--infile", help="Path to input cons.tsv file.")],
    depth: Annotated[int, typer.Option("-d", "--depth", help="Raw depth cutoff.")] = 150,
    family_sizes: Annotated[
        str, typer.Option("-f", "--family-sizes", help="Family sizes to include, comma-separated.")
    ] = DEFAULT_FAMILY_SIZES_STR,
    write_raw: Annotated[bool, typer.Option("--write-raw", help="Include raw reads in output.")] = False,
) -> None:
    """Filter consensus file by depth and family sizes."""
    from umierrorcorrect2.core.filter import filter_cons as run_filter_cons

    logger.info("Filtering consensus file")
    run_filter_cons(str(infile), depth, family_sizes, write_raw)
    logger.info("Consensus filtering complete!")


@app.command()
def downsampling(
    output: Annotated[Path, typer.Option("-o", "--output", help="Path to output directory.")],
    cons_bam: Annotated[Optional[Path], typer.Option("-c", "--cons-bam", help="Path to consensus BAM file.")] = None,
    bed_file: Annotated[
        Optional[Path], typer.Option("-bed", "--bed-file", help="Path to BED file for region annotations.")
    ] = None,
    sample_name: Annotated[Optional[str], typer.Option("-s", "--sample-name", help="Sample name.")] = None,
    fsizes: Annotated[
        str, typer.Option("-f", "--fsizes", help="Comma-separated family size thresholds to plot (e.g., '1,2,3,5').")
    ] = "1,2,3,5",
) -> None:
    """Generate downsampling analysis plots.

    Creates a saturation curve showing the number of UMI families observed at different
    sequencing depths. Multiple family size thresholds can be plotted to compare how
    different minimum family sizes affect saturation.

    Statistics are derived directly from the consensus BAM file - no separate
    stats file is required.

    A plateau in the curve indicates sequencing saturation - additional sequencing
    would not discover significantly more UMI families.
    """
    from umierrorcorrect2.downsampling import run_downsampling

    # Set up file logging
    output.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(output)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

    # Parse comma-separated family sizes
    plot_fsizes = [int(f.strip()) for f in fsizes.split(",")]

    logger.info("Generating downsampling plots")
    run_downsampling(
        str(output),
        str(cons_bam) if cons_bam else None,
        str(bed_file) if bed_file else None,
        plot_fsizes,
        sample_name,
    )
    logger.info("Downsampling analysis complete!")


@app.command()
def analysis(
    sample_sheet: Annotated[
        Path,
        typer.Option("--sample-sheet", "-s", help="Extended sample sheet CSV for analysis metadata."),
    ],
    results_dir: Annotated[
        Path,
        typer.Option("--results-dir", "-r", help="Directory containing umierrorcorrect2 output samples/."),
    ],
    output_dir: Annotated[
        Optional[Path],
        typer.Option("-o", "--output-dir", help="Output directory for analysis reports."),
    ] = None,
    family_size: Annotated[
        int,
        typer.Option("-f", "--family-size", help="Family size threshold for analysis."),
    ] = 3,
) -> None:
    """Run extended analysis on umierrorcorrect2 results.

    Performs mutation tracking, on-target calculation, and generates
    HTML reports for samples in the provided extended sample sheet.

    Examples:
        umierrorcorrect2 analysis -s samples_extended.csv -r results/samples/ -o reports/
    """
    from umierrorcorrect2.analysis import AnalysisSampleSheet, Analyzer, HTMLReporter

    # Set up file logging
    out_path = output_dir or results_dir / "analysis_reports"
    out_path.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(out_path)
    add_file_handler(log_path)
    logger.info(f"Logging analysis to {log_path}")

    try:
        with console.status("[cyan]Parsing extended sample sheet..."):
            sheet = AnalysisSampleSheet(csv_path=sample_sheet, base_path=sample_sheet.parent)
        console.print(f"[green]✓[/green] Loaded {len(sheet.samples)} samples for analysis")

        analyzer = Analyzer(family_size=family_size)
        reporter = HTMLReporter(output_dir=out_path)

        for sample in sheet.samples:
            # Each sample's core results should be in results_dir/sample_name
            sample_res_dir = results_dir / sample.name
            if not sample_res_dir.exists():
                console.print(
                    f"[yellow]Warning:[/yellow] Results directory not found for {sample.name}: {sample_res_dir}"
                )
                continue

            with console.status(f"[cyan]Analyzing {sample.name}..."):
                analyzer.analyze_sample(sample, sample_res_dir)

                report_file = out_path / f"{sample.name}_report.html"
                reporter.generate_report(sample, family_size, report_file)

            console.print(f"  [green]✓[/green] Analysis report: {report_file}")

        console.print(f"\n[bold green]Analysis complete![/bold green] Reports saved to {out_path}")

    except Exception as e:
        console.print(f"[red]Analysis failed:[/red] {e}")
        logger.exception("Analysis failed")
        raise typer.Exit(1) from e


@app.command(name="fit-model")
def fit_model(
    cons_file: Annotated[Path, typer.Option("--cons", help="Path to cons.tsv file.")],
    known_mutations_file: Annotated[
        Optional[Path],
        typer.Option("--known-mutations", help="File with known true mutation positions to exclude from fitting."),
    ] = None,
    out_file: Annotated[Path, typer.Option("-o", "--out", help="Output file for model parameters.")] = Path(
        "bgmodel.params"
    ),
    fsize: Annotated[int, typer.Option("-f", "--fsize", help="Family size cutoff.")] = 3,
) -> None:
    """Fit beta-binomial background model for variant calling."""
    from argparse import Namespace

    from umierrorcorrect2.core.fit_background_model import run_fit_bgmodel

    args = Namespace(
        cons_file=str(cons_file),
        known_mutations_file=str(known_mutations_file) if known_mutations_file else None,
        out_file=str(out_file),
        fsize=fsize,
    )

    logger.info("Fitting background model")
    run_fit_bgmodel(args)
    logger.info("Model fitting complete!")


@app.command()
def aggregate(
    results_dir: Annotated[
        Path,
        typer.Option("-r", "--results-dir", help="Directory containing sample subdirectories with pipeline outputs."),
    ],
    output_dir: Annotated[
        Optional[Path],
        typer.Option("-o", "--output-dir", help="Output directory (default: results_dir/aggregate/)."),
    ] = None,
    family_size: Annotated[
        int,
        typer.Option("-f", "--family-size", help="Consensus group size threshold to filter on."),
    ] = 3,
    regions_bed: Annotated[
        Optional[Path],
        typer.Option("-rb", "--regions-bed", help="Region BED file for re-annotating the 'Name' column."),
    ] = None,
    mutation_bed: Annotated[
        Optional[Path],
        typer.Option("-mb", "--mutation-bed", help="Mutation BED file to add is_mutation/alt_matches columns."),
    ] = None,
) -> None:
    """Aggregate consensus TSV files from multiple sample directories into a single table.

    Loads all *_cons.tsv files from sample subdirectories, filters to the specified
    consensus group size, and optionally annotates with region or mutation BED files.

    Examples:

        umierrorcorrect2 aggregate -r results/samples/ -o results/aggregate/ -f 3

        umierrorcorrect2 aggregate -r results/samples/ -mb mutations.bed
    """
    from umierrorcorrect2.analysis.post_processor import PostProcessor

    out_path = output_dir or results_dir / "aggregate"
    out_path.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(out_path)
    add_file_handler(log_path)

    try:
        pp = PostProcessor(results_dir=results_dir, family_size=family_size)
        with console.status("[cyan]Aggregating consensus data..."):
            df = pp.aggregate_cons(region_bed=regions_bed, mutation_bed=mutation_bed)

        if df.empty:
            console.print("[yellow]Warning:[/yellow] No consensus data found.")
            raise typer.Exit(1)

        out_file = out_path / "combined_cons.tsv"
        df.to_csv(out_file, sep="\t", index=False)
        n_samples = df["Sample Name"].nunique() if "Sample Name" in df.columns else "?"
        console.print(f"[green]✓[/green] Aggregated {n_samples} samples ({len(df)} rows) → {out_file}")

    except typer.Exit:
        raise
    except Exception as e:
        console.print(f"[red]Aggregate failed:[/red] {e}")
        logger.exception("Aggregate failed")
        raise typer.Exit(1) from e


@app.command()
def summarize(
    results_dir: Annotated[
        Path,
        typer.Option("-r", "--results-dir", help="Directory containing sample subdirectories with pipeline outputs."),
    ],
    output_dir: Annotated[
        Optional[Path],
        typer.Option("-o", "--output-dir", help="Output directory (default: results_dir/aggregate/)."),
    ] = None,
    family_size: Annotated[
        int,
        typer.Option("-f", "--family-size", help="Consensus group size threshold for mutation metrics."),
    ] = 3,
    regions_bed: Annotated[
        Optional[Path],
        typer.Option("-rb", "--regions-bed", help="Region BED file for on-target read numerator."),
    ] = None,
    mutation_bed: Annotated[
        Optional[Path],
        typer.Option("-mb", "--mutation-bed", help="Mutation BED file for mutation metrics and on-target numerator."),
    ] = None,
    sample_sheet: Annotated[
        Optional[Path],
        typer.Option("-s", "--samplesheet", help="Optional extended sample sheet CSV with metadata (ml_plasma etc.)."),
    ] = None,
) -> None:
    """Summarize UMI error correction results across samples.

    Computes on-target read fractions (using fastp.json for total reads) and,
    if a mutation BED is provided, per-mutation metrics (VAF, ctDNA ppm) across
    all sample subdirectories.

    Examples:

        umierrorcorrect2 summarize -r results/samples/ -o results/aggregate/

        umierrorcorrect2 summarize -r results/samples/ -mb mutations.bed -s samplesheet.csv
    """
    from umierrorcorrect2.analysis.models import AnalysisSampleSheet
    from umierrorcorrect2.analysis.post_processor import PostProcessor

    out_path = output_dir or results_dir / "aggregate"
    out_path.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(out_path)
    add_file_handler(log_path)

    try:
        pp = PostProcessor(results_dir=results_dir, family_size=family_size)

        # Load optional metadata from extended sample sheet
        ml_plasma_map: dict[str, float] | None = None
        if sample_sheet is not None:
            with console.status("[cyan]Parsing sample sheet..."):
                sheet = AnalysisSampleSheet(csv_path=sample_sheet, base_path=sample_sheet.parent)
            ml_plasma_map = {s.name: s.ml_plasma for s in sheet.samples if s.ml_plasma is not None} or None
            console.print(f"[green]✓[/green] Loaded {len(sheet.samples)} samples from sample sheet")

        # On-target fractions (always computed)
        with console.status("[cyan]Computing on-target fractions..."):
            on_target_df = pp.compute_on_target_fractions(region_bed=regions_bed, mutation_bed=mutation_bed)

        on_target_file = out_path / "on_target_summary.tsv"
        on_target_df.to_csv(on_target_file, sep="\t", index=False)
        console.print(f"[green]✓[/green] On-target summary → {on_target_file}")

        # Mutation metrics (only if mutation_bed provided)
        if mutation_bed is not None:
            with console.status("[cyan]Computing mutation metrics..."):
                cons_df = pp.aggregate_cons(mutation_bed=mutation_bed)

            if not cons_df.empty:
                metrics_df = pp.compute_mutation_metrics(cons_df, ml_plasma_map=ml_plasma_map)
                metrics_file = out_path / "mutation_metrics.tsv"
                metrics_df.to_csv(metrics_file, sep="\t", index=False)
                console.print(f"[green]✓[/green] Mutation metrics → {metrics_file}")
            else:
                console.print("[yellow]Warning:[/yellow] No consensus data found for mutation metrics.")

        console.print(f"\n[bold green]Summarize complete![/bold green] Results saved to {out_path}")

    except Exception as e:
        console.print(f"[red]Summarize failed:[/red] {e}")
        logger.exception("Summarize failed")
        raise typer.Exit(1) from e


@app.command()
def run(
    read1: Annotated[
        Optional[Path], typer.Option("-r1", "--read1", help="Path to first FASTQ file (R1) for single-sample mode.")
    ] = None,
    read2: Annotated[
        Optional[Path], typer.Option("-r2", "--read2", help="Path to second FASTQ file (R2) for single-sample mode.")
    ] = None,
    input_dir: Annotated[
        Optional[Path], typer.Option("-i", "--input-dir", help="Directory containing FASTQ files to process.")
    ] = None,
    sample_sheet: Annotated[
        Optional[Path],
        typer.Option("--sample-sheet", help="CSV/TSV sample sheet with sample_name,read1,read2 columns."),
    ] = None,
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")] = ...,  # type: ignore[assignment]
    output_dir: Annotated[Path, typer.Option("-o", "--output-dir", help="Output directory for all samples.")] = ...,  # type: ignore[assignment]
    bed_file: Annotated[
        Optional[Path], typer.Option("-rb", "--regions-bed", help="Path to BED file defining targeted regions.")
    ] = None,
    umi_length: Annotated[int, typer.Option("-ul", "--umi-length", help="Length of UMI sequence.")] = 19,
    spacer_length: Annotated[
        int, typer.Option("-sl", "--spacer-length", help="Length of spacer between UMI and read.")
    ] = 16,
    # TODO: Consider removing the user provided sample name option. Just use the basename of the input file.
    sample_name: Annotated[
        Optional[str], typer.Option("-s", "--sample-name", help="Sample name (for single-sample mode).")
    ] = None,
    threads: Annotated[int, typer.Option("-t", "--threads", help="Total number of threads to use.")] = 8,
    samples_parallel: Annotated[
        int, typer.Option("-j", "--jobs", help="Number of samples to process in parallel.")
    ] = 2,
    adapter_trimming: Annotated[
        bool,
        typer.Option(
            "--trim-adapters/--no-trim-adapters",
            help="Performs adapter trimming with fastp if used, otherwise with cutadapt.",
        ),
    ] = True,
    fastp: Annotated[bool, typer.Option("--fastp/--no-fastp", help="Enable fastp quality filtering.")] = True,
    fastp_merge_reads: Annotated[
        bool,
        typer.Option("--fastp-merge/--no-fastp-merge", help="Merge overlapping reads with fastp (paired-end only)."),
    ] = True,
    qc: Annotated[bool, typer.Option("--qc/--no-qc", help="Generate FastQC and MultiQC reports.")] = True,
    phred_score: Annotated[
        int, typer.Option("-q", "--phred-score", help="Minimum Phred quality score for fastp filtering.")
    ] = 20,
    edit_distance: Annotated[
        int, typer.Option("-d", "--edit-distance", help="Edit distance threshold for UMI clustering.")
    ] = 1,
    dual_index: Annotated[bool, typer.Option("--dual-index", help="Use dual indices (UMIs on R1 and R2).")] = False,
    reverse_index: Annotated[bool, typer.Option("--reverse-index", help="UMI is on R2 instead of R1.")] = False,
    remove_large_files: Annotated[bool, typer.Option("--remove", help="Remove original FASTQ and BAM files.")] = False,
) -> None:
    """Run the UMI Error Correct pipeline.

    Process samples through the complete pipeline: preprocessing, mapping, UMI error
    correction, consensus statistics, and variant calling.

    Input modes (exactly one required):
      --read1/-r1: Single-sample mode - process one FASTQ file directly
      --input-dir/-i: Directory mode - discover and process all FASTQ pairs
      --sample-sheet: Sample sheet mode - process samples listed in CSV/TSV

    Examples:

        # Single sample mode with single-end reads (r1) with fastp and qc enabled
        umierrorcorrect2 run -r1 sample_R1.fastq.gz -r genome.fa -o results/

        # Single sample mode with paired-end reads, without fastp (cutadapt handles adapter trimming)
        # and without qc (equivalent to the original umierrorcorrect)
        umierrorcorrect2 run -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -r genome.fa -o results/ --no-fastp --no-qc

        # Batch process all FASTQ files in a directory
        umierrorcorrect2 run -i /path/to/fastqs -r genome.fa -o results/

        # Batch process with sample sheet
        umierrorcorrect2 run --sample-sheet samples.csv -r genome.fa -o results/

        # With pre-filtering (fastp is enabled by default) but without qc and non-standard UMI configuration
        umierrorcorrect2 run -i /path/to/fastqs -r genome.fa -o results/ --no-qc -ul 12 -sl 8
    """
    from umierrorcorrect2.batch import batch_process, discover_samples, parse_sample_sheet
    from umierrorcorrect2.models.models import Sample

    # Validate input options - exactly one input mode must be provided
    input_modes = sum([read1 is not None, input_dir is not None, sample_sheet is not None])
    if input_modes == 0:
        console.print("[red]Error:[/red] One of --read1, --input-dir, or --sample-sheet must be provided.")
        raise typer.Exit(1)
    if input_modes > 1:
        console.print("[red]Error:[/red] Only one of --read1, --input-dir, or --sample-sheet can be provided.")
        raise typer.Exit(1)

    # Validate read2 is only used with read1
    if read2 is not None and read1 is None:
        console.print("[red]Error:[/red] --read2 can only be used with --read1.")
        raise typer.Exit(1)

    # Discover or parse samples based on input mode. Returns list of Sample objects.
    if read1 is not None:
        # Single-sample mode
        if not read1.exists():
            console.print(f"[red]Error:[/red] Read1 file not found: {read1}")
            raise typer.Exit(1)
        if read2 is not None and not read2.exists():
            console.print(f"[red]Error:[/red] Read2 file not found: {read2}")
            raise typer.Exit(1)

        # Derive sample name from filename if not provided
        if sample_name is None:
            # Extract sample name from R1 filename
            name = read1.stem
            for suffix in [".fastq", ".fq"]:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
            # Remove common R1 markers
            for marker in ["_R1", "R1", "_1"]:
                if marker in name:
                    name = name.split(marker)[0].rstrip("_")
                    break
            derived_sample_name = name
        else:
            derived_sample_name = sample_name

        samples = [Sample(name=derived_sample_name, read1=read1, read2=read2)]
    elif sample_sheet is not None:
        if not sample_sheet.exists():
            console.print(f"[red]Error:[/red] Sample sheet not found: {sample_sheet}")
            raise typer.Exit(1)
        try:
            with console.status("[cyan]Parsing sample sheet..."):
                samples = parse_sample_sheet(sample_sheet)
            console.print(f"[green]✓[/green] Loaded {len(samples)} sample(s)")
        except ValueError as e:
            console.print(f"[red]Error parsing sample sheet:[/red] {e}")
            raise typer.Exit(1) from None
    else:
        # input_dir mode
        if input_dir is None or not input_dir.exists():
            console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
            raise typer.Exit(1)
        with console.status("[cyan]Discovering samples..."):
            samples = discover_samples(input_dir)
        console.print(f"[green]✓[/green] Found {len(samples)} sample(s)")

    if not samples:
        console.print("[red]Error:[/red] No samples found.")
        raise typer.Exit(1)

    # Validate reference
    # TODO: Validate reference before sample search?
    if not reference.exists():
        console.print(f"[red]Error:[/red] Reference file not found: {reference}")
        raise typer.Exit(1)

    # Set up file logging
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = get_log_path(output_dir)
    add_file_handler(log_path)
    logger.info(f"Logging to {log_path}")

    # Determine processing mode for display
    if read1 is not None:
        mode_str = "Single Sample"
    elif input_dir is not None:
        mode_str = "Directory"
    else:
        mode_str = "Sample Sheet"

    # Print banner and configuration
    console.print(Text(ASCII_ART, style="bold green"))
    console.print(f"[dim]Version {__version__}[/dim]", justify="center")
    console.print()

    # Build configuration table
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_column("Key", style="dim")
    config_table.add_column("Value")

    config_table.add_row("Mode", mode_str)
    config_table.add_row("Samples", str(len(samples)))

    if mode_str == "Single Sample":
        config_table.add_row("Read1", str(read1))
        if read2:
            config_table.add_row("Read2", str(read2))
    elif mode_str == "Directory":
        config_table.add_row("Input", str(input_dir))
    else:
        config_table.add_row("Sample sheet", str(sample_sheet))

    config_table.add_row("Reference", str(reference))
    config_table.add_row("Output", str(output_dir))

    if len(samples) > 1:
        config_table.add_row("Threads", f"{threads} ({samples_parallel} parallel)")
    else:
        config_table.add_row("Threads", str(threads))

    # Build features list
    features = []
    if fastp:
        features.append("fastp")
    if fastp_merge_reads and fastp:
        features.append("merge reads")
    if adapter_trimming:
        features.append(f"trim adapters ({'fastp' if fastp else 'cutadapt'})")
    if qc:
        features.append("QC reports")
    if features:
        config_table.add_row("Features", ", ".join(features))

    config_table.add_row("UMI length", str(umi_length))
    config_table.add_row("Spacer length", str(spacer_length))

    console.print(Panel(config_table, title="[bold]Configuration[/bold]", border_style="blue"))
    console.print()

    logger.info(f"Starting processing of {len(samples)} sample(s)")

    # Run batch processing
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
        fastp=fastp,
        merge_reads=fastp_merge_reads,
        phred_score=phred_score,
        run_qc=qc,
        adapter_trimming=adapter_trimming,
        remove_large_files=remove_large_files,
    )

    # Exit with non-zero code if any sample failed
    failed_count = sum(1 for r in results if not r.success)
    if failed_count > 0:
        logger.warning(f"{failed_count} sample(s) failed to process")
        raise typer.Exit(1)

    logger.info("Processing complete!")


def main_cli() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main_cli()
