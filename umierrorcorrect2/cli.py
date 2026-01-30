#!/usr/bin/env python3
"""Unified CLI for UMI Error Correct using Typer."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
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
    reference: Annotated[Path, typer.Option("-r", "--reference", help="Path to reference genome FASTA.")] = ...,
    output_dir: Annotated[Path, typer.Option("-o", "--output-dir", help="Output directory for all samples.")] = ...,
    bed_file: Annotated[
        Optional[Path], typer.Option("-rb", "--regions-bed", help="Path to BED file defining targeted regions.")
    ] = None,
    # TODO: Mutation bed is currently not used, but will be required for reporting in the future.
    mut_bed: Annotated[
        Optional[Path],
        typer.Option("-mb", "--mutations-bed", help="Path to BED file defining patient-specific mutation positions."),
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
            samples = parse_sample_sheet(sample_sheet)
        except ValueError as e:
            console.print(f"[red]Error parsing sample sheet:[/red] {e}")
            raise typer.Exit(1) from None
    else:
        # input_dir mode
        if input_dir is None or not input_dir.exists():
            console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
            raise typer.Exit(1)
        samples = discover_samples(input_dir)

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

    # Print summary
    console.print(Text(ASCII_ART, style="bold green"))
    console.print(f"  Version: {__version__}")
    console.print("\n")
    console.print(f"  Mode: {mode_str}")
    console.print(f"  Samples: {len(samples)}")

    if mode_str == "Single Sample":
        console.print(f"  Read1: {read1}")
        console.print(f"  Read2: {read2}")
    elif mode_str == "Directory":
        console.print(f"  Input directory: {input_dir}")
    else:
        console.print(f"  Sample sheet: {sample_sheet}")

    console.print(f"  Reference: {reference}")
    console.print(f"  Output: {output_dir}")
    if len(samples) > 1:
        console.print(f"  Threads: {threads} ({samples_parallel} samples in parallel)")
    else:
        console.print(f"  Threads: {threads}")
    if fastp:
        console.print("  Preprocessing: fastp")
    if fastp_merge_reads and fastp:
        console.print("  Merge reads: enabled")
    if adapter_trimming and fastp:
        console.print("  Adapter trimming: fastp")
    if adapter_trimming and not fastp:
        console.print("  Adapter trimming: cutadapt")
    if qc:
        console.print("  QC reports: enabled")
    console.print(f"  UMI length: {umi_length}")
    console.print(f"  Spacer length: {spacer_length}")
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
