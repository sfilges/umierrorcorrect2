#!/usr/bin/env python3
"""
Pipeline orchestration for umierrorcorrect.

This module provides the main entry point for running the complete UMI error
correction pipeline on a single sample.

Pipeline Flow:
    FASTQ -> preprocess -> align_bwa -> umi_error_correct -> get_consensus_statistics -> call_variants -> VCF

External Dependencies:
    - bwa: Required for sequence mapping (must have indexed reference genome)
    - gzip/pigz: For compression

:Authors: Tobias Osterlund, Stefan Filges
"""

from pathlib import Path

from umierrorcorrect2.align import align_bwa, check_bwa_index
from umierrorcorrect2.call_variants import run_call_variants
from umierrorcorrect2.core.check_args import check_args_fastq
from umierrorcorrect2.core.logging_config import get_logger
from umierrorcorrect2.core.utils import get_sample_name
from umierrorcorrect2.get_consensus_statistics import run_get_consensus_statistics
from umierrorcorrect2.preprocess import run_preprocessing
from umierrorcorrect2.umi_error_correct import run_umi_errorcorrect

logger = get_logger(__name__)


def run_pipeline(args):
    """Run the complete UMI error correction pipeline on a single sample.

    Orchestrates the full pipeline: preprocessing FASTQ files to extract UMIs,
    mapping reads to reference genome, performing UMI-based error correction
    to generate consensus sequences, computing statistics, and calling variants.

    Args:
        args: Namespace object containing pipeline parameters:
            - read1 (str): Path to R1 FASTQ file
            - read2 (str, optional): Path to R2 FASTQ file for paired-end
            - reference_file (str): Path to BWA-indexed reference genome
            - output_path (str): Output directory for all results
            - sample_name (str, optional): Sample name (derived from read1 if not provided)
            - umi_length (str): Length of UMI sequence
            - spacer_length (str): Nucleotides between UMI and read start
            - num_threads (str): Number of threads for mapping
            - bed_file (str, optional): BED file defining target regions
            - edit_distance_threshold (int): Edit distance for UMI clustering
            - dual_index (bool): UMIs present on both R1 and R2
            - reverse_index (bool): UMI on R2 instead of R1

    Returns:
        None. Results are written to the output directory:
            - {sample}_consensus_reads.bam: Consensus BAM file
            - {sample}.cons: Tab-separated consensus statistics
            - {sample}.vcf: Called variants
            - {sample}.hist: UMI family size histogram
    """
    # Set sample name if not provided
    if not args.sample_name:
        args.sample_name = get_sample_name(args.read1, args.mode)

    # Validate arguments
    args = check_args_fastq(args)
    check_bwa_index(args.reference_file)

    # -----------------------------------------
    # Run preprocessing
    # -----------------------------------------
    from umierrorcorrect2.models.models import PreprocessConfig

    # Get fastp config from args (may be None if not using batch mode or fastp disabled)
    fastp_config = getattr(args, "fastp_config", None)

    preprocess_config = PreprocessConfig(
        read1=Path(args.read1),
        read2=Path(args.read2) if args.read2 else None,
        output_path=Path(args.output_path),
        umi_length=int(args.umi_length),
        spacer_length=int(args.spacer_length),
        sample_name=args.sample_name,
        num_threads=int(args.num_threads),
        dual_index=args.dual_index,
        reverse_index=args.reverse_index,
        adapter_trimming=args.adapter_trimming,
        adapter_sequence=getattr(args, "adapter_sequence", "illumina"),
        force=getattr(args, "force", False),
        tmpdir=getattr(args, "tmpdir", None),
        fastp_config=fastp_config,
    )
    fastq_files, nseqs = run_preprocessing(preprocess_config)
    logger.info(f"Files: {' '.join(fastq_files)}, number of reads: {nseqs}")

    # -----------------------------------------
    # Run mapping
    # -----------------------------------------
    bam_file = align_bwa(
        args.num_threads, args.reference_file, fastq_files, args.output_path, args.sample_name, args.remove_large_files
    )
    args.bam_file = bam_file
    args.regions_from_tag = False

    # -----------------------------------------
    # Run UMI error correction
    # -----------------------------------------
    from umierrorcorrect2.models.models import UMIErrorCorrectConfig

    umi_config = UMIErrorCorrectConfig(
        reference_file=Path(args.reference_file),
        output_path=Path(args.output_path),
        bam_file=Path(bam_file),
        sample_name=args.sample_name,
        bed_file=Path(args.bed_file) if args.bed_file else None,
        position_threshold=getattr(args, "position_threshold", 20),
        edit_distance_threshold=args.edit_distance_threshold,
        num_threads=int(args.num_threads) if args.num_threads else None,
        include_singletons=getattr(args, "include_singletons", True),
        consensus_method=getattr(args, "consensus_method", "position"),
        regions_from_bed=getattr(args, "regions_from_bed", False),
        regions_from_tag=False,
        indel_frequency_threshold=getattr(args, "indel_frequency_threshold", 60.0),
        consensus_frequency_threshold=getattr(args, "consensus_frequency_threshold", 60.0),
        output_json=getattr(args, "output_json", False),
        remove_large_files=getattr(args, "remove_large_files", False),
    )
    run_umi_errorcorrect(umi_config)
    output_path = Path(args.output_path)
    cons_bam = str(output_path / f"{args.sample_name}_consensus_reads.bam")

    # -----------------------------------------
    # Run consensus statistics
    # -----------------------------------------
    run_get_consensus_statistics(args.output_path, cons_bam, None, False, args.sample_name)
    args.cons_file = None

    # -----------------------------------------
    # Run variant calling
    # -----------------------------------------
    run_call_variants(args)
    logger.info("Finished UMI Error Correct")
