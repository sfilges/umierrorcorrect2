#!/usr/bin/env python3
"""
UMI error correct, run_umierrorcorrect.py - Run the pipeline

:Author: Tobias Osterlund

Purpose
-------

Run the pipeline

"""

import logging
from pathlib import Path

from umierrorcorrect.call_variants import run_call_variants
from umierrorcorrect.core.check_args import check_args_fastq, get_sample_name
from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics
from umierrorcorrect.preprocess import run_preprocessing
from umierrorcorrect.run_mapping import check_bwa_index, run_mapping
from umierrorcorrect.umi_error_correct import run_umi_errorcorrect


def main(args):
    if not args.sample_name:
        args.sample_name = get_sample_name(args.read1, args.mode)
    args = check_args_fastq(args)
    check_bwa_index(args.reference_file)
    # args=check_args_bam(args)
    fastq_files, nseqs = run_preprocessing(args)  # run preprocessing
    logging.info(f"Files: {' '.join(fastq_files)}, number of reads: {nseqs}")
    bam_file = run_mapping(
        args.num_threads, args.reference_file, fastq_files, args.output_path, args.sample_name, args.remove_large_files
    )  # run mapping
    args.bam_file = bam_file
    # print(args.bam_file)
    args.regions_from_tag = False
    run_umi_errorcorrect(args)  # run umi errorcorrect
    output_path = Path(args.output_path)
    cons_bam = str(output_path / f"{args.sample_name}_consensus_reads.bam")
    stat_filename = str(output_path / f"{args.sample_name}.hist")
    run_get_consensus_statistics(args.output_path, cons_bam, stat_filename, False, args.sample_name)
    args.cons_file = None
    # args.params_file=None
    run_call_variants(args)
    logging.info("Finished UMI Error Correct")
