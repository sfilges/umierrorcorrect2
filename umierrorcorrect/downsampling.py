#!/usr/bin/env python3
"""Downsampling analysis functionality for UMI error correction.

This module provides functions for generating downsampling analysis plots
and statistics to assess sequencing depth and UMI family coverage.
"""

import random
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES, HISTOGRAM_SUFFIX
from umierrorcorrect.core.logging_config import get_logger
from umierrorcorrect.get_consensus_statistics import (
    get_stat,
    RegionConsensusStats,
    calculate_target_coverage,
    get_overall_statistics,
)

logger = get_logger(__name__)


def plot_downsampling(results_tot, fsize, plot_filename):
    """Generate a downsampling plot showing depth vs UMI family count.

    Args:
        results_tot: List of dictionaries containing downsampled statistics.
        fsize: Family size cutoff to use for plotting.
        plot_filename: Path to save the output plot.
    """
    x = []
    y = []
    for r in results_tot[0]:
        h = results_tot[0][r]
        x.append(h.total_reads[int(fsize)])
        y.append(h.umis[int(fsize)])
    plt.plot(x, y, "o-")
    plt.xlabel("Depth")
    plt.ylabel("Number of UMI families")
    plt.title("Downsampling plot")
    plt.box(False)
    plt.xlim(0, max(x) + 40000)
    plt.ylim(0, max(y) + 1000)
    plt.savefig(plot_filename)


def save_downsampled_table(all_results, tot_results, out_filename):
    """Save downsampled statistics to a tab-separated file.

    Args:
        all_results: List of dictionaries with per-region downsampled statistics.
        tot_results: List of dictionaries with total downsampled statistics.
        out_filename: Path to save the output file.
    """
    with Path(out_filename).open("w") as g:
        for r in tot_results[0]:
            text = tot_results[0][r].write_stats()
            lines = text.split("\n")
            for line in lines:
                g.write("downsampled" + str(r) + "\t" + line + "\n")
        for region in all_results:
            for r in region:
                text = region[r].write_stats()
                lines = text.split("\n")
                for line in lines:
                    g.write("downsampled" + str(r) + "\t" + line + "\n")


def downsample_reads_per_region(hist, _fraction, fsizes, onlyNamed=True):
    """Downsample reads and calculate statistics at various rates.

    Args:
        hist: List of region_cons_stat objects with histogram data.
        _fraction: Unused parameter (downsample rates are hardcoded).
        fsizes: List of family sizes to calculate statistics for.
        onlyNamed: If True, only process regions with names.

    Returns:
        List of dictionaries mapping downsample rates to region_cons_stat objects.
    """
    all_results = []
    downsample_rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for h in hist:
        run_analysis = True
        if onlyNamed and h.name == "":
            run_analysis = False
        if run_analysis:
            # hist: list of family sizes (e.g., [5, 3, 10, 2] means 4 families with 5, 3, 10, and 2 reads)
            # singletons: number of singletons
            num_families = len(h.hist)
            # tmpnames: list of family indices (e.g., [0, 1, 2, 3])
            tmpnames = np.array(range(0, num_families))

            # singnames: list of singleton indices (e.g., [3, 4])
            singnames = list(range(num_families, num_families + h.singletons))

            # Convert to np.array
            times = np.array(h.hist)

            # Expand the temp names by the number of reads per family
            # Result: [0,0,0,0,0, 1,1,1, 2,2,2,2,2,2,2,2,2,2]
            # (5 reads from family 0, 3 from family 1, 10 from family 2)
            reads = np.repeat(tmpnames, times, axis=0)  # expand to one entry per read
            reads = list(reads) + singnames
            results = {}
            for r in downsample_rates:
                # At 50% (r=0.5), we'd sample 10 random reads out of 20 total. Let's say we get:
                # ds_reads = [0, 0, 0, 1, 2, 2, 2, 2, 2, 3]
                ds_reads = random.sample(list(reads), round(r * len(reads)))  # noqa: S311 - downsample
                new_hist = Counter(ds_reads).values()  # collapse to one entry per UMI family
                # Counter: {0: 3, 1: 1, 2: 5, 3: 1}
                # new_hist = [3, 1, 5, 1] (sorted descending: [5, 3, 1, 1])
                new_hist = sorted(new_hist, reverse=True)  # sort
                new_singletons = list(new_hist).count(1)  # count singletons in new
                new_stat = region_cons_stat(h.regionid, h.pos, h.name, new_singletons, h.fsizes)
                new_stat.add_histogram(new_hist, fsizes)
                results[r] = new_stat
            all_results.append(results)
    return all_results


def run_downsampling(output_path, consensus_filename, stat_filename, fsize, samplename=None):
    """Run downsampling analysis and generate output files.

    Args:
        output_path: Directory for output files.
        consensus_filename: Path to consensus BAM file (auto-detected if None).
        stat_filename: Path to histogram file (auto-detected if None).
        fsize: Family size cutoff for downsampling plot.
        samplename: Sample name for output files (auto-detected if None).
    """
    logger.info("Getting consensus statistics")
    out_path = Path(output_path)

    if not consensus_filename:
        consensus_filename = str(list(out_path.glob("*_consensus_reads.bam"))[0])
    if not samplename:
        samplename = Path(consensus_filename).name.replace("_consensus_reads.bam", "")
    if not stat_filename:
        stat_filename = str(out_path / f"{samplename}{HISTOGRAM_SUFFIX}")

    hist = get_stat(consensus_filename, stat_filename)
    fsizes = list(DEFAULT_FAMILY_SIZES)[1:]  # Exclude 0, which is handled separately
    tot_results = region_cons_stat("All", "all_regions", "", 0, fsizes)
    for h in hist:
        # print(h.hist)
        tot_results.hist = tot_results.hist + h.hist
        tot_results.singletons += h.singletons

    downsample_rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    tot = downsample_reads_per_region([tot_results], downsample_rates, fsizes, False)
    all_results = downsample_reads_per_region(hist, downsample_rates, fsizes, True)
    filename = str(out_path / f"{samplename}_downsampled_coverage.txt")
    save_downsampled_table(all_results, tot, filename)
    filename = str(out_path / f"{samplename}_downsampled_plot.png")
    plot_downsampling(tot, fsize, filename)
