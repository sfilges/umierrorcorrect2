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

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES
from umierrorcorrect.core.logging_config import get_logger
from umierrorcorrect.get_consensus_statistics import (
    RegionConsensusStats,
    get_stat,
)

logger = get_logger(__name__)


def plot_downsampling(results_tot, fsizes, plot_filename):
    """Generate a downsampling plot showing total sequencing depth vs UMI family count.

    Plots a separate line for each family size threshold, allowing comparison of
    how different thresholds affect saturation curves.

    Args:
        results_tot: List of dictionaries containing downsampled statistics.
        fsizes: List of family size thresholds to plot (e.g., [1, 2, 3, 5]).
        plot_filename: Path to save the output plot.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get total sequencing depth for X-axis (same for all fsizes)
    x = []
    for r in sorted(results_tot[0].keys()):
        h = results_tot[0][r]
        x.append(h.total_reads[0])  # Total reads regardless of family size

    max_y = 0
    for fsize in fsizes:
        y = []
        for r in sorted(results_tot[0].keys()):
            h = results_tot[0][r]
            y.append(h.umis[int(fsize)])
        ax.plot(x, y, "o-", label=f"Family size ≥{fsize}")
        max_y = max(max_y, max(y) if y else 0)

    ax.set_xlabel("Total Sequencing Depth (reads)")
    ax.set_ylabel("Number of UMI Families")
    ax.set_title("Downsampling Saturation Analysis")
    ax.legend(loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, max(x) * 1.05 if x else 1)
    ax.set_ylim(0, max_y * 1.1 if max_y > 0 else 1)

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=150)
    plt.close(fig)


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
        _fraction: List of downsample rates.
        fsizes: List of family sizes to calculate statistics for.
        onlyNamed: If True, only process regions with names.

    Returns:
        List of dictionaries mapping downsample rates to region_cons_stat objects.
    """
    all_results = []
    for h in hist:
        run_analysis = True
        if onlyNamed and h.name == "":
            run_analysis = False
        if run_analysis:
            # family_sizes: list of family sizes (e.g., [5, 3, 10, 2] means 4 families with 5, 3, 10, and 2 reads)
            # singletons: number of singletons
            num_families = len(h.family_sizes)
            # tmpnames: list of family indices (e.g., [0, 1, 2, 3])
            tmpnames = np.array(range(0, num_families))

            # singnames: list of singleton indices (e.g., [3, 4])
            singnames = list(range(num_families, num_families + h.singletons))

            # Convert to np.array
            times = np.array(h.family_sizes)

            # Expand the temp names by the number of reads per family
            # Result: [0,0,0,0,0, 1,1,1, 2,2,2,2,2,2,2,2,2,2]
            # (5 reads from family 0, 3 from family 1, 10 from family 2)
            reads = np.repeat(tmpnames, times, axis=0)  # expand to one entry per read
            reads = list(reads) + singnames
            results = {}
            for r in _fraction:
                # At 50% (r=0.5), we'd sample 10 random reads out of 20 total. Let's say we get:
                # ds_reads = [0, 0, 0, 1, 2, 2, 2, 2, 2, 3]
                ds_reads = random.sample(list(reads), round(r * len(reads)))  # noqa: S311 - downsample
                new_hist = Counter(ds_reads).values()  # collapse to one entry per UMI family
                # Counter: {0: 3, 1: 1, 2: 5, 3: 1}
                # new_hist = [3, 1, 5, 1] (sorted descending: [5, 3, 1, 1])
                new_hist = sorted(new_hist, reverse=True)  # sort
                new_singletons = list(new_hist).count(1)  # count singletons in new
                # Separate consensus families (≥2 reads) from singletons to avoid double-counting
                # RegionConsensusStats expects singletons passed separately to constructor,
                # and only consensus families passed to add_family_sizes()
                consensus_hist = [x for x in new_hist if x > 1]
                new_stat = RegionConsensusStats(h.regionid, h.pos, h.name, new_singletons, h.fsizes)
                new_stat.add_family_sizes(consensus_hist, fsizes)
                results[r] = new_stat
            all_results.append(results)
    return all_results


def run_downsampling(output_path, consensus_bam, bed_file=None, plot_fsizes=None, samplename=None):
    """Run downsampling analysis and generate output files.

    Args:
        output_path: Directory for output files.
        consensus_bam: Path to consensus BAM file (auto-detected if None).
        bed_file: Optional path to BED file for region annotations.
        plot_fsizes: List of family size thresholds to plot (default: [1, 2, 3, 5]).
        samplename: Sample name for output files (auto-detected if None).
    """
    if plot_fsizes is None:
        plot_fsizes = [1, 2, 3, 5]

    logger.info("Running downsampling analysis")
    out_path = Path(output_path)

    if not consensus_bam:
        bam_files = list(out_path.glob("*_consensus_reads.bam"))
        if not bam_files:
            raise FileNotFoundError(f"No consensus BAM file found in {out_path}")
        consensus_bam = str(bam_files[0])
    if not samplename:
        samplename = Path(consensus_bam).name.replace("_consensus_reads.bam", "")

    # Get stats directly from BAM (no stats file needed)
    region_stats_list = get_stat(consensus_bam, bed_file)
    fsizes = list(DEFAULT_FAMILY_SIZES)[1:]  # Exclude 0, which is handled separately
    tot_results = RegionConsensusStats("All", "all_regions", "", 0, fsizes)

    for h in region_stats_list:
        tot_results.family_sizes = tot_results.family_sizes + h.family_sizes
        tot_results.singletons += h.singletons

    downsample_rates = [x * 0.1 for x in range(1, 11)]  # 0.1 to 1.0
    tot = downsample_reads_per_region([tot_results], downsample_rates, fsizes, False)
    all_results = downsample_reads_per_region(region_stats_list, downsample_rates, fsizes, True)

    filename = str(out_path / f"{samplename}_downsampled_coverage.txt")
    save_downsampled_table(all_results, tot, filename)

    filename = str(out_path / f"{samplename}_downsampled_plot.png")
    logger.info(f"Generating downsampling plot with family size thresholds: {plot_fsizes}")
    plot_downsampling(tot, plot_fsizes, filename)
    logger.info(f"Downsampling plot saved to {filename}")
