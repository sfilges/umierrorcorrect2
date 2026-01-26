#!/usr/bin/env python3
"""Consensus statistics functionality - derives all stats from consensus BAM.

This module provides functions for calculating and reporting statistics
from consensus BAM files. All statistics are derived directly from the BAM,
eliminating the need for separate stats files.
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES, HISTOGRAM_SUFFIX
from umierrorcorrect.core.get_regions_from_bed import read_bed
from umierrorcorrect.core.logging_config import get_logger

logger = get_logger(__name__)


@dataclass
class RegionStats:
    """Statistics for a genomic region, derived from consensus BAM."""

    chrom: str
    start: int
    end: int
    name: str = ""
    consensus_counts: list[int] = field(default_factory=list)
    singleton_count: int = 0

    @property
    def position(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def total_consensus(self) -> int:
        return len(self.consensus_counts)

    @property
    def total_raw_reads(self) -> int:
        """Total original reads (sum of family sizes + singletons)."""
        return sum(self.consensus_counts) + self.singleton_count

    def umis_at_threshold(self, min_family_size: int) -> int:
        """Count UMIs with family size >= threshold."""
        if min_family_size <= 1:
            return len(self.consensus_counts) + self.singleton_count
        return sum(1 for c in self.consensus_counts if c >= min_family_size)

    def reads_at_threshold(self, min_family_size: int) -> int:
        """Count reads from families >= threshold."""
        if min_family_size <= 1:
            return sum(self.consensus_counts) + self.singleton_count
        return sum(c for c in self.consensus_counts if c >= min_family_size)


def parse_consensus_read_name(qname: str) -> tuple[str, int]:
    """Parse consensus BAM read name.

    Examples:
        Consensus_read_0_TTGTAAAGCATGAAATAGC_Count=144
        Consensus_read_0_TTGTAAAGCATGAAATAGC_a_Count=144
        Singleton_read_0_AAGAAAAGCAAGAAATGGT_Count=1

    Returns:
        (read_type, count) where read_type is 'Consensus' or 'Singleton'
    """
    parts = qname.split("_")
    read_type = parts[0]
    count = int(qname.split("=")[-1])
    return read_type, count


def get_stats_from_bam(consensus_bam: str | Path, bed_file: str | Path | None = None) -> list[RegionStats]:
    """Extract all statistics directly from consensus BAM.

    Groups reads by (chrom, start_position) to automatically merge
    regions that were split during processing.

    Args:
        consensus_bam: Path to consensus BAM file.
        bed_file: Optional BED file for region name annotations.

    Returns:
        List of RegionStats objects, sorted by position.
    """
    # Load BED annotations if provided
    bed_regions: dict[str, list[tuple[int, int, str]]] = {}
    if bed_file:
        bed_regions = read_bed(str(bed_file))

    # Key: (chrom, start_pos) -> RegionStats
    regions: dict[tuple[str, int], RegionStats] = {}

    with pysam.AlignmentFile(str(consensus_bam), "rb") as f:
        for read in f.fetch():
            if read.reference_name is None:
                continue

            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end or start + 1

            read_type, count = parse_consensus_read_name(read.query_name)

            key = (chrom, start)
            if key not in regions:
                # Look up region name from BED if available
                name = _get_bed_annotation(bed_regions, chrom, start, end)
                regions[key] = RegionStats(chrom=chrom, start=start, end=end, name=name)

            # Update end position to max seen
            if end > regions[key].end:
                regions[key].end = end

            if read_type == "Consensus":
                regions[key].consensus_counts.append(count)
            else:  # Singleton
                regions[key].singleton_count += 1

    # Sort by chromosome and position
    sorted_regions = sorted(regions.values(), key=lambda r: (r.chrom, r.start))
    return sorted_regions


def _get_bed_annotation(bed_regions: dict, chrom: str, start: int, end: int) -> str:
    """Get region name from BED annotations."""
    if chrom not in bed_regions:
        return ""
    for bed_start, bed_end, name in bed_regions[chrom]:
        if bed_start <= start <= bed_end or bed_start <= end <= bed_end:
            return name
    return ""


class RegionConsensusStats:
    """Statistics for a specific genomic region.

    This class provides backward compatibility with the output format
    while deriving all data from the consensus BAM.

    Attributes:
        regionid (str): Identifier for the region.
        pos (str): Position information.
        name (str): Name of the region.
        singletons (int): Number of singleton reads.
        family_sizes (list): List of family sizes.
        total_reads (dict): Total reads at different family size thresholds.
        umis (dict): UMI counts at different family size thresholds.
        fsizes (list): Family size thresholds.
    """

    def __init__(self, regionid, pos, name, singletons, fsizes):
        """Initialize RegionConsensusStats.

        Args:
            regionid (str): Region identifier.
            pos (str): Position string.
            name (str): Region name.
            singletons (int): Count of singleton reads.
            fsizes (list): List of family size thresholds.
        """
        self.regionid = regionid
        self.pos = pos
        self.name = name
        self.singletons = singletons
        self.family_sizes = []
        self.total_reads = {}
        self.umis = {}
        for fsize in fsizes:
            self.total_reads[fsize] = 0
            self.umis[fsize] = 0
            self.total_reads[0] = self.singletons
            self.umis[0] = self.singletons
        self.total_reads[1] = self.singletons
        self.umis[1] = self.singletons
        self.fsizes = fsizes

    @classmethod
    def from_region_stats(cls, region: RegionStats, fsizes: list[int]) -> "RegionConsensusStats":
        """Create RegionConsensusStats from a RegionStats object.

        Args:
            region: RegionStats object from BAM parsing.
            fsizes: List of family size thresholds.

        Returns:
            RegionConsensusStats with calculated statistics.
        """
        stat = cls(
            regionid=f"{region.chrom}:{region.start}",
            pos=region.position,
            name=region.name,
            singletons=region.singleton_count,
            fsizes=fsizes,
        )
        stat.add_family_sizes(region.consensus_counts, fsizes)
        return stat

    def add_family_sizes(self, family_sizes, fsizes):
        """Add family sizes to the statistics.

        Updates total_reads and umis counts based on family sizes and thresholds.

        Args:
            family_sizes (list): List of family sizes to add.
            fsizes (list): List of family size thresholds.
        """
        # Threshold 0 represents raw read statistics (no UMI deduplication)
        # Total reads and UMIs are both equal to the sum of all reads in the families
        self.total_reads[0] += sum(family_sizes)
        self.umis[0] += sum(family_sizes)

        # Thresholds >= 1 represent unique molecular statistics
        for fsize in fsizes:
            tmp = [x for x in family_sizes if x >= fsize]
            self.total_reads[fsize] += sum(tmp)
            self.umis[fsize] += len(tmp)
        self.family_sizes.extend(family_sizes)

    def write_stats(self):
        """Format statistics for output.

        Returns:
            str: Tab-separated string of statistics for writing to file.
        """
        lines = []
        r0 = self.total_reads[0]
        u0 = self.umis[0]
        line = "\t".join([str(self.regionid), self.pos, self.name, "0", "1.0", str(r0), str(u0)])
        lines.append(line)
        for fsize in self.fsizes:
            fraction = 0 if r0 == 0 else self.total_reads[fsize] / r0
            line = "\t".join(
                [
                    str(self.regionid),
                    self.pos,
                    self.name,
                    str(fsize),
                    str(1.0 * fraction),
                    str(self.total_reads[fsize]),
                    str(self.umis[fsize]),
                ]
            )
            lines.append(line)
        return "\n".join(lines)


def get_stat(consensus_bam: Path | str, bed_file: Path | str | None = None) -> list[RegionConsensusStats]:
    """Get consensus statistics from BAM file.

    This is the main entry point, replacing the old version that needed stats file.
    All statistics are now derived directly from the consensus BAM.

    Args:
        consensus_bam: Path to consensus BAM file.
        bed_file: Optional BED file for region name annotations.

    Returns:
        List of RegionConsensusStats objects.
    """
    fsizes = list(DEFAULT_FAMILY_SIZES)[1:]  # Exclude 0, which is handled separately
    regions = get_stats_from_bam(consensus_bam, bed_file)
    return [RegionConsensusStats.from_region_stats(r, fsizes) for r in regions]


def write_stats_file(stats: list[RegionConsensusStats], output_path: Path, sample_name: str) -> Path:
    """Write statistics to a TSV file for backward compatibility.

    Args:
        stats: List of RegionConsensusStats objects.
        output_path: Output directory.
        sample_name: Sample name for the output file.

    Returns:
        Path to the written stats file.
    """
    stats_file = output_path / f"{sample_name}{HISTOGRAM_SUFFIX}"
    with stats_file.open("w") as f:
        for stat in stats:
            # Write one line per region with consensus and singleton counts
            f.write(
                f"{stat.regionid}\t{stat.pos}\t{stat.name}\t"
                f"consensus_reads: {len(stat.family_sizes)}\t"
                f"singletons: {stat.singletons}\n"
            )
    return stats_file


def calculate_target_coverage(stats, fsizes=None):
    """Calculate target coverage statistics for on-target and off-target reads.

    Args:
        stats (list): List of RegionConsensusStats objects.
        fsizes (list, optional): List of family size thresholds. Defaults to None.

    Returns:
        str: Formatted string of target coverage statistics with columns:
            family_size, on_target, off_target, total, on_target_fraction, off_target_fraction
    """
    if fsizes is None:
        fsizes = stats[0].fsizes if stats else list(DEFAULT_FAMILY_SIZES)[1:]

    reads_all = {}
    reads_target = {}
    reads_offtarget = {}
    fsizes_calc = fsizes.copy()
    fsizes_calc.insert(0, 0)
    for fsize in fsizes_calc:
        reads_all[fsize] = 0
        reads_target[fsize] = 0
        reads_offtarget[fsize] = 0
    for region in stats:
        for fsize in fsizes_calc:
            reads_all[fsize] += region.umis[fsize]
            if region.name:
                reads_target[fsize] += region.umis[fsize]
            else:
                reads_offtarget[fsize] += region.umis[fsize]
    lines = []
    for fsize in fsizes_calc:
        if reads_all[fsize] > 0:
            on_target_frac = reads_target[fsize] / reads_all[fsize]
            off_target_frac = reads_offtarget[fsize] / reads_all[fsize]
            lines.append(
                f"{fsize}\t{reads_target[fsize]}\t{reads_offtarget[fsize]}\t{reads_all[fsize]}\t{on_target_frac}\t{off_target_frac}"
            )
        else:
            lines.append(f"{fsize}\t{reads_target[fsize]}\t{reads_offtarget[fsize]}\t{reads_all[fsize]}\t0\t0")

    return "\n".join(lines)


def get_overall_statistics(hist, fsizes=None):
    """Aggregate statistics across all regions.

    Args:
        hist (list): List of RegionConsensusStats objects.
        fsizes (list, optional): List of family size thresholds. Defaults to None.

    Returns:
        RegionConsensusStats: Aggregated statistics object.
    """
    if fsizes is None:
        fsizes = hist[0].fsizes if hist else list(DEFAULT_FAMILY_SIZES)[1:]

    histall = RegionConsensusStats("All", "all_regions", "", 0, fsizes)
    fsizesnew = fsizes.copy()
    histall.fsizes = fsizes
    fsizesnew.insert(0, 0)
    for fsize in fsizesnew:
        histall.total_reads[fsize] = 0
        histall.umis[fsize] = 0

    for region in hist:
        for fsize in fsizesnew:
            histall.total_reads[fsize] += region.total_reads[fsize]
            histall.umis[fsize] += region.umis[fsize]
    return histall


def run_get_consensus_statistics(
    output_path: str, consensus_filename: str | None, bed_file: str | None, output_raw: bool, samplename: str | None
):
    """Execute the consensus statistics calculation pipeline.

    Args:
        output_path: Directory for output files.
        consensus_filename: Path to the consensus BAM file.
        bed_file: Optional path to BED file for region annotations.
        output_raw: Whether to output raw group counts.
        samplename: Name of the sample.
    """
    logger.info("Getting consensus statistics")
    out_path = Path(output_path)

    if not consensus_filename:
        bam_files = list(out_path.glob("*_consensus_reads.bam"))
        if not bam_files:
            raise FileNotFoundError(f"No consensus BAM file found in {out_path}")
        consensus_filename = str(bam_files[0])

    if not samplename:
        samplename = Path(consensus_filename).name.replace("_consensus_reads.bam", "")

    # Get stats directly from BAM
    hist = get_stat(consensus_filename, bed_file)
    fsizes = list(DEFAULT_FAMILY_SIZES)[1:]  # Exclude 0, which is handled separately
    histall = get_overall_statistics(hist, fsizes)

    # Write summary statistics
    outfilename = out_path / f"{samplename}_summary_statistics.txt"
    logger.info(f"Writing consensus statistics to {outfilename}")
    with outfilename.open("w") as g:
        g.write(histall.write_stats() + "\n")
        for stat in hist:
            g.write(stat.write_stats() + "\n")

    # Write target coverage (only meaningful with a BED file for region annotations)
    if bed_file:
        outfilename = out_path / f"{samplename}_target_coverage.txt"
        with outfilename.open("w") as g:
            g.write(calculate_target_coverage(hist, fsizes))

    # Write raw group counts if requested
    if output_raw:
        largehist = []
        outfilename = out_path / f"{samplename}_consensus_group_counts.txt"
        for h in hist:
            largehist = largehist + h.family_sizes
            largehist = largehist + [1] * h.singletons
        hist_counts = Counter(largehist)
        with outfilename.open("w") as g:
            for size in sorted(hist_counts):
                g.write(str(size) + "\t" + str(hist_counts[size]) + "\n")

    logger.info("Finished consensus statistics")


def main(output_path, consensus_filename, bed_file, output_raw, samplename):
    run_get_consensus_statistics(output_path, consensus_filename, bed_file, output_raw, samplename)
