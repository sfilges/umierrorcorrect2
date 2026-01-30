#!/usr/bin/env python3
import warnings
from pathlib import Path

import pysam

# Type alias for region dictionary
RegionDict = dict[str, list[tuple[int, int, str]]]


def read_bed(bedfile: str) -> RegionDict:
    """Read a BED file and return regions organized by contig."""
    regions: RegionDict = {}
    with Path(bedfile).open() as f:
        for line in f:
            line = line.strip()
            parts = line.split()
            if len(parts) >= 4:
                contig, start, end, name = parts[0:4]
                if contig not in regions:
                    regions[contig] = []
                regions[contig].append((int(start), int(end), name))
    return regions


def sort_regions(regions: RegionDict) -> RegionDict:
    """Sort regions by start position within each contig."""
    newregions: RegionDict = {}
    for contig in regions:
        newregions[contig] = sorted(regions[contig], key=lambda tup: tup[0])
    return newregions


def merge_regions(regions: RegionDict, pos_threshold: int) -> RegionDict:
    """Merge overlapping or adjacent regions within threshold distance."""
    newregions: RegionDict = {}
    for contig in regions:
        newregions[contig] = []
        starts, ends, names = zip(*regions[contig])
        current_start = starts[0]
        current_end = ends[0]
        current_name: list[str] = []
        current_name.append(names[0])

        for current_index, start in enumerate(starts):
            if current_end + pos_threshold < start:
                newregions[contig].append(
                    (current_start - pos_threshold, current_end + pos_threshold, ",".join(current_name))
                )
                current_start = start
                current_name = []
            current_end = ends[current_index]
            if names[current_index] not in current_name:
                current_name.append(names[current_index])

        newregions[contig].append(
            (current_start - pos_threshold, current_end + pos_threshold, ",".join(current_name))
        )  # save last entry
    return newregions


def get_first_annotation(regions: list[tuple[int, int, str]], pos: int) -> str:
    """Get first matching annotation for a position from a list of regions.

    .. deprecated::
        This function is not used in production code. Consider using
        `get_all_annotations` instead which returns all matching annotations.

    Args:
        regions: List of (start, end, name) tuples defining annotated regions.
        pos: Position to look up.

    Returns:
        Name of the first matching region, or empty string if no match.
    """
    warnings.warn(
        "get_first_annotation is deprecated and unused. Use get_all_annotations instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    for start, end, name in regions:
        if pos >= start and pos <= end:
            return name
    return ""


def get_all_annotations(regions: list[tuple[int, int, str]], pos: int) -> str:
    """Get all annotations for a position (comma-separated if multiple).

    Args:
        regions: List of (start, end, name) tuples defining annotated regions.
        pos: Position to look up.

    Returns:
        Comma-separated string of all matching region names, or empty string if no match.
    """
    annotation: list[str] = []
    for start, end, name in regions:
        if pos >= start and pos <= end:
            annotation.append(name)
    return ",".join(annotation)


def get_overlap(annotation_regions: list[tuple[int, int, str]], start: int, end: int) -> str:
    """Check if a region overlaps with any annotation region."""
    for annotation_start, annotation_end, annotation_name in annotation_regions:
        if annotation_start <= end and start <= annotation_end:  # test for overlap
            return annotation_name
    return ""


def expand_regions_from_bed(regions: RegionDict, bamfile: str) -> RegionDict:
    """Expand regions based on actual read positions in BAM file."""
    newregions: RegionDict = {}
    with pysam.AlignmentFile(bamfile, "rb") as f:
        for contig in regions:
            newregions[contig] = []
            for annotation_start, annotation_end, annotation_name in regions[contig]:
                minpos = annotation_start
                maxpos = annotation_end
                reads = f.pileup(contig, annotation_start, annotation_end)
                for pileupColumn in reads:
                    for r in pileupColumn.pileups:
                        pos = r.alignment.reference_start
                        if pos < minpos:
                            minpos = pos
                        if pos > maxpos:
                            maxpos = pos
                newregions[contig].append((minpos, maxpos, annotation_name))
    return newregions
