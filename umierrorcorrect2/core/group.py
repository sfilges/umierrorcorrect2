#!/usr/bin/env python3
from collections import Counter

import pysam

from umierrorcorrect2.core.get_regions_from_bed import expand_regions_from_bed, merge_regions, read_bed, sort_regions


def get_chromosome_list_from_bam(f: pysam.AlignmentFile) -> list[str]:
    """Get list of chromosomes (contigs) from BAM file that have mapped reads."""
    contiglist = []
    for chrx in f.get_index_statistics():
        if chrx.total > 0:
            contiglist.append(chrx.contig)
    return contiglist


def group_by_position(
    f: pysam.AlignmentFile, chrx: str, pos_threshold: int
) -> tuple[dict[int, Counter[str]], dict[int, int]]:
    """
    Group reads by position and barcode.

    Args:
        f: Open pysam AlignmentFile.
        chrx: Chromosome/contig name.
        pos_threshold: Distance threshold for grouping.

    Returns:
        A tuple (regions, ends):
            - regions: maps start_position -> barcode -> count.
            - ends: maps start_position -> end_position.
    """
    regions: dict[int, Counter[str]] = {}
    ends: dict[int, int] = {}
    reads = f.fetch(chrx)
    current_pos = -pos_threshold
    current_end = -(2 * pos_threshold) + 1
    current_region_qnames: dict[str, set[str]] = {}  # barcode -> set of qnames to treat mate pairs as one fragment

    for line in reads:
        pos = line.pos
        if pos > current_end + pos_threshold:
            # finalize old region
            if current_pos != -pos_threshold:
                ends[current_pos] = current_end + 1
                regions[current_pos] = Counter()
                # Count unique fragments (QNAMEs) per barcode
                for bc, qnames in current_region_qnames.items():
                    regions[current_pos][bc] = len(qnames)
                current_region_qnames = {}

            current_pos = pos
            current_end = pos
            current_region_qnames = {}

        barcode = line.qname.rsplit(":", 1)[-1]
        if barcode not in current_region_qnames:
            current_region_qnames[barcode] = set()
        # Add query name to set to ensure R1 and R2 from same fragment only count once
        current_region_qnames[barcode].add(line.qname)
        current_end = pos

    # finalize last region
    if current_pos != -pos_threshold:
        ends[current_pos] = current_end + 1
        regions[current_pos] = Counter()
        for bc, qnames in current_region_qnames.items():
            regions[current_pos][bc] = len(qnames)

    return (regions, ends)


def count_umis_in_region(f: pysam.AlignmentFile, chrx: str, pos_start: int, pos_end: int) -> Counter[str]:
    """
    Count UMIs (barcodes) in a specific genomic region.

    Args:
        f: Open pysam AlignmentFile.
        chrx: Chromosome/contig name.
        pos_start: Start position (0-based).
        pos_end: End position (0-based).

    Returns:
        Counter mapping barcode -> count of unique fragments.
    """
    region_qnames: dict[str, set[str]] = {}  # barcode -> set of qnames to treat mate pairs as one fragment
    reads = f.fetch(chrx, pos_start, pos_end)
    for line in reads:
        # pos = line.pos
        barcode = line.qname.rsplit(":", 1)[-1]
        if barcode not in region_qnames:
            region_qnames[barcode] = set()
        # Add query name to set to ensure R1 and R2 from same fragment only count once
        region_qnames[barcode].add(line.qname)

    region: Counter[str] = Counter()
    for barcode in region_qnames:
        # Final family size is the number of unique original fragments (QNAMEs)
        region[barcode] = len(region_qnames[barcode])
    return region


def get_max_number_of_barcodes(regions: dict[int, Counter[str]], pos: int) -> int:
    """Get the count of the most common barcode at a specific position."""
    if len(regions[pos]) > 0:
        umi, count = regions[pos].most_common(1)[0]
        return count
    else:
        return 0


def remove_singleton_regions(
    regions: dict[str, dict[int, Counter[str]]], cutoff: int
) -> dict[str, dict[int, Counter[str]]]:
    """
    Remove regions where the most common barcode has fewer reads than cutoff.

    Args:
        regions: dict mapping chromosome -> position -> barcode -> count.
        cutoff: Minimum count threshold.

    Returns:
        Filtered regions dictionary.
    """
    newregions: dict[str, dict[int, Counter[str]]] = {}
    for chrx in regions:
        newregions[chrx] = {
            x: y for x, y in regions[chrx].items() if get_max_number_of_barcodes(regions[chrx], x) > cutoff
        }
    return newregions


def readBam(
    infile: str, position_threshold: int
) -> tuple[dict[str, dict[int, Counter[str]]], dict[str, dict[int, int]]]:
    """
    Read grouped BAM-file (UMI-groups-sorted) and extract sequences from each UMI-group.

    Saves one representative read from each group.

    Args:
        infile: Path to the input BAM file.
        position_threshold: Threshold for position grouping.

    Returns:
        Tuple (regions, chrends):
            - regions: dict mapping contig -> start_pos -> barcode -> count.
            - chrends: dict mapping contig -> start_pos -> end_pos.
    """
    with pysam.AlignmentFile(infile, "rb") as f:
        chrs = get_chromosome_list_from_bam(f)
        chrregions: dict[str, dict[int, Counter[str]]] = {}
        chrends: dict[str, dict[int, int]] = {}
        for chrx in chrs:
            regions, ends = group_by_position(f, chrx, position_threshold)
            chrregions[chrx] = regions
            chrends[chrx] = ends
            # chrregions[chrx]=group_by_position(f,chrx,position_threshold)

        regions_filtered = remove_singleton_regions(chrregions, 2)
        # for chrx in chrs:
        #     print(chrx,regions[chrx])

        # for chrx in chrs:
        #     regions2=regions[chrx]
        #     for rr in regions2:
        #         print(chrx,rr,regions2[rr].most_common(10))
    return (regions_filtered, chrends)


def read_bam_from_bed(
    infile: str, bedfile: str, position_threshold: int
) -> tuple[dict[str, dict[int, Counter[str]]], dict[str, dict[int, int]]]:
    """
    Read BAM file restricted to regions defined in a BED file.

    Args:
        infile: Path to BAM file.
        bedfile: Path to BED file.
        position_threshold: merging threshold.

    Returns:
        Tuple (regions, chrends) similar to readBam.
    """
    chrregions: dict[str, dict[int, Counter[str]]] = {}
    chrends: dict[str, dict[int, int]] = {}
    regions_bed = read_bed(bedfile)
    regions_bed = sort_regions(regions_bed)
    regions_bed = merge_regions(regions_bed, position_threshold)
    regions_bed = expand_regions_from_bed(regions_bed, infile)
    newregions = []
    for contig in regions_bed:
        for a, b, c in regions_bed[contig]:
            newregions.append((contig, a, b, c))

    with pysam.AlignmentFile(infile, "rb") as f:
        chrs = get_chromosome_list_from_bam(f)
        for contig, start, end, _name in newregions:
            if contig in chrs:
                if contig not in chrregions:
                    chrregions[contig] = {}
                chrregions[contig][start] = count_umis_in_region(f, contig, start, end)
                if contig not in chrends:
                    chrends[contig] = {}
                chrends[contig][start] = end
        regions_filtered = remove_singleton_regions(chrregions, 2)
    return (regions_filtered, chrends)


def read_bam_from_tag(
    infile: str,
) -> tuple[dict[str, dict[str, Counter[str]]], dict[str, dict[str, int]], dict[str, dict[str, int]]]:
    """
    Read BAM file groupings from the 'UG' tag.

    Args:
        infile: Path to BAM file.

    Returns:
        Tuple (regions, starts, ends):
            - regions: dict[contig][utag] -> Counter[barcode]
            - starts: dict[contig][utag] -> start_pos
            - ends: dict[contig][utag] -> end_pos
    """
    regions_qnames: dict[str, dict[str, dict[str, set[str]]]] = {}  # barcode -> set of qnames
    starts: dict[str, dict[str, int]] = {}
    ends: dict[str, dict[str, int]] = {}
    with pysam.AlignmentFile(infile, "rb") as f:
        reads = f.fetch()
        for r in reads:
            contig = r.reference_name
            if contig not in regions_qnames:
                regions_qnames[contig] = {}
                starts[contig] = {}
                ends[contig] = {}
            try:
                utag = r.get_tag("UG")
            except KeyError:
                print("UG tag not present in BAM file")
                utag = ""
            if utag not in regions_qnames[contig]:
                regions_qnames[contig][utag] = {}
                starts[contig][utag] = r.reference_start
                ends[contig][utag] = r.reference_end
            barcode = r.qname.rsplit(":", 1)[-1]
            if barcode not in regions_qnames[contig][utag]:
                regions_qnames[contig][utag][barcode] = set()
            # Add query name to set to ensure R1 and R2 from same fragment only count once
            regions_qnames[contig][utag][barcode].add(r.qname)

            if r.reference_start < starts[contig][utag]:
                starts[contig][utag] = r.reference_start
            if r.reference_end > ends[contig][utag]:
                ends[contig][utag] = r.reference_end

    regions: dict[str, dict[str, Counter[str]]] = {}
    for contig in regions_qnames:
        regions[contig] = {}
        for utag in regions_qnames[contig]:
            regions[contig][utag] = Counter()
            # Final family size is the number of unique original fragments (QNAMEs)
            for barcode in regions_qnames[contig][utag]:
                regions[contig][utag][barcode] = len(regions_qnames[contig][utag][barcode])
    return (regions, starts, ends)
