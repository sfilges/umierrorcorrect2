#!/usr/bin/env python3
from pathlib import Path

import pysam

from umierrorcorrect2.core.constants import DEFAULT_FAMILY_SIZES_STR


def filter_cons(
    filename: str, raw_depth_cutoff: int = 150, fsizes: str = DEFAULT_FAMILY_SIZES_STR, writeraw: bool = False
) -> None:
    """Filter consensus file by depth and family sizes.

    Args:
        filename: Path to input consensus file (tsv).
        raw_depth_cutoff: Minimum raw depth to keep a consensus group.
        fsizes: Comma-separated string of family sizes to keep.
        writeraw: Whether to write raw reads to the output file.
    """
    outfilename = filename.replace("_cons.tsv", "_filtered_cons.tsv")
    fs = fsizes.split(",")
    with Path(filename).open() as f, Path(outfilename).open("w") as g:
        header = f.readline()
        g.write(header)
        passdepth = False
        for line in f:
            parts = line.split("\t")
            if parts[3] not in "":
                fsize = parts[13]
                if fsize == "0":
                    depth = int(parts[12])
                    if depth >= raw_depth_cutoff:
                        if writeraw:
                            g.write(line)
                        passdepth = True
                    else:
                        passdepth = False
                elif passdepth is True and fsize in fs:
                    g.write(line)


def filter_bam(infilename: str, outfilename: str, consensus_cutoff: int) -> None:
    """Filter BAM file by removing reads below consensus depth threshold.

    Args:
        infilename: Path to input BAM file.
        outfilename: Path to output BAM file.
        consensus_cutoff: Minimum consensus family size to keep.
    """
    consensus_cutoff = int(consensus_cutoff)
    with pysam.AlignmentFile(infilename, "rb") as f, pysam.AlignmentFile(outfilename, "wb", template=f) as g:
        reads = f.fetch()
        for read in reads:
            size = int(read.qname.rsplit("=", 1)[-1])  # type: ignore

            if size >= consensus_cutoff:
                g.write(read)
