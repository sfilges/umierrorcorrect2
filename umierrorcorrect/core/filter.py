#!/usr/bin/env python3
from pathlib import Path

import pysam

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES_STR


def filter_cons(filename, raw_depth_cutoff=150, fsizes=DEFAULT_FAMILY_SIZES_STR, writeraw=False):
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


def filter_bam(infilename, outfilename, consensus_cutoff):
    consensus_cutoff = int(consensus_cutoff)
    with pysam.AlignmentFile(infilename, "rb") as f, pysam.AlignmentFile(outfilename, "wb", template=f) as g:
        reads = f.fetch()
        for read in reads:
            size = int(read.qname.rsplit("=", 1)[-1])
            if size >= consensus_cutoff:
                g.write(read)
