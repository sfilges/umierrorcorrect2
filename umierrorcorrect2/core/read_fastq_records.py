#!/usr/bin/env python3

from collections.abc import Generator
from typing import TextIO


def read_fastq(infile: TextIO) -> Generator[tuple[str, str, str], None, None]:
    """Read one fastq record at a time using a generator.

    Args:
        infile: Open file handle for reading FASTQ records.

    Yields:
        Tuple of (name, sequence, quality) for each record.
    """
    for line in infile:
        name = line.rstrip()
        seq = infile.readline().rstrip()
        infile.readline()  # skip + line
        qual = infile.readline().rstrip()
        yield (name, seq, qual)


def read_fastq_paired_end(r1file: TextIO, r2file: TextIO) -> Generator[tuple[str, str, str, str, str, str], None, None]:
    """Read paired-end FASTQ records from two files.

    Args:
        r1file: Open file handle for R1 FASTQ file.
        r2file: Open file handle for R2 FASTQ file.

    Yields:
        Tuple of (name1, seq1, qual1, name2, seq2, qual2) for each read pair.
    """
    for line1, line2 in zip(r1file, r2file):
        name1 = line1.rstrip()
        name2 = line2.rstrip()
        seq1 = r1file.readline().rstrip()
        seq2 = r2file.readline().rstrip()
        r1file.readline()
        r2file.readline()
        qual1 = r1file.readline().rstrip()
        qual2 = r2file.readline().rstrip()
        yield (name1, seq1, qual1, name2, seq2, qual2)
