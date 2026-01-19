#!/usr/bin/env python3
"""Type aliases for umierrorcorrect."""

from typing import TypeAlias

import pysam

# Basic types
Barcode: TypeAlias = str
Position: TypeAlias = int
Contig: TypeAlias = str
RegionId: TypeAlias = str

# Sequence types
Sequence: TypeAlias = str
QualityString: TypeAlias = str
CigarString: TypeAlias = str

# Dictionary types for UMI handling
BarcodeDict: TypeAlias = dict[Barcode, int]
AdjacencyMatrix: TypeAlias = dict[Barcode, list[Barcode]]
ClusterList: TypeAlias = list[list[Barcode]]

# Region types
Region: TypeAlias = tuple[Contig, Position, Position]
BedRegions: TypeAlias = dict[Contig, list[tuple[Position, Position, str]]]

# BAM-related types
AlignedRead: TypeAlias = pysam.AlignedSegment
PositionMatrix: TypeAlias = dict[Barcode, list[AlignedRead]]
SingletonMatrix: TypeAlias = dict[Barcode, AlignedRead]

# Consensus types
ConsensusDict: TypeAlias = dict[Position, dict[str, list[int] | dict[str, int]]]
ConsInfoDict: TypeAlias = dict[Position, dict[int, dict[str, int]]]

# Family size type
FamilySizes: TypeAlias = tuple[int, ...] | list[int]
