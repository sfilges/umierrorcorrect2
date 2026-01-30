#!/usr/bin/env python3
"""Constants and type aliases used throughout the umierrorcorrect package."""

from typing import TypeAlias

import pysam

# =============================================================================
# Type Aliases
# =============================================================================
# TODO: Currently these are not used anywhere. Should they be removed?
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

# =============================================================================
# Phred Score Constants
# =============================================================================
MAX_PHRED_SCORE = 60
"""Maximum phred quality score to assign to consensus bases."""

PHRED_TABLE_SIZE = 94
"""Size of phred lookup tables (ASCII 33-126)."""

DEFAULT_MAPPING_QUALITY = 60
"""Default mapping quality for consensus reads."""

# =============================================================================
# Consensus Generation Constants
# =============================================================================
COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS = 50
"""Coverage threshold above which simple majority voting is used instead of
probability-based consensus calculation."""

DEFAULT_INDEL_FREQUENCY_THRESHOLD = 60.0
"""Default frequency threshold (%) for including indels in consensus."""

DEFAULT_CONSENSUS_FREQUENCY_THRESHOLD = 60.0
"""Default frequency threshold (%) for consensus base calls."""

# =============================================================================
# UMI Clustering Constants
# =============================================================================
DEFAULT_EDIT_DISTANCE_THRESHOLD = 1
"""Default edit distance threshold for UMI clustering."""

SUBSTRING_OPTIMIZATION_THRESHOLD = 30
"""Number of UMIs above which substring-based optimization is used for clustering."""

# =============================================================================
# Family Size Constants
# =============================================================================
DEFAULT_FAMILY_SIZES = (0, 1, 2, 3, 4, 5, 7, 10, 20, 30)
"""Default family size thresholds for consensus statistics."""

SINGLETON_FAMILY_SIZES = (0, 1)
"""Family sizes used for singleton reads."""

# =============================================================================
# File Suffixes
# =============================================================================
CONSENSUS_BAM_SUFFIX = "_consensus_reads.bam"
"""Suffix for consensus BAM output files."""

CONSENSUS_TSV_SUFFIX = "_cons.tsv"
"""Suffix for consensus TSV output files."""

HISTOGRAM_SUFFIX = "_consensus_stats.tsv"
"""Suffix for consensus statistics files (previously .hist)."""

DEFAULT_FAMILY_SIZES_STR = ",".join(str(x) for x in DEFAULT_FAMILY_SIZES)
"""String version of DEFAULT_FAMILY_SIZES for CLI defaults."""

# =============================================================================
# BAM Tags
# =============================================================================
READ_GROUP_TAG = "L1"
"""Default read group tag for consensus reads."""

NM_TAG = "NM"
"""Edit distance tag name."""

RG_TAG = "RG"
"""Read group tag name."""

ASCII_ART = r"""
▗▖ ▗▖▗▖  ▗▖▗▄▄▄▖    ▗▄▄▄▖▗▄▄▖ ▗▄▄▖  ▗▄▖ ▗▄▄▖      ▗▄▄▖ ▗▄▖ ▗▄▄▖ ▗▄▄▖ ▗▄▄▄▖ ▗▄▄▖▗▄▄▄▖
▐▌ ▐▌▐▛▚▞▜▌  █      ▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌    ▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▌   ▐▌     █
▐▌ ▐▌▐▌  ▐▌  █      ▐▛▀▀▘▐▛▀▚▖▐▛▀▚▖▐▌ ▐▌▐▛▀▚▖    ▐▌   ▐▌ ▐▌▐▛▀▚▖▐▛▀▚▖▐▛▀▀▘▐▌     █
▝▚▄▞▘▐▌  ▐▌▗▄█▄▖    ▐▙▄▄▖▐▌ ▐▌▐▌ ▐▌▝▚▄▞▘▐▌ ▐▌    ▝▚▄▄▖▝▚▄▞▘▐▌ ▐▌▐▌ ▐▌▐▙▄▄▖▝▚▄▄▖  █
"""
"""ASCII art for the CLI banner."""

ASCII_ART_2 = r"""
▌ ▌▙▗▌▜▘ ▛▀▘             ▞▀▖               ▐
▌ ▌▌▘▌▐  ▙▄ ▙▀▖▙▀▖▞▀▖▙▀▖ ▌  ▞▀▖▙▀▖▙▀▖▞▀▖▞▀▖▜▀
▌ ▌▌ ▌▐  ▌  ▌  ▌  ▌ ▌▌   ▌ ▖▌ ▌▌  ▌  ▛▀ ▌ ▖▐ ▖
▝▀ ▘ ▘▀▘ ▀▀▘▘  ▘  ▝▀ ▘   ▝▀ ▝▀ ▘  ▘  ▝▀▘▝▀  ▀
"""
"""Alternative ASCII art for the CLI banner."""
