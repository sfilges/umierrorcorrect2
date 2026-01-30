"""Shared pytest fixtures for umierrorcorrect tests."""

import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir():
    """Return path to the main test_data directory."""
    return Path(__file__).parent / "test_data"


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)


@pytest.fixture
def sample_barcode_counts():
    """Sample UMI barcode counts for clustering tests.

    Contains:
    - ACGTACGTACGT: 10 reads (most abundant)
    - ACGTACGTACGA: 2 reads (1 edit distance from first - should cluster)
    - ACGTACGTACGG: 3 reads (1 edit distance from first - should cluster)
    - GGGGGGGGGGGG: 5 reads (very different - should not cluster)
    - GGGGGGGGGGGA: 1 read (1 edit distance from GGGG... - should cluster)
    """
    return {
        "ACGTACGTACGT": 10,
        "ACGTACGTACGA": 2,
        "ACGTACGTACGG": 3,
        "GGGGGGGGGGGG": 5,
        "GGGGGGGGGGGA": 1,
    }


@pytest.fixture
def sample_barcode_counts_large():
    """Larger barcode set to test substring optimization (>30 barcodes)."""
    import random

    random.seed(42)
    bases = "ACGT"
    barcodes = {}

    # Generate 50 random 12-mer barcodes
    for _ in range(50):
        barcode = "".join(random.choice(bases) for _ in range(12))  # noqa: S311
        barcodes[barcode] = random.randint(1, 100)  # noqa: S311

    # Add some known clusters (1 edit distance apart)
    barcodes["AAAAAAAAAAAA"] = 50
    barcodes["AAAAAAAAAAAC"] = 10
    barcodes["AAAAAAAAAAAG"] = 5

    return barcodes


@pytest.fixture
def sample_bed_regions():
    """Sample BED regions for testing.

    Returns a dict with contig as key and list of (start, end, name) tuples.
    """
    return {
        "chr1": [
            (100, 200, "region1"),
            (300, 400, "region2"),
            (450, 550, "region3"),
        ],
        "chr2": [
            (1000, 1100, "region4"),
        ],
    }


@pytest.fixture
def sample_bed_regions_overlapping():
    """Sample BED regions that overlap or are adjacent."""
    return {
        "chr1": [
            (100, 200, "region1"),
            (190, 300, "region2"),  # Overlaps with region1
            (500, 600, "region3"),
        ],
    }


@pytest.fixture
def sample_bed_regions_unsorted():
    """Sample BED regions that are not sorted by position."""
    return {
        "chr1": [
            (300, 400, "region2"),
            (100, 200, "region1"),
            (450, 550, "region3"),
        ],
    }


@pytest.fixture
def temp_bed_file(temp_output_dir):
    """Create a temporary BED file for testing."""
    bed_path = temp_output_dir / "test_regions.bed"
    bed_content = """chr1\t100\t200\tregion1
chr1\t300\t400\tregion2
chr2\t1000\t1100\tregion3
"""
    bed_path.write_text(bed_content)
    return bed_path


@pytest.fixture
def sample_consensus_positions():
    """Sample position data for consensus testing.

    Each position is a dict mapping nucleotides to lists of phred scores.
    """
    return [
        # Position 0: Strong A consensus
        {"A": [30, 30, 30, 30], "C": [], "G": [], "T": []},
        # Position 1: Mixed with A winning
        {"A": [30, 30, 25], "C": [20], "G": [], "T": []},
        # Position 2: Strong T consensus
        {"A": [], "C": [], "G": [], "T": [35, 35, 35, 35]},
        # Position 3: Low quality, ambiguous
        {"A": [10, 10], "C": [10, 10], "G": [], "T": []},
    ]


@pytest.fixture
def sample_umi_groups():
    """Sample UMI groups with read counts for consensus testing."""
    return {
        "ACGTACGTACGT": 10,
        "TGCATGCATGCA": 5,
        "AAAACCCCGGGG": 3,
    }
