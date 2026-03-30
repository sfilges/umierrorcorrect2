"""Tests for PostProcessor mutation detection and VAF calculation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from umierrorcorrect2.analysis.post_processor import PostProcessor, _get_alt_allele_count

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CONS_COLS = [
    "Sample Name",
    "Contig",
    "Position",
    "Name",
    "Reference",
    "A",
    "C",
    "G",
    "T",
    "I",
    "D",
    "N",
    "Coverage",
    "Consensus group size",
    "Max Non-ref Allele Count",
    "Max Non-ref Allele Frequency",
    "Max Non-ref Allele",
]


def _make_cons_row(
    sample: str = "S1",
    contig: str = "chr3",
    position: int = 179218304,
    ref: str = "A",
    a: int = 35168,
    c: int = 0,
    g: int = 0,
    t: int = 0,
    i: int = 0,
    d: int = 0,
    n: int = 0,
    fsize: int = 3,
    max_nonref_allele: str = "",
    max_nonref_count: int = 0,
    max_nonref_freq: float = 0.0,
) -> dict:
    coverage = a + c + g + t + i + d
    return {
        "Sample Name": sample,
        "Contig": contig,
        "Position": position,
        "Name": "",
        "Reference": ref,
        "A": a,
        "C": c,
        "G": g,
        "T": t,
        "I": i,
        "D": d,
        "N": n,
        "Coverage": coverage,
        "Consensus group size": fsize,
        "Max Non-ref Allele Count": max_nonref_count,
        "Max Non-ref Allele Frequency": max_nonref_freq,
        "Max Non-ref Allele": max_nonref_allele,
    }


def _make_mutation_bed(tmp_path: Path, lines: list[str]) -> Path:
    bed = tmp_path / "mutations.bed"
    bed.write_text("# chrom\tstart\tend\tname\tref\talt\n" + "\n".join(lines) + "\n")
    return bed


# ---------------------------------------------------------------------------
# _get_alt_allele_count unit tests
# ---------------------------------------------------------------------------


def test_get_alt_allele_count_returns_correct_base():
    row = pd.Series({"expected_alt": "C", "A": 100, "C": 49, "G": 64, "T": 0, "I": 0, "D": 0})
    assert _get_alt_allele_count(row) == 49


def test_get_alt_allele_count_returns_zero_when_no_alt():
    row = pd.Series({"expected_alt": "C", "A": 100, "C": 0, "G": 64, "T": 0, "I": 0, "D": 0})
    assert _get_alt_allele_count(row) == 0


def test_get_alt_allele_count_returns_zero_for_nan_expected_alt():
    row = pd.Series({"expected_alt": float("nan"), "A": 100, "C": 49, "G": 64})
    assert _get_alt_allele_count(row) == 0


def test_get_alt_allele_count_unknown_alt_returns_zero():
    row = pd.Series({"expected_alt": "X", "A": 100, "C": 49})
    assert _get_alt_allele_count(row) == 0


# ---------------------------------------------------------------------------
# _join_mutations integration tests
# ---------------------------------------------------------------------------


def test_two_mutations_same_locus_both_detected(tmp_path):
    """E545A (→C) and E545G (→G) at the same position must both be alt_matches=True."""
    # Position 179218304: ref=A, C=49, G=64 — G is the max non-ref, but C also present
    row = _make_cons_row(a=35168, c=49, g=64, max_nonref_allele="G", max_nonref_count=64)
    cons_df = pd.DataFrame([row, row])  # same position appears twice (one per mutation)

    bed = _make_mutation_bed(
        tmp_path,
        [
            "chr3\t179218304\t179218304\tPIK3CA_p.E545A\tA\tC",
            "chr3\t179218304\t179218304\tPIK3CA_p.E545G\tA\tG",
        ],
    )

    pp = PostProcessor(tmp_path)
    result = pp._join_mutations(cons_df, bed)

    e545a = result[result["mutation_name"] == "PIK3CA_p.E545A"].iloc[0]
    e545g = result[result["mutation_name"] == "PIK3CA_p.E545G"].iloc[0]

    assert e545a["alt_matches"]
    assert e545a["alt_allele_count"] == 49

    assert e545g["alt_matches"]
    assert e545g["alt_allele_count"] == 64


def test_mutation_detected_even_when_not_max_nonref(tmp_path):
    """If expected alt has reads but is not the max non-ref allele, it should still be detected."""
    row = _make_cons_row(a=1000, c=5, g=50, max_nonref_allele="G", max_nonref_count=50)
    cons_df = pd.DataFrame([row])

    bed = _make_mutation_bed(
        tmp_path,
        [
            "chr3\t179218304\t179218304\tMUT1\tA\tC",
        ],
    )

    pp = PostProcessor(tmp_path)
    result = pp._join_mutations(cons_df, bed)

    assert result.iloc[0]["alt_matches"]
    assert result.iloc[0]["alt_allele_count"] == 5


def test_mutation_not_detected_when_no_alt_reads(tmp_path):
    """If the expected alt has zero reads, alt_matches must be False."""
    row = _make_cons_row(a=1000, c=0, g=0, max_nonref_allele="", max_nonref_count=0)
    cons_df = pd.DataFrame([row])

    bed = _make_mutation_bed(
        tmp_path,
        [
            "chr3\t179218304\t179218304\tMUT1\tA\tC",
        ],
    )

    pp = PostProcessor(tmp_path)
    result = pp._join_mutations(cons_df, bed)

    assert not result.iloc[0]["alt_matches"]
    assert result.iloc[0]["alt_allele_count"] == 0


def test_normal_single_mutation_detected(tmp_path):
    """Standard case: single mutation, expected alt is the max non-ref allele."""
    row = _make_cons_row(a=30000, c=0, g=77, max_nonref_allele="G", max_nonref_count=77)
    cons_df = pd.DataFrame([row])

    bed = _make_mutation_bed(
        tmp_path,
        [
            "chr3\t179218304\t179218304\tPIK3CA_p.Q546E\tA\tG",
        ],
    )

    pp = PostProcessor(tmp_path)
    result = pp._join_mutations(cons_df, bed)

    assert result.iloc[0]["alt_matches"]
    assert result.iloc[0]["alt_allele_count"] == 77


# ---------------------------------------------------------------------------
# compute_mutation_metrics VAF correctness
# ---------------------------------------------------------------------------


def _make_joined_df(alt_allele_count: int, coverage: int, alt_matches: bool) -> pd.DataFrame:
    """Build a minimal post-join DataFrame for compute_mutation_metrics."""
    return pd.DataFrame(
        [
            {
                "Sample Name": "S1",
                "Contig": "chr3",
                "Position": 179218304,
                "Reference": "A",
                "A": coverage - alt_allele_count,
                "C": 0,
                "G": 0,
                "T": 0,
                "I": 0,
                "D": 0,
                "N": 0,
                "Coverage": coverage,
                "Consensus group size": 3,
                "Max Non-ref Allele Count": alt_allele_count,
                "Max Non-ref Allele Frequency": alt_allele_count / coverage if coverage else 0.0,
                "Max Non-ref Allele": "C",
                "mutation_name": "MUT1",
                "expected_alt": "C",
                "is_mutation": True,
                "alt_allele_count": alt_allele_count,
                "alt_matches": alt_matches,
            }
        ]
    )


def test_vaf_uses_alt_allele_count(tmp_path):
    """VAF must be computed from alt_allele_count, not Max Non-ref Allele Count."""
    # alt_allele_count=49, Coverage=35281
    df = _make_joined_df(alt_allele_count=49, coverage=35281, alt_matches=True)
    # Deliberately set Max Non-ref Allele Count to a different value to confirm independence
    df["Max Non-ref Allele Count"] = 999

    pp = PostProcessor(tmp_path)
    result = pp.compute_mutation_metrics(df)

    expected_vaf = 49 / 35281 * 100
    assert abs(result.iloc[0]["VAF (%)"] - expected_vaf) < 1e-6
    assert abs(result.iloc[0]["ctDNA ppm"] - 49 / 35281 * 1e6) < 1e-3


def test_vaf_is_zero_when_no_alt_reads(tmp_path):
    """VAF must be 0.0 when alt_allele_count=0."""
    df = _make_joined_df(alt_allele_count=0, coverage=35000, alt_matches=False)

    pp = PostProcessor(tmp_path)
    result = pp.compute_mutation_metrics(df)

    assert result.iloc[0]["VAF (%)"] == 0.0
    assert result.iloc[0]["ctDNA ppm"] == 0.0


def test_two_mutations_same_locus_correct_vafs(tmp_path):
    """End-to-end: E545A and E545G at same position get independent correct VAFs."""
    # Build joined DataFrame as _join_mutations would produce
    base = {
        "Sample Name": "S1",
        "Contig": "chr3",
        "Position": 179218304,
        "Reference": "A",
        "A": 35168,
        "C": 49,
        "G": 64,
        "T": 0,
        "I": 0,
        "D": 0,
        "N": 0,
        "Coverage": 35281,
        "Consensus group size": 3,
        "Max Non-ref Allele Count": 64,
        "Max Non-ref Allele Frequency": 64 / 35281,
        "Max Non-ref Allele": "G",
        "is_mutation": True,
    }
    rows = [
        {**base, "mutation_name": "PIK3CA_p.E545A", "expected_alt": "C", "alt_allele_count": 49, "alt_matches": True},
        {**base, "mutation_name": "PIK3CA_p.E545G", "expected_alt": "G", "alt_allele_count": 64, "alt_matches": True},
    ]
    df = pd.DataFrame(rows)

    pp = PostProcessor(tmp_path)
    result = pp.compute_mutation_metrics(df)

    e545a = result[result["mutation_name"] == "PIK3CA_p.E545A"].iloc[0]
    e545g = result[result["mutation_name"] == "PIK3CA_p.E545G"].iloc[0]

    assert abs(e545a["VAF (%)"] - 49 / 35281 * 100) < 1e-6
    assert abs(e545g["VAF (%)"] - 64 / 35281 * 100) < 1e-6
