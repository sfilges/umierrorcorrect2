from pathlib import Path

import pytest
from umierrorcorrect2.analysis.models import load_mutations


def test_load_mutations_valid(tmp_path):
    """Test loading a valid mutation BED file."""
    bed_content = (
        "# chrom  start  end  name  ref  alt\n"
        "chr3\t179218303\t179218303\tPIK3CA_p.E545K\tG\tA\n"
        "chr17\t7577538\t7577538\tTP53_R248W\tC\tT\n"
    )
    bed_path = tmp_path / "mutations.bed"
    bed_path.write_text(bed_content)

    mutations = load_mutations(bed_path)

    assert len(mutations) == 2
    assert mutations[0].chromosome == "chr3"
    assert mutations[0].position == 179218303
    assert mutations[0].name == "PIK3CA_p.E545K"
    assert mutations[0].reference == "G"
    assert mutations[0].alternate == "A"

    assert mutations[1].name == "TP53_R248W"


def test_load_mutations_with_comments_and_empty_lines(tmp_path):
    """Test that comments and empty lines are ignored."""
    bed_content = "\n# This is a comment\nchr3\t100\t100\tMUT1\tA\tT\n   \nchr4\t200\t200\tMUT2\tG\tC\n"
    bed_path = tmp_path / "test.bed"
    bed_path.write_text(bed_content)

    mutations = load_mutations(bed_path)
    assert len(mutations) == 2
    assert mutations[0].name == "MUT1"
    assert mutations[1].name == "MUT2"


def test_load_mutations_invalid_format(tmp_path):
    """Test that malformed lines raise ValueError."""
    bed_content = "chr3\t100\t100\tONLY_FOUR_FIELDS"
    bed_path = tmp_path / "bad.bed"
    bed_path.write_text(bed_content)

    with pytest.raises(ValueError, match="requires 6 fields"):
        load_mutations(bed_path)


def test_load_mutations_with_extra_columns(tmp_path):
    """Test that extra columns (beyond the required 6) are ignored."""
    bed_content = "chr3\t100\t100\tMUT1\tA\tT\tEXTRA_COL1\tEXTRA_COL2\n"
    bed_path = tmp_path / "extra.bed"
    bed_path.write_text(bed_content)

    mutations = load_mutations(bed_path)
    assert len(mutations) == 1
    assert mutations[0].name == "MUT1"
    assert mutations[0].alternate == "T"


def test_load_mutations_nonexistent_file():
    """Test that loading a non-existent file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_mutations(Path("nonexistent.bed"))
