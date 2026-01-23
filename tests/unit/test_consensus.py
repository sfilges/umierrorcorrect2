"""Unit tests for umierrorcorrect.core.consensus module."""

import pytest
from umierrorcorrect.core.consensus import (
    ConsensusRead,
    _phred_to_error,
    _phred_to_prob,
    calc_consensus,
    calc_consensus_probabilities,
    get_ascii,
    get_most_common_allele,
    get_phred,
    get_position_coverage,
)


class TestPhredConversions:
    """Tests for phred score conversion functions."""

    def test_get_phred_basic(self):
        """Test basic phred score extraction from ASCII."""
        # ASCII '!' is 33, which is phred 0
        assert get_phred("!") == 0
        # ASCII 'I' is 73, which is phred 40
        assert get_phred("I") == 40
        # ASCII '~' is 126, which is phred 93
        assert get_phred("~") == 93

    def test_get_ascii_basic(self):
        """Test basic ASCII character generation from phred score."""
        # Phred 0 should give '!'
        assert get_ascii(0) == "!"
        # Phred 40 should give 'I'
        assert get_ascii(40) == "I"
        # Phred 93 should give '~'
        assert get_ascii(93) == "~"

    def test_phred_ascii_roundtrip(self):
        """Test that phred -> ascii -> phred roundtrips correctly."""
        for phred in range(0, 94):
            assert get_phred(get_ascii(phred)) == phred

    def test_phred_to_prob_values(self):
        """Test phred to probability of correctness conversion."""
        # Phred 0: 10^0 = 1, so error = 1, prob = 0
        assert _phred_to_prob(0) == pytest.approx(0.0, abs=1e-10)

        # Phred 10: error = 0.1, prob = 0.9
        assert _phred_to_prob(10) == pytest.approx(0.9, abs=1e-10)

        # Phred 20: error = 0.01, prob = 0.99
        assert _phred_to_prob(20) == pytest.approx(0.99, abs=1e-10)

        # Phred 30: error = 0.001, prob = 0.999
        assert _phred_to_prob(30) == pytest.approx(0.999, abs=1e-10)

    def test_phred_to_error_values(self):
        """Test phred to error probability conversion."""
        # Phred 0: 10^0 = 1
        assert _phred_to_error(0) == pytest.approx(1.0, abs=1e-10)

        # Phred 10: 10^-1 = 0.1
        assert _phred_to_error(10) == pytest.approx(0.1, abs=1e-10)

        # Phred 20: 10^-2 = 0.01
        assert _phred_to_error(20) == pytest.approx(0.01, abs=1e-10)

        # Phred 30: 10^-3 = 0.001
        assert _phred_to_error(30) == pytest.approx(0.001, abs=1e-10)

    def test_phred_to_prob_high_values(self):
        """Test that high phred values are capped at max."""
        # Values >= 94 should return the same as 93
        assert _phred_to_prob(94) == _phred_to_prob(93)
        assert _phred_to_prob(100) == _phred_to_prob(93)

    def test_phred_to_error_high_values(self):
        """Test that high phred error values are capped at max."""
        assert _phred_to_error(94) == _phred_to_error(93)
        assert _phred_to_error(100) == _phred_to_error(93)


class TestConsensusRead:
    """Tests for ConsensusRead class."""

    def test_init(self):
        """Test ConsensusRead initialization."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        assert cr.contig == "chr1"
        assert cr.start_pos == 100
        assert cr.count == 10
        assert "ACGT" in cr.name
        assert "Count=10" in cr.name

    def test_add_base(self):
        """Test adding bases to consensus read."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        cr.add_base("A", "I")
        cr.add_base("C", "I")
        cr.add_base("G", "I")
        cr.add_base("T", "I")

        assert cr.seq == "ACGT"
        assert cr.qual == "IIII"
        assert len(cr.cigarstring) == 4

    def test_add_insertion(self):
        """Test adding insertions to consensus read."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        cr.add_base("A", "I")
        cr.add_insertion("CG")  # 2-base insertion
        cr.add_base("T", "I")

        assert cr.seq == "ACGT"
        assert cr.nmtag == 2  # Insertion length
        assert cr.indel_read == 1

    def test_add_deletion(self):
        """Test adding deletions to consensus read."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        cr.add_base("A", "I")
        cr.add_deletion(3)  # 3-base deletion
        cr.add_base("T", "I")

        assert cr.nmtag == 3  # Deletion length
        assert cr.indel_read == -1
        # Cigar should contain deletion markers
        assert "2" in cr.cigarstring

    def test_get_cigar(self):
        """Test CIGAR generation."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        cr.add_base("A", "I")
        cr.add_base("C", "I")
        cr.add_base("G", "I")
        cr.add_base("T", "I")

        cigar = cr.get_cigar()
        # Should be 4M (4 matches)
        assert cigar == ((0, 4),)

    def test_get_cigar_with_insertion(self):
        """Test CIGAR generation with insertion."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        cr.add_base("A", "I")
        cr.add_insertion("CG")  # 2-base insertion
        cr.add_base("T", "I")

        cigar = cr.get_cigar()
        # Should be 1M2I1M (1 match, 2 insertion, 1 match)
        assert cigar == ((0, 1), (1, 2), (0, 1))

    def test_seq_property_lazy_evaluation(self):
        """Test that seq property lazily joins parts."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGT", 10)

        cr.add_base("A", "I")
        cr.add_base("C", "I")

        # First access joins the parts
        seq1 = cr.seq
        assert seq1 == "AC"

        # Add more bases
        cr.add_base("G", "I")

        # Second access should reflect new bases
        seq2 = cr.seq
        assert seq2 == "ACG"


class TestCalcConsensus:
    """Tests for consensus calculation functions."""

    def test_calc_consensus_single_base(self):
        """Test consensus calculation with single base type."""
        cons_pos = {
            "A": [30, 30, 30],  # 3 reads with A, phred 30
            "C": [],
            "G": [],
            "T": [],
        }

        # A should have highest probability
        prob_a = calc_consensus("A", cons_pos)
        prob_c = calc_consensus("C", cons_pos)

        assert prob_a > prob_c

    def test_calc_consensus_mixed_bases(self):
        """Test consensus calculation with mixed bases."""
        cons_pos = {
            "A": [30, 30, 30],  # 3 reads with A
            "C": [30],  # 1 read with C
            "G": [],
            "T": [],
        }

        prob_a = calc_consensus("A", cons_pos)
        prob_c = calc_consensus("C", cons_pos)

        # A has more support, should have higher probability
        assert prob_a > prob_c

    def test_calc_consensus_probabilities_unanimous(self):
        """Test probability calculation with unanimous support."""
        cons_pos = {
            "A": [30, 30, 30, 30],
            "C": [],
            "G": [],
            "T": [],
        }

        base, phred = calc_consensus_probabilities(cons_pos)

        assert base == "A"
        assert phred == 60  # Max phred for unanimous

    def test_calc_consensus_probabilities_mixed(self):
        """Test probability calculation with mixed support."""
        cons_pos = {
            "A": [30, 30, 30],
            "C": [30],
            "G": [],
            "T": [],
        }

        base, phred = calc_consensus_probabilities(cons_pos)

        # A should win (more support)
        assert base == "A"
        # Phred should be positive but less than max
        assert 0 < phred <= 60


class TestGetMostCommonAllele:
    """Tests for get_most_common_allele function."""

    def test_single_base_type(self):
        """Test with only one base type."""
        cons_pos = {
            "A": [30, 30, 30],
            "C": [],
            "G": [],
            "T": [],
        }

        allele, percent = get_most_common_allele(cons_pos)

        assert allele == "A"
        assert percent == 100.0

    def test_mixed_bases(self):
        """Test with mixed bases."""
        cons_pos = {
            "A": [30, 30, 30],  # 3 reads
            "C": [30],  # 1 read
            "G": [],
            "T": [],
        }

        allele, percent = get_most_common_allele(cons_pos)

        assert allele == "A"
        assert percent == 75.0  # 3/4

    def test_with_insertion(self):
        """Test with insertion present."""
        cons_pos = {
            "A": [30, 30],
            "C": [],
            "G": [],
            "T": [],
            "I": {"CG": 3},  # Insertion of "CG" in 3 reads
        }

        allele, percent = get_most_common_allele(cons_pos)

        # Insertion should be most common
        assert allele == "ICG"
        assert percent == 60.0  # 3/5

    def test_with_deletion(self):
        """Test with deletion present."""
        cons_pos = {
            "A": [30],
            "C": [],
            "G": [],
            "T": [],
            "D": {2: 4},  # 2-base deletion in 4 reads
        }

        allele, percent = get_most_common_allele(cons_pos)

        assert allele == "D2"
        assert percent == 80.0  # 4/5


class TestGetPositionCoverage:
    """Tests for get_position_coverage function."""

    def test_bases_only(self):
        """Test coverage calculation with only bases."""
        cons_pos = {
            "A": [30, 30, 30],
            "C": [30],
            "G": [],
            "T": [],
        }

        coverage = get_position_coverage(cons_pos)
        assert coverage == 4

    def test_with_deletion(self):
        """Test coverage calculation includes deletions."""
        cons_pos = {
            "A": [30, 30],
            "C": [],
            "G": [],
            "T": [],
            "D": {2: 3},  # 3 reads with deletion
        }

        coverage = get_position_coverage(cons_pos)
        assert coverage == 5  # 2 A's + 3 deletions

    def test_with_insertion(self):
        """Test coverage calculation excludes insertions."""
        cons_pos = {
            "A": [30, 30, 30],
            "C": [],
            "G": [],
            "T": [],
            "I": {"CG": 2},  # Insertions not counted
        }

        coverage = get_position_coverage(cons_pos)
        assert coverage == 3  # Only the 3 A's

    def test_empty(self):
        """Test coverage calculation with empty position."""
        cons_pos = {
            "A": [],
            "C": [],
            "G": [],
            "T": [],
        }

        coverage = get_position_coverage(cons_pos)
        assert coverage == 0


class TestIntegration:
    """Integration tests for consensus module."""

    def test_phred_workflow(self):
        """Test full workflow: phred -> probability -> back."""
        # Start with high quality base
        original_qual = "I"  # Phred 40

        # Get phred score
        phred = get_phred(original_qual)
        assert phred == 40

        # Convert to probability
        prob = _phred_to_prob(phred)
        assert prob > 0.999  # Very high probability

        # Convert back to ASCII
        recovered_qual = get_ascii(phred)
        assert recovered_qual == original_qual

    def test_consensus_read_full_workflow(self):
        """Test building a complete consensus read."""
        cr = ConsensusRead("chr1", "region1", 100, "ACGTACGT", 5)

        # Build a simple read
        for base in "ACGTACGT":
            cr.add_base(base, "I")

        assert cr.seq == "ACGTACGT"
        assert len(cr.qual) == 8
        assert cr.count == 5

        # Get CIGAR
        cigar = cr.get_cigar()
        assert cigar == ((0, 8),)  # 8M
