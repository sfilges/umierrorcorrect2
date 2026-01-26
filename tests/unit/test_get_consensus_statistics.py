"""Unit tests for umierrorcorrect.get_consensus_statistics module."""

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES
from umierrorcorrect.get_consensus_statistics import (
    RegionConsensusStats,
    RegionStats,
    calculate_target_coverage,
    get_overall_statistics,
    parse_consensus_read_name,
)


class TestRegionConsensusStats:
    """Tests for RegionConsensusStats class."""

    def test_init_basic(self):
        """Test basic initialization."""
        fsizes = [1, 2, 3, 5, 10]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 5, fsizes)

        assert stat.regionid == "1"
        assert stat.pos == "chr1:100-200"
        assert stat.name == "gene1"
        assert stat.singletons == 5
        assert stat.family_sizes == []
        assert stat.fsizes == fsizes

    def test_init_singleton_counts(self):
        """Test that singletons are counted correctly at initialization."""
        fsizes = [1, 2, 3, 5, 10]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 5, fsizes)

        # Singletons should contribute to threshold 0 and 1
        assert stat.total_reads[0] == 5
        assert stat.umis[0] == 5
        assert stat.total_reads[1] == 5
        assert stat.umis[1] == 5

        # Higher thresholds should start at 0
        assert stat.total_reads[2] == 0
        assert stat.umis[2] == 0

    def test_init_no_singletons(self):
        """Test initialization with no singletons."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        assert stat.singletons == 0
        assert stat.total_reads[0] == 0
        assert stat.umis[0] == 0

    def test_add_family_sizes_basic(self):
        """Test adding family sizes updates statistics correctly."""
        fsizes = [1, 2, 3, 5]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        # Add families of sizes [5, 3, 2, 2, 1]
        sizes = [5, 3, 2, 2, 1]
        stat.add_family_sizes(sizes, fsizes)

        # Total reads at threshold 0: 5+3+2+2+1 = 13
        assert stat.total_reads[0] == 13
        # UMIs at threshold 0: Should equal total reads (13) because fsize=0 implies "raw reads"
        assert stat.umis[0] == 13

        # At threshold 1: all families pass
        assert stat.total_reads[1] == 13
        assert stat.umis[1] == 5

        # At threshold 2: families [5,3,2,2] pass = 12 reads, 4 families
        assert stat.total_reads[2] == 12
        assert stat.umis[2] == 4

        # At threshold 3: families [5,3] pass = 8 reads, 2 families
        assert stat.total_reads[3] == 8
        assert stat.umis[3] == 2

        # At threshold 5: families [5] pass = 5 reads, 1 family
        assert stat.total_reads[5] == 5
        assert stat.umis[5] == 1

    def test_add_family_sizes_extends_list(self):
        """Test that family sizes are added to the list."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        stat.add_family_sizes([5, 3], fsizes)
        assert stat.family_sizes == [5, 3]

        stat.add_family_sizes([2, 1], fsizes)
        assert stat.family_sizes == [5, 3, 2, 1]

    def test_add_family_sizes_multiple_calls(self):
        """Test that multiple add_family_sizes calls accumulate correctly."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        stat.add_family_sizes([5, 3], fsizes)
        stat.add_family_sizes([2, 1], fsizes)

        # Total: 5+3+2+1 = 11 reads
        # UMIs at threshold 0 should equal total reads (11)
        assert stat.total_reads[0] == 11
        assert stat.umis[0] == 11

    def test_write_stats_format(self):
        """Test that write_stats produces correct format."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 2, fsizes)
        stat.add_family_sizes([5, 3], fsizes)

        output = stat.write_stats()
        lines = output.split("\n")

        # Should have one line per threshold: 0, 1, 2, 3
        assert len(lines) == 4

        # Check first line (threshold 0)
        parts = lines[0].split("\t")
        assert parts[0] == "1"  # regionid
        assert parts[1] == "chr1:100-200"  # pos
        assert parts[2] == "gene1"  # name
        assert parts[3] == "0"  # threshold


class TestCalculateTargetCoverage:
    """Tests for calculate_target_coverage function."""

    def test_does_not_mutate_input(self):
        """Test that the function does not mutate the input fsizes list."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]
        original_fsizes = fsizes.copy()

        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 5, fsizes)
        stat.add_family_sizes([10, 5, 3], fsizes)

        calculate_target_coverage([stat])

        # fsizes should be unchanged
        assert fsizes == original_fsizes

    def test_named_vs_unnamed_regions(self):
        """Test that named and unnamed regions are tracked separately."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]

        # Named region (on-target)
        named = RegionConsensusStats("1", "chr1:100-200", "gene1", 0, fsizes)
        named.add_family_sizes([10], fsizes)

        # Unnamed region (off-target)
        unnamed = RegionConsensusStats("2", "chr1:300-400", "", 0, fsizes)
        unnamed.add_family_sizes([10], fsizes)

        output = calculate_target_coverage([named, unnamed])
        lines = output.split("\n")

        # First line is threshold 0
        # Format: family_size, on_target, off_target, total, on_target_frac, off_target_frac
        parts = lines[0].split("\t")
        on_target = int(parts[1])
        off_target = int(parts[2])
        total = int(parts[3])

        # At threshold 0, on_target should be reads in named regions (10)
        # off_target should be reads in unnamed regions (10)
        # and total should be all reads (20)
        assert on_target == 10
        assert off_target == 10
        assert total == 20

    def test_empty_stats(self):
        """Test with empty statistics list."""
        output = calculate_target_coverage([])
        lines = output.split("\n")

        # Should still produce output for all thresholds
        assert len(lines) == len(DEFAULT_FAMILY_SIZES)


class TestGetOverallStatistics:
    """Tests for get_overall_statistics function."""

    def test_aggregates_across_regions(self):
        """Test that statistics are correctly aggregated across regions."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]

        stat1 = RegionConsensusStats("1", "chr1:100-200", "gene1", 2, fsizes)
        stat1.add_family_sizes([5, 3], fsizes)

        stat2 = RegionConsensusStats("2", "chr1:300-400", "gene2", 3, fsizes)
        stat2.add_family_sizes([4, 2], fsizes)

        overall = get_overall_statistics([stat1, stat2])

        # Total reads at threshold 0:
        # stat1: 5+3 + 2 singletons = 10
        # stat2: 4+2 + 3 singletons = 9
        # overall: 19
        assert overall.total_reads[0] == 19

        # UMIs at threshold 0:
        # stat1: 5+3 + 2 singletons = 10
        # stat2: 4+2 + 3 singletons = 9
        # overall: 19
        assert overall.umis[0] == 19

    def test_overall_regionid(self):
        """Test that overall statistics have correct regionid."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 0, fsizes)

        overall = get_overall_statistics([stat])

        assert overall.regionid == "All"
        assert overall.pos == "all_regions"
        assert overall.name == ""

    def test_empty_regions(self):
        """Test with empty region list."""
        overall = get_overall_statistics([])

        assert overall.total_reads[0] == 0
        assert overall.umis[0] == 0


class TestRegionStats:
    """Tests for the RegionStats dataclass."""

    def test_basic_properties(self):
        """Test basic RegionStats properties."""
        stats = RegionStats(
            chrom="chr1", start=100, end=200, name="gene1", consensus_counts=[5, 3, 10], singleton_count=2
        )

        assert stats.chrom == "chr1"
        assert stats.start == 100
        assert stats.end == 200
        assert stats.name == "gene1"
        assert stats.consensus_counts == [5, 3, 10]
        assert stats.singleton_count == 2

    def test_position_property(self):
        """Test the position property."""
        stats = RegionStats(chrom="chr1", start=100, end=200)
        assert stats.position == "chr1:100-200"

    def test_total_consensus(self):
        """Test total_consensus property."""
        stats = RegionStats(chrom="chr1", start=100, end=200, consensus_counts=[5, 3, 10])
        assert stats.total_consensus == 3

    def test_total_raw_reads(self):
        """Test total_raw_reads property."""
        stats = RegionStats(chrom="chr1", start=100, end=200, consensus_counts=[5, 3, 10], singleton_count=2)
        # 5 + 3 + 10 + 2 = 20
        assert stats.total_raw_reads == 20

    def test_umis_at_threshold(self):
        """Test umis_at_threshold method."""
        stats = RegionStats(chrom="chr1", start=100, end=200, consensus_counts=[5, 3, 10, 2, 1], singleton_count=2)

        # Threshold 1: all consensus (5) + singletons (2) = 7
        assert stats.umis_at_threshold(1) == 7

        # Threshold 2: [5, 3, 10, 2] = 4
        assert stats.umis_at_threshold(2) == 4

        # Threshold 5: [5, 10] = 2
        assert stats.umis_at_threshold(5) == 2

        # Threshold 10: [10] = 1
        assert stats.umis_at_threshold(10) == 1

    def test_reads_at_threshold(self):
        """Test reads_at_threshold method."""
        stats = RegionStats(chrom="chr1", start=100, end=200, consensus_counts=[5, 3, 10], singleton_count=2)

        # Threshold 1: 5 + 3 + 10 + 2 = 20
        assert stats.reads_at_threshold(1) == 20

        # Threshold 3: 5 + 3 + 10 = 18
        assert stats.reads_at_threshold(3) == 18

        # Threshold 5: 5 + 10 = 15
        assert stats.reads_at_threshold(5) == 15


class TestParseConsensusReadName:
    """Tests for parse_consensus_read_name function."""

    def test_consensus_read(self):
        """Test parsing a consensus read name."""
        read_type, count = parse_consensus_read_name("Consensus_read_0_TTGTAAAGCATGAAATAGC_Count=144")
        assert read_type == "Consensus"
        assert count == 144

    def test_consensus_read_with_subcluster(self):
        """Test parsing a consensus read name with subcluster tag."""
        read_type, count = parse_consensus_read_name("Consensus_read_0_TTGTAAAGCATGAAATAGC_a_Count=144")
        assert read_type == "Consensus"
        assert count == 144

    def test_singleton_read(self):
        """Test parsing a singleton read name."""
        read_type, count = parse_consensus_read_name("Singleton_read_0_AAGAAAAGCAAGAAATGGT_Count=1")
        assert read_type == "Singleton"
        assert count == 1


class TestRegionConsensusStatsFromRegionStats:
    """Tests for RegionConsensusStats.from_region_stats factory method."""

    def test_from_region_stats_basic(self):
        """Test creating RegionConsensusStats from RegionStats."""
        region = RegionStats(chrom="chr1", start=100, end=200, name="gene1", consensus_counts=[5, 3], singleton_count=2)
        fsizes = [1, 2, 3, 5]

        stat = RegionConsensusStats.from_region_stats(region, fsizes)

        assert stat.regionid == "chr1:100"
        assert stat.pos == "chr1:100-200"
        assert stat.name == "gene1"
        assert stat.singletons == 2
        assert stat.family_sizes == [5, 3]

    def test_from_region_stats_counts(self):
        """Test that counts are calculated correctly."""
        region = RegionStats(chrom="chr1", start=100, end=200, consensus_counts=[5, 3, 2], singleton_count=2)
        fsizes = [1, 2, 3, 5]

        stat = RegionConsensusStats.from_region_stats(region, fsizes)

        # Total reads at threshold 0: 5+3+2 + 2 singletons = 12
        assert stat.total_reads[0] == 12

        # At threshold 2: families [5,3,2] pass = 10 reads
        assert stat.total_reads[2] == 10
        assert stat.umis[2] == 3

        # At threshold 3: families [5,3] pass = 8 reads
        assert stat.total_reads[3] == 8
        assert stat.umis[3] == 2
