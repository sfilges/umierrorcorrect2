"""Unit tests for umierrorcorrect.core.get_regions_from_bed module."""

import pytest
from umierrorcorrect.core.get_regions_from_bed import (
    get_all_annotations,
    get_first_annotation,
    get_overlap,
    merge_regions,
    read_bed,
    sort_regions,
)


class TestReadBed:
    """Tests for read_bed function."""

    def test_read_valid_bed(self, temp_bed_file):
        """Test reading a valid BED file."""
        regions = read_bed(temp_bed_file)

        assert "chr1" in regions
        assert "chr2" in regions
        assert len(regions["chr1"]) == 2
        assert len(regions["chr2"]) == 1

        # Check first region
        assert regions["chr1"][0] == (100, 200, "region1")
        assert regions["chr1"][1] == (300, 400, "region2")
        assert regions["chr2"][0] == (1000, 1100, "region3")

    def test_read_empty_bed(self, temp_output_dir):
        """Test reading an empty BED file."""
        empty_bed = temp_output_dir / "empty.bed"
        empty_bed.write_text("")
        regions = read_bed(empty_bed)

        assert regions == {}

    def test_read_bed_with_extra_columns(self, temp_output_dir):
        """Test reading BED file with extra columns (>4)."""
        bed_path = temp_output_dir / "extra_cols.bed"
        bed_content = "chr1\t100\t200\tregion1\t0\t+\n"
        bed_path.write_text(bed_content)
        regions = read_bed(bed_path)

        assert len(regions["chr1"]) == 1
        # Should only use first 4 columns
        assert regions["chr1"][0] == (100, 200, "region1")

    def test_read_bed_skips_short_lines(self, temp_output_dir):
        """Test that lines with <4 columns are skipped."""
        bed_path = temp_output_dir / "short_lines.bed"
        bed_content = """chr1\t100\t200\tregion1
chr1\t300\t400
chr2\t500\t600\tregion2
"""
        bed_path.write_text(bed_content)
        regions = read_bed(bed_path)

        # Second line should be skipped
        assert len(regions["chr1"]) == 1
        assert len(regions["chr2"]) == 1


class TestSortRegions:
    """Tests for sort_regions function."""

    def test_sort_unsorted_regions(self, sample_bed_regions_unsorted):
        """Test sorting unsorted regions."""
        sorted_regions = sort_regions(sample_bed_regions_unsorted)

        # Regions should now be sorted by start position
        starts = [r[0] for r in sorted_regions["chr1"]]
        assert starts == [100, 300, 450]

    def test_sort_already_sorted(self, sample_bed_regions):
        """Sorting already sorted regions should not change order."""
        sorted_regions = sort_regions(sample_bed_regions)

        starts = [r[0] for r in sorted_regions["chr1"]]
        assert starts == [100, 300, 450]

    def test_sort_preserves_contigs(self, sample_bed_regions):
        """Sorting should preserve all contigs."""
        sorted_regions = sort_regions(sample_bed_regions)

        assert set(sorted_regions.keys()) == set(sample_bed_regions.keys())


class TestMergeRegions:
    """Tests for merge_regions function."""

    def test_merge_non_overlapping(self, sample_bed_regions):
        """Non-overlapping regions should remain separate."""
        merged = merge_regions(sample_bed_regions, pos_threshold=10)

        # With threshold 10, regions 100-200 and 300-400 are too far apart
        # Regions 300-400 and 450-550 are also too far apart
        assert len(merged["chr1"]) == 3

    def test_merge_with_threshold(self):
        """Test merging with position threshold."""
        regions = {
            "chr1": [
                (100, 200, "region1"),
                (210, 300, "region2"),  # Within 10bp of region1's end + threshold
            ],
        }
        merged = merge_regions(regions, pos_threshold=10)

        # With threshold 10: region1 ends at 200+10=210, region2 starts at 210
        # They should merge
        assert len(merged["chr1"]) == 1
        # Names should be joined
        assert "region1" in merged["chr1"][0][2]
        assert "region2" in merged["chr1"][0][2]

    def test_merge_expands_boundaries(self):
        """Merged regions should have expanded boundaries by threshold."""
        regions = {
            "chr1": [
                (100, 200, "region1"),
            ],
        }
        merged = merge_regions(regions, pos_threshold=10)

        # Start should be 100 - 10 = 90
        # End should be 200 + 10 = 210
        assert merged["chr1"][0][0] == 90
        assert merged["chr1"][0][1] == 210

    def test_merge_overlapping_regions(self, sample_bed_regions_overlapping):
        """Overlapping regions should be merged."""
        merged = merge_regions(sample_bed_regions_overlapping, pos_threshold=0)

        # region1 (100-200) and region2 (190-300) overlap
        # They should merge into one region
        # region3 (500-600) is separate
        assert len(merged["chr1"]) == 2

    def test_merge_preserves_multiple_contigs(self, sample_bed_regions):
        """Merging should work independently for each contig."""
        merged = merge_regions(sample_bed_regions, pos_threshold=10)

        assert "chr1" in merged
        assert "chr2" in merged
        assert len(merged["chr2"]) == 1


class TestGetFirstAnnotation:
    """Tests for get_first_annotation function (deprecated)."""

    def test_position_in_region(self, sample_bed_regions):
        """Position inside a region should return its name."""
        regions = sample_bed_regions["chr1"]

        with pytest.warns(DeprecationWarning):
            assert get_first_annotation(regions, 150) == "region1"
        with pytest.warns(DeprecationWarning):
            assert get_first_annotation(regions, 350) == "region2"
        with pytest.warns(DeprecationWarning):
            assert get_first_annotation(regions, 500) == "region3"

    def test_position_at_boundaries(self, sample_bed_regions):
        """Position at region boundaries should be included."""
        regions = sample_bed_regions["chr1"]

        with pytest.warns(DeprecationWarning):
            # At start boundary
            assert get_first_annotation(regions, 100) == "region1"
        with pytest.warns(DeprecationWarning):
            # At end boundary
            assert get_first_annotation(regions, 200) == "region1"

    def test_position_outside_regions(self, sample_bed_regions):
        """Position outside all regions should return empty string."""
        regions = sample_bed_regions["chr1"]

        with pytest.warns(DeprecationWarning):
            assert get_first_annotation(regions, 50) == ""
        with pytest.warns(DeprecationWarning):
            assert get_first_annotation(regions, 250) == ""
        with pytest.warns(DeprecationWarning):
            assert get_first_annotation(regions, 1000) == ""

    def test_empty_regions(self):
        """Empty regions list should return empty string."""
        with pytest.warns(DeprecationWarning):
            assert get_first_annotation([], 100) == ""


class TestGetAllAnnotations:
    """Tests for get_all_annotations function (returns all overlapping annotations)."""

    def test_single_annotation(self, sample_bed_regions):
        """Position in one region should return that region's name."""
        regions = sample_bed_regions["chr1"]
        result = get_all_annotations(regions, 150)

        assert result == "region1"

    def test_multiple_annotations(self, sample_bed_regions_overlapping):
        """Position in overlapping regions should return all names."""
        regions = sample_bed_regions_overlapping["chr1"]

        # Position 195 is in both region1 (100-200) and region2 (190-300)
        result = get_all_annotations(regions, 195)

        assert "region1" in result
        assert "region2" in result

    def test_no_annotation(self, sample_bed_regions):
        """Position outside all regions should return empty string."""
        regions = sample_bed_regions["chr1"]
        result = get_all_annotations(regions, 250)

        assert result == ""


class TestGetOverlap:
    """Tests for get_overlap function."""

    def test_query_overlaps_region(self, sample_bed_regions):
        """Query range overlapping a region should return region name."""
        regions = sample_bed_regions["chr1"]

        # Query 150-180 overlaps region1 (100-200)
        assert get_overlap(regions, 150, 180) == "region1"

    def test_query_spans_boundary(self, sample_bed_regions):
        """Query spanning region boundary should still return overlap."""
        regions = sample_bed_regions["chr1"]

        # Query 180-220 starts in region1, ends outside
        assert get_overlap(regions, 180, 220) == "region1"

    def test_query_contains_region(self, sample_bed_regions):
        """Query containing entire region should overlap."""
        regions = sample_bed_regions["chr1"]

        # Query 50-250 contains region1 (100-200)
        assert get_overlap(regions, 50, 250) == "region1"

    def test_query_inside_region(self, sample_bed_regions):
        """Query inside region should overlap."""
        regions = sample_bed_regions["chr1"]

        # Query 120-180 is inside region1 (100-200)
        assert get_overlap(regions, 120, 180) == "region1"

    def test_query_no_overlap(self, sample_bed_regions):
        """Query not overlapping any region should return empty string."""
        regions = sample_bed_regions["chr1"]

        # Query 220-280 is between regions
        assert get_overlap(regions, 220, 280) == ""

    def test_query_adjacent_no_overlap(self, sample_bed_regions):
        """Query adjacent but not overlapping should return empty string."""
        regions = sample_bed_regions["chr1"]

        # Query 201-299 starts just after region1, ends just before region2
        assert get_overlap(regions, 201, 299) == ""

    def test_empty_regions(self):
        """Empty regions list should return empty string."""
        assert get_overlap([], 100, 200) == ""


class TestIntegration:
    """Integration tests for BED file processing pipeline."""

    def test_full_pipeline(self, temp_bed_file):
        """Test the full read -> sort -> merge pipeline."""
        # Read
        regions = read_bed(temp_bed_file)
        assert len(regions) == 2

        # Sort
        regions = sort_regions(regions)
        starts = [r[0] for r in regions["chr1"]]
        assert starts == sorted(starts)

        # Merge with small threshold
        merged = merge_regions(regions, pos_threshold=5)

        # Verify merged regions are valid
        for contig in merged:
            for start, end, name in merged[contig]:
                assert start < end
                assert len(name) > 0
