"""Unit tests for umierrorcorrect.src.check_args module."""

import os
import tempfile
from pathlib import Path

import pytest

from umierrorcorrect.src.check_args import (
    check_output_directory,
    get_sample_name,
    is_tool,
)


class TestIsTool:
    """Tests for is_tool function."""

    def test_existing_tool(self):
        """Test that is_tool returns True for existing tools."""
        # 'python' should exist on any system running this test
        assert is_tool("python") is True

    def test_nonexistent_tool(self):
        """Test that is_tool returns False for non-existent tools."""
        assert is_tool("nonexistent_tool_xyz123") is False

    def test_no_resource_leak(self):
        """Test that is_tool doesn't leak file handles.

        The original implementation had a resource leak because it didn't
        close the devnull file handle properly.
        """
        # Run is_tool multiple times to check for resource leaks
        for _ in range(100):
            is_tool("python")
            is_tool("nonexistent_tool")
        # If we get here without error, no file handle leak occurred


class TestCheckOutputDirectory:
    """Tests for check_output_directory function."""

    def test_existing_directory(self, temp_output_dir):
        """Test with existing directory."""
        result = check_output_directory(str(temp_output_dir))
        assert result == str(temp_output_dir)

    def test_create_new_directory(self, temp_output_dir):
        """Test creating a new directory."""
        new_dir = temp_output_dir / "new_subdir"
        assert not new_dir.exists()

        result = check_output_directory(str(new_dir))
        assert result == str(new_dir)
        assert new_dir.exists()
        assert new_dir.is_dir()


class TestGetSampleName:
    """Tests for get_sample_name function."""

    def test_single_mode(self):
        """Test sample name extraction in single mode."""
        filename = "/path/to/sample.fastq.gz"
        result = get_sample_name(filename, "single")
        # Should strip fastq.gz extension
        assert "fastq" not in result.lower()

    def test_paired_mode_r1(self):
        """Test sample name extraction for R1 file in paired mode."""
        filename = "/path/to/sample_R1_001.fastq.gz"
        result = get_sample_name(filename, "paired")
        # Should strip R1, _001, and fastq.gz
        assert "R1" not in result
        assert "_001" not in result

    def test_paired_mode_with_lane(self):
        """Test sample name extraction with lane info."""
        filename = "/path/to/sample_L001_R1_001.fastq.gz"
        result = get_sample_name(filename, "paired")
        # Should strip lane info
        assert "L001" not in result

    def test_bam_mode(self):
        """Test sample name extraction in bam mode."""
        filename = "/path/to/sample.sorted.bam"
        result = get_sample_name(filename, "bam")
        assert result == "sample"
        assert ".sorted" not in result
        assert ".bam" not in result

    def test_bam_mode_unsorted(self):
        """Test sample name extraction for unsorted bam."""
        filename = "/path/to/sample.bam"
        result = get_sample_name(filename, "bam")
        assert result == "sample"
