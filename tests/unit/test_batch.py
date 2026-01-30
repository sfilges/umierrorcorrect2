"""Unit tests for umierrorcorrect.batch module."""

import pytest
from umierrorcorrect2.batch import (
    Sample,
    discover_samples,
    parse_sample_sheet,
)
from umierrorcorrect2.core.check_args import is_tool


class TestSampleDataclass:
    """Tests for the Sample Pydantic model."""

    def test_sample_creation_paired(self, temp_output_dir):
        """Test creating a paired-end sample with file validation."""
        r1_path = temp_output_dir / "R1.fastq.gz"
        r2_path = temp_output_dir / "R2.fastq.gz"
        r1_path.touch()
        r2_path.touch()

        sample = Sample(
            name="test_sample",
            read1=r1_path,
            read2=r2_path,
        )
        assert sample.name == "test_sample"
        assert sample.read1 == r1_path
        assert sample.read2 == r2_path

    def test_sample_creation_single(self, temp_output_dir):
        """Test creating a single-end sample with file validation."""
        r1_path = temp_output_dir / "R1.fastq.gz"
        r1_path.touch()

        sample = Sample(
            name="test_sample",
            read1=r1_path,
        )
        assert sample.name == "test_sample"
        assert sample.read1 == r1_path
        assert sample.read2 is None

    def test_sample_validation_missing_read1(self, temp_output_dir):
        """Test that Sample raises error for non-existent read1."""
        with pytest.raises(ValueError, match="Read1 file not found"):
            Sample(
                name="test_sample",
                read1=temp_output_dir / "nonexistent.fastq.gz",
            )

    def test_sample_validation_missing_read2(self, temp_output_dir):
        """Test that Sample raises error for non-existent read2."""
        r1_path = temp_output_dir / "R1.fastq.gz"
        r1_path.touch()

        with pytest.raises(ValueError, match="Read2 file not found"):
            Sample(
                name="test_sample",
                read1=r1_path,
                read2=temp_output_dir / "nonexistent.fastq.gz",
            )


class TestDiscoverSamples:
    """Tests for discover_samples function."""

    def test_discover_r1_r2_pattern(self, temp_output_dir):
        """Test discovering samples with _R1/_R2 naming pattern."""
        # Create mock FASTQ files
        (temp_output_dir / "sample1_R1_001.fastq.gz").touch()
        (temp_output_dir / "sample1_R2_001.fastq.gz").touch()
        (temp_output_dir / "sample2_R1_001.fastq.gz").touch()
        (temp_output_dir / "sample2_R2_001.fastq.gz").touch()

        samples = discover_samples(temp_output_dir)

        assert len(samples) == 2
        sample_names = [s.name for s in samples]
        assert "sample1" in sample_names
        assert "sample2" in sample_names

        for sample in samples:
            assert sample.read1.exists()
            assert sample.read2 is not None
            assert sample.read2.exists()

    def test_discover_1_2_pattern(self, temp_output_dir):
        """Test discovering samples with _1/_2 naming pattern."""
        # Create mock FASTQ files
        (temp_output_dir / "sampleA_1.fastq.gz").touch()
        (temp_output_dir / "sampleA_2.fastq.gz").touch()

        samples = discover_samples(temp_output_dir)

        assert len(samples) == 1
        assert samples[0].name == "sampleA"
        assert samples[0].read1.exists()
        assert samples[0].read2 is not None
        assert samples[0].read2.exists()

    def test_discover_single_end(self, temp_output_dir):
        """Test discovering single-end samples (R1 only)."""
        # Create mock FASTQ files with no R2
        (temp_output_dir / "single_R1.fastq.gz").touch()

        samples = discover_samples(temp_output_dir)

        assert len(samples) == 1
        assert samples[0].name == "single"
        assert samples[0].read1.exists()
        assert samples[0].read2 is None

    def test_discover_empty_directory(self, temp_output_dir):
        """Test discovering samples in empty directory."""
        samples = discover_samples(temp_output_dir)
        assert len(samples) == 0

    def test_discover_sorted_by_name(self, temp_output_dir):
        """Test that discovered samples are sorted by name."""
        # Create samples out of order
        (temp_output_dir / "zebra_R1.fastq.gz").touch()
        (temp_output_dir / "alpha_R1.fastq.gz").touch()
        (temp_output_dir / "middle_R1.fastq.gz").touch()

        samples = discover_samples(temp_output_dir)

        assert len(samples) == 3
        assert samples[0].name == "alpha"
        assert samples[1].name == "middle"
        assert samples[2].name == "zebra"

    def test_discover_uncompressed_fastq(self, temp_output_dir):
        """Test discovering uncompressed FASTQ files."""
        (temp_output_dir / "uncomp_R1.fastq").touch()
        (temp_output_dir / "uncomp_R2.fastq").touch()

        samples = discover_samples(temp_output_dir)

        assert len(samples) == 1
        assert samples[0].name == "uncomp"

    def test_discover_recursive(self, temp_output_dir):
        """Test discovering samples recursively in subdirectories."""
        subdir = temp_output_dir / "subdir" / "nested"
        subdir.mkdir(parents=True)
        (subdir / "nested_R1.fastq.gz").touch()
        (subdir / "nested_R2.fastq.gz").touch()

        samples = discover_samples(temp_output_dir)

        assert len(samples) == 1
        assert samples[0].name == "nested"
        assert samples[0].read1 == subdir / "nested_R1.fastq.gz"
        assert samples[0].read2 == subdir / "nested_R2.fastq.gz"


class TestParseSampleSheet:
    """Tests for parse_sample_sheet function."""

    def test_parse_csv_sample_sheet(self, temp_output_dir):
        """Test parsing a CSV sample sheet."""
        # Create sample FASTQ files
        r1_path = temp_output_dir / "sample1_R1.fastq.gz"
        r2_path = temp_output_dir / "sample1_R2.fastq.gz"
        r1_path.touch()
        r2_path.touch()

        # Create sample sheet
        sheet_path = temp_output_dir / "samples.csv"
        sheet_path.write_text(f"sample_name,read1,read2\nsample1,{r1_path},{r2_path}\n")

        samples = parse_sample_sheet(sheet_path)

        assert len(samples) == 1
        assert samples[0].name == "sample1"
        assert samples[0].read1 == r1_path
        assert samples[0].read2 == r2_path

    def test_parse_tsv_sample_sheet(self, temp_output_dir):
        """Test parsing a TSV sample sheet."""
        # Create sample FASTQ files
        r1_path = temp_output_dir / "sample1_R1.fastq.gz"
        r1_path.touch()

        # Create sample sheet
        sheet_path = temp_output_dir / "samples.tsv"
        sheet_path.write_text(f"sample_name\tread1\tread2\nsample1\t{r1_path}\t\n")

        samples = parse_sample_sheet(sheet_path)

        assert len(samples) == 1
        assert samples[0].name == "sample1"
        assert samples[0].read1 == r1_path
        assert samples[0].read2 is None

    def test_parse_alternative_column_names(self, temp_output_dir):
        """Test parsing sample sheet with alternative column names."""
        r1_path = temp_output_dir / "s1_R1.fastq.gz"
        r2_path = temp_output_dir / "s1_R2.fastq.gz"
        r1_path.touch()
        r2_path.touch()

        # Use alternative column names (sample, r1, r2)
        sheet_path = temp_output_dir / "samples.csv"
        sheet_path.write_text(f"sample,r1,r2\ns1,{r1_path},{r2_path}\n")

        samples = parse_sample_sheet(sheet_path)

        assert len(samples) == 1
        assert samples[0].name == "s1"

    def test_parse_case_insensitive_headers(self, temp_output_dir):
        """Test that column headers are case-insensitive."""
        r1_path = temp_output_dir / "test_R1.fastq.gz"
        r1_path.touch()

        sheet_path = temp_output_dir / "samples.csv"
        sheet_path.write_text(f"Sample_Name,Read1,Read2\ntest,{r1_path},\n")

        samples = parse_sample_sheet(sheet_path)

        assert len(samples) == 1
        assert samples[0].name == "test"

    def test_parse_missing_sample_name_column(self, temp_output_dir):
        """Test error when sample_name column is missing."""
        sheet_path = temp_output_dir / "bad_samples.csv"
        sheet_path.write_text("read1,read2\npath1,path2\n")

        with pytest.raises(ValueError, match="sample_name"):
            parse_sample_sheet(sheet_path)

    def test_parse_missing_read1_column(self, temp_output_dir):
        """Test error when read1 column is missing."""
        sheet_path = temp_output_dir / "bad_samples.csv"
        sheet_path.write_text("sample_name,read2\nsample1,path2\n")

        with pytest.raises(ValueError, match="read1"):
            parse_sample_sheet(sheet_path)

    def test_parse_nonexistent_read1_file(self, temp_output_dir):
        """Test error when read1 file doesn't exist."""
        sheet_path = temp_output_dir / "samples.csv"
        sheet_path.write_text("sample_name,read1\nsample1,/nonexistent/path.fastq.gz\n")

        with pytest.raises(ValueError, match="not found"):
            parse_sample_sheet(sheet_path)

    def test_parse_multiple_samples(self, temp_output_dir):
        """Test parsing sample sheet with multiple samples."""
        # Create sample files
        for i in range(3):
            (temp_output_dir / f"sample{i}_R1.fastq.gz").touch()
            (temp_output_dir / f"sample{i}_R2.fastq.gz").touch()

        # Create sample sheet
        sheet_path = temp_output_dir / "samples.csv"
        lines = ["sample_name,read1,read2"]
        for i in range(3):
            r1 = temp_output_dir / f"sample{i}_R1.fastq.gz"
            r2 = temp_output_dir / f"sample{i}_R2.fastq.gz"
            lines.append(f"sample{i},{r1},{r2}")
        sheet_path.write_text("\n".join(lines) + "\n")

        samples = parse_sample_sheet(sheet_path)

        assert len(samples) == 3


class TestIsTool:
    """Tests for is_tool function."""

    def test_existing_tool(self):
        """Test that is_tool returns True for existing tools."""
        # 'ls' should exist on any Unix-like system running this test
        assert is_tool("ls") is True

    def test_nonexistent_tool(self):
        """Test that is_tool returns False for non-existent tools."""
        assert is_tool("nonexistent_tool_xyz123") is False
