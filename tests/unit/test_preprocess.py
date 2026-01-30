from unittest.mock import patch

import pytest
from umierrorcorrect2.models.models import FastpConfig
from umierrorcorrect2.preprocess import run_fastp


@pytest.fixture
def mock_is_tool():
    with patch("umierrorcorrect2.preprocess.is_tool", return_value=True) as m:
        yield m


@pytest.fixture
def mock_subprocess():
    with patch("umierrorcorrect2.preprocess.subprocess.run") as m:
        m.return_value.returncode = 0
        m.return_value.stderr = "fastp output"
        yield m


def test_run_fastp_merged_no_keep(mock_is_tool, mock_subprocess, tmp_path):
    """Test merged reads default behavior: NO unmerged outputs."""
    read1 = tmp_path / "R1.fastq.gz"
    read2 = tmp_path / "R2.fastq.gz"
    output_dir = tmp_path / "out"
    sample_name = "sample1"
    config = FastpConfig(merge_reads=True, keep_unmerged=False, umi_enabled=True, umi_length=12)

    result = run_fastp(read1, read2, output_dir, sample_name, config)

    assert result.merged_reads == output_dir / "sample1_umis_in_header.fastq.gz"
    assert result.filtered_read1 is None
    assert result.filtered_read2 is None

    # Verify args
    args = mock_subprocess.call_args[0][0]
    assert "--merge" in args
    assert "--merged_out" in args
    assert "--out1" not in args
    assert "--out2" not in args


def test_run_fastp_merged_keep_unmerged(mock_is_tool, mock_subprocess, tmp_path):
    """Test merged reads with keep_unmerged=True: output unmerged files."""
    read1 = tmp_path / "R1.fastq.gz"
    read2 = tmp_path / "R2.fastq.gz"
    output_dir = tmp_path / "out"
    sample_name = "sample1"
    config = FastpConfig(merge_reads=True, keep_unmerged=True, umi_enabled=True, umi_length=12)

    result = run_fastp(read1, read2, output_dir, sample_name, config)

    assert result.merged_reads == output_dir / "sample1_umis_in_header.fastq.gz"
    assert result.filtered_read1 == output_dir / "sample1_R1_unmerged_umis_in_header.fastq.gz"
    assert result.filtered_read2 == output_dir / "sample1_R2_unmerged_umis_in_header.fastq.gz"

    args = mock_subprocess.call_args[0][0]
    assert "--merge" in args
    assert "--out1" in args
    assert "--out2" in args


def test_run_fastp_pair_no_merge(mock_is_tool, mock_subprocess, tmp_path):
    """Test paired reads without merging."""
    read1 = tmp_path / "R1.fastq.gz"
    read2 = tmp_path / "R2.fastq.gz"
    output_dir = tmp_path / "out"
    sample_name = "sample1"
    config = FastpConfig(merge_reads=False, umi_enabled=True, umi_length=12)

    result = run_fastp(read1, read2, output_dir, sample_name, config)

    assert result.merged_reads is None
    assert result.filtered_read1 == output_dir / "sample1_R1_umis_in_header.fastq.gz"
    assert result.filtered_read2 == output_dir / "sample1_R2_umis_in_header.fastq.gz"

    args = mock_subprocess.call_args[0][0]
    assert "--merge" not in args
    assert "--out1" in args
    assert "--out2" in args


def test_run_fastp_single(mock_is_tool, mock_subprocess, tmp_path):
    """Test single end reads."""
    read1 = tmp_path / "R1.fastq.gz"
    read2 = None
    output_dir = tmp_path / "out"
    sample_name = "sample1"
    config = FastpConfig(merge_reads=False, umi_enabled=True, umi_length=12)

    result = run_fastp(read1, read2, output_dir, sample_name, config)

    assert result.merged_reads is None
    assert result.filtered_read1 == output_dir / "sample1_umis_in_header.fastq.gz"
    assert result.filtered_read2 is None

    args = mock_subprocess.call_args[0][0]
    assert "--in2" not in args
    assert "--out1" in args
