import subprocess

import pytest


@pytest.mark.integration
@pytest.mark.requires_bwa
def test_pipeline_end_to_end(temp_output_dir, test_data_dir):
    """
    Runs the full pipeline on a small subset of real data.
    Verifies that the main output files are created and not empty.
    """
    # Inputs
    r1 = test_data_dir / "integration_R1.fastq.gz"
    r2 = test_data_dir / "integration_R2.fastq.gz"
    ref = test_data_dir / "ref.fa"
    bed = test_data_dir / "test_regions.bed"

    # Ensure inputs exist
    assert r1.exists(), "R1 file missing"
    assert r2.exists(), "R2 file missing"
    assert ref.exists(), "Reference file missing"

    # Output directory is provided by fixture

    # Construct command
    # We use the CLI entry point 'umierrorcorrect2' assuming it's installed in the environment
    # or accessible via 'python -m umierrorcorrect'

    cmd = [
        "umierrorcorrect2",
        "batch",
        "-r1",
        str(r1),
        "-r2",
        str(r2),
        "-r",
        str(ref),
        "-o",
        str(temp_output_dir),
        "-rb",
        str(bed),
        "-ul",
        "19",  # UMI length
        "-sl",
        "16",  # Spacer length
        "--no-fastp",  # Skip fastp for speed/simplicity if not strictly needed, or keep if installed
        "--threads",
        "2",
    ]

    # Run pipeline
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check for success
    if result.returncode != 0:
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)

    assert result.returncode == 0, f"Pipeline failed with return code {result.returncode}"

    # Verify outputs
    # The pipeline creates a 'samples' subdirectory
    sample_dir = temp_output_dir / "samples" / "integration"
    assert sample_dir.exists(), f"Sample output directory not created at {sample_dir}"

    # Expected files
    cons_bam = sample_dir / "integration_consensus_reads.bam"
    cons_stats = sample_dir / "integration_cons.tsv"
    vcf = sample_dir / "integration.vcf"

    assert cons_bam.exists(), "Consensus BAM missing"
    assert cons_stats.exists(), "Consensus stats file missing"
    assert vcf.exists(), "VCF file missing"

    # Basic content checks
    assert cons_bam.stat().st_size > 0, "Consensus BAM is empty"
    assert cons_stats.stat().st_size > 0, "Consensus stats file is empty"
    # VCF might be empty of variants, but should have a header
    assert vcf.stat().st_size > 0, "VCF file is empty"

    print("Pipeline finished successfully and outputs verified.")
