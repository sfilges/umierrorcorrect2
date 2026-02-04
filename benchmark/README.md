# UMIErrorCorrect Benchmarks

Benchmark scripts comparing UMIErrorCorrect v1 (original) vs v2 (modernized).

## Prerequisites

- **hyperfine**: `brew install hyperfine`
- **bwa**: For alignment
- **samtools**: For BAM processing
- **fastp**: For read filtering/merging
- **fastqc**: For quality control reports
- **multiqc**: For aggregated QC reports
- **GNU parallel**: For v1 wrapper parallel processing
- **umierrorcorrect v1**: Installed with `run_umierrorcorrect.py` in PATH
- **umierrorcorrect2 v2**: Globally installed with `umierrorcorrect2` in PATH

## Scripts

### `benchmark_umi.sh` - Full benchmark suite

Comprehensive benchmarking with multiple runs and statistical analysis.

```bash
./benchmark_umi.sh [options]

Options:
  -t THREADS   Number of threads (default: 4)
  -r RUNS      Number of benchmark runs (default: 3)
  -o OUTDIR    Output directory (default: ./results)
  -s           Skip full pipeline benchmarks (only run UMI-only)
  -p           Don't preserve results (default: preserve with timestamp)
  -h           Show help
```

**Benchmarks included:**
1. Full pipeline WITH fastp + QC (v1 wrapper vs v2)
2. Full pipeline WITHOUT fastp/QC (direct v1 vs v2 comparison)
3. UMI consensus generation only

### `umi_pipeline_parallel.sh` - v1 wrapper with fastp

Wrapper script around v1 that adds fastp preprocessing, parallel processing,
and QC reports. Used for fair comparison with v2's built-in fastp support.

```bash
./umi_pipeline_parallel.sh [options]

Key options:
  -i, --input-dir      Input directory with FASTQ files
  -r, --reference      Reference genome
  -b, --bed            BED file for regions
  -u, --umi_length     UMI length (default: 19)
  -s, --spacer_length  Spacer length (default: 16)
  -t, --threads        Number of parallel jobs (default: nproc/2)
  -f, --no_filtering   Skip fastp filtering
  --merge              Merge overlapping reads
  --skip-fastqc        Skip FastQC
  --skip-multiqc       Skip MultiQC
```

### `quick_benchmark.sh` - Rapid iteration

Quick single-run comparison for development.

```bash
# Without fastp (direct comparison)
./quick_benchmark.sh

# With fastp preprocessing
./quick_benchmark.sh --with-fastp
```

### `analyze_results.py` - Result analysis

Parse hyperfine JSON output and generate comparison reports.

```bash
# Analyze all results in directory
python analyze_results.py benchmark_results/

# Analyze single file
python analyze_results.py benchmark_results/03_umi_consensus_only.json
```

## Test Data

The benchmarks use test data from `/Users/ctosimsen/Documents/GitHub/test_data/`:
- `LOD-PiK3CA-100-ng-WT-6_R1_001.fastq.gz` (wild-type, smaller)
- `LOD-PiK3CA-100-ng-0-02-VAF-5_R1_001.fastq.gz` (0.02% VAF, larger)

## Output

Results are saved to `results/` (gitignored):
```
results/
├── work/                      # Temporary work directory (cleared between runs)
├── archive_YYYYMMDD_HHMMSS/   # Preserved results for mutation analysis
│   ├── work/
│   │   ├── v1_fastp/          # v1 output: BAMs, consensus, variants
│   │   ├── v2_fastp/          # v2 output: BAMs, consensus, variants
│   │   ├── v1_nofastp/
│   │   └── v2_nofastp/
│   ├── *.json                 # Hyperfine timing results
│   └── *.md                   # Markdown reports
├── *.json                     # Latest timing results
├── *.md                       # Latest markdown tables
└── SUMMARY.md                 # Combined report
```

By default, after benchmarking completes, a **final production run** is executed
to generate pipeline outputs for mutation analysis. This is separate from the
timed runs to ensure clean, complete output files.

Use `-p` to skip the production run if you only want timing data.

## Benchmark Methodology

### Fair Comparison

Both versions are tested with:
- Same UMI length (19 bp)
- Same spacer length (16 bp)
- Same thread count
- Same reference genome (hg38)
- Same input FASTQ files
- Single job mode (`-j 1` for v2)

### What's Being Measured

| Benchmark | v1 | v2 |
|-----------|----|----|
| Full pipeline (fastp + QC) | `umi_pipeline_parallel.sh --merge` | `umierrorcorrect2 run -j 1` |
| Full pipeline (no fastp/QC) | `run_umierrorcorrect.py` | `umierrorcorrect2 run -j 1 --no-fastp --no-qc` |
| UMI consensus only | `umi_error_correct.py` | `umierrorcorrect2 consensus` |

### Warmup

All benchmarks include 1 warmup run to avoid cold-cache effects.

## Interpreting Results

Hyperfine reports:
- **Mean**: Average execution time
- **Stddev**: Standard deviation (consistency)
- **Min/Max**: Range of execution times
- **Speedup**: Relative performance comparison

Lower times and smaller stddev indicate better performance.
