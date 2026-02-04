#!/bin/bash
#
# Benchmarking script for UMIErrorCorrect v1 vs v2
# Uses hyperfine for reliable timing measurements
#
# Usage: ./benchmark_umi.sh [options]
#   -t THREADS   Number of threads (default: 4)
#   -r RUNS      Number of benchmark runs (default: 3)
#   -o OUTDIR    Output directory (default: ./benchmark_results)
#   -s           Skip full pipeline benchmarks (only run UMI-only)
#   -h           Show help

set -euo pipefail

#######################################
# Configuration
#######################################

# Paths - adjust these as needed
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V1_WRAPPER="$SCRIPT_DIR/umi_pipeline_parallel.sh"

# Reference and test data
REF="/Users/ctosimsen/Documents/data/genomes/hg38/hg38.fa"
TEST_DIR="/Users/ctosimsen/Documents/GitHub/test_data"

# Test samples (using the smaller WT sample for faster benchmarks)
SAMPLE1_R1="$TEST_DIR/LOD-PiK3CA-100-ng-WT-6_R1_001.fastq.gz"
SAMPLE1_R2="$TEST_DIR/LOD-PiK3CA-100-ng-WT-6_R2_001.fastq.gz"
SAMPLE1_NAME="LOD-PiK3CA-100-ng-WT-6"

# Alternatively, use the larger VAF sample for more realistic benchmarks
# SAMPLE1_R1="$TEST_DIR/LOD-PiK3CA-100-ng-0-02-VAF-5_R1_001.fastq.gz"
# SAMPLE1_R2="$TEST_DIR/LOD-PiK3CA-100-ng-0-02-VAF-5_R2_001.fastq.gz"
# SAMPLE1_NAME="LOD-PiK3CA-100-ng-0-02-VAF-5"

# UMI parameters
UMI_LENGTH=19
SPACER_LENGTH=16

# Defaults
THREADS=4
RUNS=3
OUT_DIR="$SCRIPT_DIR/results"
SKIP_FULL=false
PRESERVE=true  # Preserve final results for mutation analysis

#######################################
# Parse arguments
#######################################

display_help() {
    echo "Usage: $(basename "$0") [options]"
    echo ""
    echo "Benchmark UMIErrorCorrect v1 vs v2"
    echo ""
    echo "Options:"
    echo "  -t THREADS   Number of threads (default: $THREADS)"
    echo "  -r RUNS      Number of benchmark runs (default: $RUNS)"
    echo "  -o OUTDIR    Output directory (default: $OUT_DIR)"
    echo "  -s           Skip full pipeline benchmarks (only run UMI-only)"
    echo "  -p           Don't preserve results (default: preserve with timestamp)"
    echo "  -h           Show this help"
    echo ""
    exit 0
}

while getopts "t:r:o:sph" opt; do
    case $opt in
        t) THREADS=$OPTARG ;;
        r) RUNS=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        s) SKIP_FULL=true ;;
        p) PRESERVE=false ;;
        h) display_help ;;
        *) display_help ;;
    esac
done

#######################################
# Setup
#######################################

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."

    local missing=()

    command -v hyperfine &>/dev/null || missing+=("hyperfine")
    command -v bwa &>/dev/null || missing+=("bwa")
    command -v samtools &>/dev/null || missing+=("samtools")
    command -v run_umierrorcorrect.py &>/dev/null || missing+=("umierrorcorrect v1")

    # Check v2
    command -v umierrorcorrect2 &>/dev/null || missing+=("umierrorcorrect2 v2")

    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing dependencies: ${missing[*]}"
        exit 1
    fi

    log_info "All dependencies found"
}

# Check input files
check_inputs() {
    log_info "Checking input files..."

    [[ -f "$REF" ]] || { log_error "Reference not found: $REF"; exit 1; }
    [[ -f "$SAMPLE1_R1" ]] || { log_error "R1 not found: $SAMPLE1_R1"; exit 1; }
    [[ -f "$SAMPLE1_R2" ]] || { log_error "R2 not found: $SAMPLE1_R2"; exit 1; }

    log_info "All input files found"
}

# Create BWA index if needed
create_bwa_index() {
    if [[ ! -f "${REF}.bwt" ]]; then
        log_info "Creating BWA index for reference genome..."
        bwa index "$REF"
    else
        log_info "BWA index already exists"
    fi
}

#######################################
# Benchmark functions
#######################################

# Benchmark 1: Full pipeline WITH fastp preprocessing
benchmark_full_with_fastp() {
    log_info "=== Benchmark: Full Pipeline WITH fastp ==="

    local result_json="$OUT_DIR/01_full_pipeline_with_fastp.json"
    local v1_out="$OUT_DIR/work/v1_fastp"
    local v2_out="$OUT_DIR/work/v2_fastp"

    # v1 uses umi_pipeline_parallel.sh wrapper with fastp
    # v2 uses built-in fastp preprocessing

    hyperfine \
        --warmup 1 \
        --runs "$RUNS" \
        --export-json "$result_json" \
        --export-markdown "$OUT_DIR/01_full_pipeline_with_fastp.md" \
        --prepare "rm -rf '$v1_out' '$v2_out'; mkdir -p '$v1_out'; cp '$SAMPLE1_R1' '$SAMPLE1_R2' '$v1_out/'" \
        --command-name "v1-with-fastp-qc" \
        "bash '$V1_WRAPPER' \
            -i '$v1_out' \
            -r '$REF' \
            -u $UMI_LENGTH \
            -s $SPACER_LENGTH \
            -t $THREADS \
            --merge 2>/dev/null" \
        --command-name "v2-with-fastp-qc" \
        "umierrorcorrect2 run \
            -r1 '$SAMPLE1_R1' \
            -r2 '$SAMPLE1_R2' \
            -r '$REF' \
            -ul $UMI_LENGTH \
            -sl $SPACER_LENGTH \
            -t $THREADS \
            -j 1 \
            -o '$v2_out' 2>/dev/null"

    log_info "Results saved to $result_json"
}

# Benchmark 2: Full pipeline WITHOUT fastp (raw umierrorcorrect)
benchmark_full_no_fastp() {
    log_info "=== Benchmark: Full Pipeline WITHOUT fastp ==="

    local result_json="$OUT_DIR/02_full_pipeline_no_fastp.json"
    local v1_out="$OUT_DIR/work/v1_nofastp"
    local v2_out="$OUT_DIR/work/v2_nofastp"

    hyperfine \
        --warmup 1 \
        --runs "$RUNS" \
        --export-json "$result_json" \
        --export-markdown "$OUT_DIR/02_full_pipeline_no_fastp.md" \
        --prepare "rm -rf '$v1_out' '$v2_out'" \
        --command-name "v1-no-fastp" \
        "run_umierrorcorrect.py \
            -r1 '$SAMPLE1_R1' \
            -r2 '$SAMPLE1_R2' \
            -r '$REF' \
            -ul $UMI_LENGTH \
            -sl $SPACER_LENGTH \
            -t $THREADS \
            -mode paired \
            -o '$v1_out' 2>/dev/null" \
        --command-name "v2-no-fastp" \
        "umierrorcorrect2 run \
            -r1 '$SAMPLE1_R1' \
            -r2 '$SAMPLE1_R2' \
            -r '$REF' \
            -ul $UMI_LENGTH \
            -sl $SPACER_LENGTH \
            -t $THREADS \
            -j 1 \
            -o '$v2_out' \
            --no-fastp \
            --no-qc 2>/dev/null"

    log_info "Results saved to $result_json"
}

# Prepare BAM files for UMI-only benchmarking
prepare_bam_files() {
    log_info "Preparing BAM files for UMI-only benchmark..."

    local v1_prep="$OUT_DIR/work/v1_prep"
    local v2_prep="$OUT_DIR/work/v2_prep"

    # Prepare v1 BAM (run full pipeline once)
    if [[ ! -d "$v1_prep" ]] || [[ -z "$(find "$v1_prep" -name '*.sorted.bam' 2>/dev/null)" ]]; then
        log_info "Running v1 preprocessing + mapping..."
        rm -rf "$v1_prep"
        run_umierrorcorrect.py \
            -r1 "$SAMPLE1_R1" \
            -r2 "$SAMPLE1_R2" \
            -r "$REF" \
            -ul $UMI_LENGTH \
            -sl $SPACER_LENGTH \
            -t $THREADS \
            -mode paired \
            -o "$v1_prep" 2>/dev/null
    fi

    # Prepare v2 BAM (run preprocess + mapping)
    if [[ ! -d "$v2_prep" ]] || [[ -z "$(find "$v2_prep" -name '*.sorted.bam' 2>/dev/null)" ]]; then
        log_info "Running v2 preprocessing + mapping..."
        rm -rf "$v2_prep"
        mkdir -p "$v2_prep"

        umierrorcorrect2 preprocess \
            -r1 "$SAMPLE1_R1" \
            -r2 "$SAMPLE1_R2" \
            -ul $UMI_LENGTH \
            -sl $SPACER_LENGTH \
            -o "$v2_prep" \
            -t $THREADS \
            --no-fastp 2>/dev/null

        umierrorcorrect2 mapping \
            -i "$v2_prep" \
            -r "$REF" \
            -o "$v2_prep" \
            -t $THREADS 2>/dev/null
    fi

    # Find BAM files
    V1_BAM=$(find "$v1_prep" -name "*.sorted.bam" | head -1)
    V2_BAM=$(find "$v2_prep" -name "*.sorted.bam" | head -1)

    if [[ -z "$V1_BAM" ]] || [[ -z "$V2_BAM" ]]; then
        log_error "Could not find BAM files for UMI-only benchmark"
        log_error "V1 BAM: $V1_BAM"
        log_error "V2 BAM: $V2_BAM"
        exit 1
    fi

    log_info "V1 BAM: $V1_BAM"
    log_info "V2 BAM: $V2_BAM"

    # Export for use in benchmark
    export V1_BAM V2_BAM
}

# Benchmark 3: UMI consensus generation only
benchmark_umi_only() {
    log_info "=== Benchmark: UMI Consensus Only ==="

    prepare_bam_files

    local result_json="$OUT_DIR/03_umi_consensus_only.json"
    local v1_cons="$OUT_DIR/work/v1_cons_bench"
    local v2_cons="$OUT_DIR/work/v2_cons_bench"

    hyperfine \
        --warmup 1 \
        --runs "$RUNS" \
        --export-json "$result_json" \
        --export-markdown "$OUT_DIR/03_umi_consensus_only.md" \
        --prepare "rm -rf '$v1_cons' '$v2_cons'" \
        --command-name "v1-umi-consensus" \
        "umi_error_correct.py \
            -b '$V1_BAM' \
            -r '$REF' \
            -o '$v1_cons' 2>/dev/null" \
        --command-name "v2-umi-consensus" \
        "umierrorcorrect2 consensus \
            -b '$V2_BAM' \
            -r '$REF' \
            -o '$v2_cons' 2>/dev/null"

    log_info "Results saved to $result_json"
}

# Run final production run and preserve outputs for mutation analysis
preserve_results() {
    local timestamp=$(date +%Y%m%d_%H%M%S)
    local archive_dir="$OUT_DIR/archive_${timestamp}"

    log_info "=== Running Final Production Run for Mutation Analysis ==="

    mkdir -p "$archive_dir/work"

    # Run v1 with fastp (final production run)
    log_info "Running v1 (with fastp) - production run..."
    local v1_out="$archive_dir/work/v1_fastp"
    mkdir -p "$v1_out"
    cp "$SAMPLE1_R1" "$SAMPLE1_R2" "$v1_out/"
    bash "$V1_WRAPPER" \
        -i "$v1_out" \
        -r "$REF" \
        -u $UMI_LENGTH \
        -s $SPACER_LENGTH \
        -t $THREADS \
        --merge 2>/dev/null || log_warn "v1 production run had warnings"

    # Run v2 with fastp (final production run)
    log_info "Running v2 (with fastp) - production run..."
    local v2_out="$archive_dir/work/v2_fastp"
    umierrorcorrect2 run \
        -r1 "$SAMPLE1_R1" \
        -r2 "$SAMPLE1_R2" \
        -r "$REF" \
        -ul $UMI_LENGTH \
        -sl $SPACER_LENGTH \
        -t $THREADS \
        -j 1 \
        -o "$v2_out" 2>/dev/null || log_warn "v2 production run had warnings"

    # Copy benchmark results (JSON, MD) to archive
    cp "$OUT_DIR"/*.json "$archive_dir/" 2>/dev/null || true
    cp "$OUT_DIR"/*.md "$archive_dir/" 2>/dev/null || true

    log_info "Results preserved to: $archive_dir"
    log_info "Pipeline outputs for mutation analysis:"
    ls -la "$archive_dir/work/" 2>/dev/null || true
}

# Generate summary report
generate_summary() {
    log_info "=== Generating Summary Report ==="

    local summary="$OUT_DIR/SUMMARY.md"

    cat > "$summary" << EOF
# UMIErrorCorrect Benchmark Summary

**Date:** $(date)
**Threads:** $THREADS
**Runs per benchmark:** $RUNS
**Sample:** $SAMPLE1_NAME
**UMI Length:** $UMI_LENGTH
**Spacer Length:** $SPACER_LENGTH

## System Info
- **OS:** $(uname -s) $(uname -r)
- **CPU:** $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "N/A")
- **Cores:** $(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo "N/A")

## Results

EOF

    # Append individual benchmark results
    for md in "$OUT_DIR"/*.md; do
        if [[ -f "$md" ]] && [[ "$md" != "$summary" ]]; then
            echo "### $(basename "$md" .md)" >> "$summary"
            echo "" >> "$summary"
            cat "$md" >> "$summary"
            echo "" >> "$summary"
        fi
    done

    log_info "Summary saved to $summary"
}

#######################################
# Main
#######################################

main() {
    echo ""
    log_info "UMIErrorCorrect v1 vs v2 Benchmark"
    log_info "=================================="
    echo ""

    check_dependencies
    check_inputs

    mkdir -p "$OUT_DIR/work"

    create_bwa_index

    if [[ "$SKIP_FULL" == false ]]; then
        benchmark_full_with_fastp
        benchmark_full_no_fastp
    fi

    #benchmark_umi_only

    generate_summary

    # Preserve results for later mutation analysis
    if [[ "$PRESERVE" == true ]]; then
        preserve_results
    fi

    echo ""
    log_info "All benchmarks complete!"
    log_info "Results directory: $OUT_DIR"
    echo ""
}

main "$@"
