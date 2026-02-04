#!/bin/bash
#
# Quick benchmark for rapid iteration during development
# Runs fewer iterations and tests both full pipeline and UMI consensus step
#
# Usage: ./quick_benchmark.sh [--with-fastp]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V2_DIR="$(dirname "$SCRIPT_DIR")"
V1_WRAPPER="$SCRIPT_DIR/umi_pipeline_parallel.sh"

REF="$V2_DIR/data/Homo_sapiens_PIK3CA_sequence.fa"
BED="$V2_DIR/data/test_regions.bed"
TEST_DIR="/Users/ctosimsen/Documents/GitHub/test_data"
R1="$TEST_DIR/LOD-PiK3CA-100-ng-WT-6_R1_001.fastq.gz"
R2="$TEST_DIR/LOD-PiK3CA-100-ng-WT-6_R2_001.fastq.gz"
OUT_DIR="$SCRIPT_DIR/quick_results"
THREADS=4

WITH_FASTP=false
[[ "${1:-}" == "--with-fastp" ]] && WITH_FASTP=true

echo "Quick UMI Benchmark"
echo "==================="
echo ""

# Create BWA index if needed
if [[ ! -f "${REF}.bwt" ]]; then
    echo "Creating BWA index..."
    bwa index "$REF"
fi

mkdir -p "$OUT_DIR/work"

if $WITH_FASTP; then
    #######################################
    # Full pipeline WITH fastp (1 run each)
    #######################################
    echo "Testing full pipeline WITH fastp (1 run each)..."
    echo ""

    # Time v1 with fastp wrapper
    echo "Running v1 (with fastp wrapper)..."
    rm -rf "$OUT_DIR/work/v1_fastp"
    mkdir -p "$OUT_DIR/work/v1_fastp"
    cp "$R1" "$R2" "$OUT_DIR/work/v1_fastp/"
    V1_START=$(date +%s.%N)
    bash "$V1_WRAPPER" \
        -i "$OUT_DIR/work/v1_fastp" \
        -r "$REF" \
        -b "$BED" \
        -u 19 -s 16 -t $THREADS \
        --merge --skip-fastqc --skip-multiqc 2>/dev/null
    V1_END=$(date +%s.%N)
    V1_TIME=$(echo "$V1_END - $V1_START" | bc)
    echo "  v1 (fastp): ${V1_TIME}s"

    # Time v2 with fastp
    echo "Running v2 (with fastp)..."
    rm -rf "$OUT_DIR/work/v2_fastp"
    V2_START=$(date +%s.%N)
    umierrorcorrect2 run -j 1 \
        -r1 "$R1" -r2 "$R2" -r "$REF" -rb "$BED" \
        -ul 19 -sl 16 -t $THREADS \
        -o "$OUT_DIR/work/v2_fastp" \
        --no-qc 2>/dev/null
    V2_END=$(date +%s.%N)
    V2_TIME=$(echo "$V2_END - $V2_START" | bc)
    echo "  v2 (fastp): ${V2_TIME}s"

else
    #######################################
    # Full pipeline NO fastp (1 run each)
    #######################################
    echo "Testing full pipeline (1 run each, no fastp)..."
    echo ""

    # Time v1
    echo "Running v1..."
    rm -rf "$OUT_DIR/work/v1"
    V1_START=$(date +%s.%N)
    run_umierrorcorrect.py \
        -r1 "$R1" -r2 "$R2" -r "$REF" -bed "$BED" \
        -ul 19 -sl 16 -t $THREADS -mode paired \
        -o "$OUT_DIR/work/v1" 2>/dev/null
    V1_END=$(date +%s.%N)
    V1_TIME=$(echo "$V1_END - $V1_START" | bc)
    echo "  v1: ${V1_TIME}s"

    # Time v2
    echo "Running v2..."
    rm -rf "$OUT_DIR/work/v2"
    V2_START=$(date +%s.%N)
    umierrorcorrect2 run -j 1 \
        -r1 "$R1" -r2 "$R2" -r "$REF" -rb "$BED" \
        -ul 19 -sl 16 -t $THREADS \
        -o "$OUT_DIR/work/v2" \
        --no-fastp --no-qc 2>/dev/null
    V2_END=$(date +%s.%N)
    V2_TIME=$(echo "$V2_END - $V2_START" | bc)
    echo "  v2: ${V2_TIME}s"
fi

# Calculate speedup
SPEEDUP=$(echo "scale=2; $V1_TIME / $V2_TIME" | bc)
echo ""
echo "Speedup (v1/v2): ${SPEEDUP}x"

#######################################
# UMI consensus only (using hyperfine)
#######################################

echo ""
echo "UMI Consensus benchmark (3 runs)..."
echo ""

V1_BAM=$(find "$OUT_DIR/work/v1" -name "*.sorted.bam" | head -1)
V2_BAM=$(find "$OUT_DIR/work/v2" -name "*.sorted.bam" | head -1)

if [[ -n "$V1_BAM" ]] && [[ -n "$V2_BAM" ]]; then
    hyperfine \
        --warmup 1 \
        --runs 3 \
        --prepare "rm -rf '$OUT_DIR/work/v1_cons' '$OUT_DIR/work/v2_cons'" \
        --command-name "v1-consensus" \
        "umi_error_correct.py -b '$V1_BAM' -r '$REF' -o '$OUT_DIR/work/v1_cons' 2>/dev/null" \
        --command-name "v2-consensus" \
        "umierrorcorrect2 consensus -b '$V2_BAM' -r '$REF' -o '$OUT_DIR/work/v2_cons' 2>/dev/null"
else
    echo "Could not find BAM files for consensus benchmark"
fi

echo ""
echo "Done! Results in $OUT_DIR"
