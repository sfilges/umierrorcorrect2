#!/usr/bin/env python3
"""Benchmark script for comparing hamming distance implementations.

Run with: python scripts/benchmark_hamming.py

Requires numba for full comparison: pip install numba
"""

import random
import time
from operator import ne

# Configuration
UMI_LENGTH = 19
NUM_UMIS = 10_000
NUM_COMPARISONS = 1_000_000
NUCLEOTIDES = "ACGT"


def generate_random_umi(length: int = UMI_LENGTH) -> str:
    """Generate a random UMI sequence."""
    return "".join(random.choice(NUCLEOTIDES) for _ in range(length))


def generate_umi_pairs(num_umis: int, num_comparisons: int) -> list[tuple[str, str]]:
    """Generate random UMI pairs for benchmarking."""
    umis = [generate_random_umi() for _ in range(num_umis)]
    pairs = [(random.choice(umis), random.choice(umis)) for _ in range(num_comparisons)]
    return pairs


# Implementation 1: Original (baseline)
def hamming_original(a: str, b: str) -> int:
    """Original implementation from the codebase."""
    return sum(i != j for i, j in zip(a, b))


# Implementation 2: Optimized native with map
def hamming_native_map(a: str, b: str) -> int:
    """Optimized native Python using map."""
    return sum(map(ne, a, b))


# Implementation 3: Numba JIT (if available)
try:
    from numba import njit

    @njit(cache=True)
    def _hamming_numba_core(a: bytes, b: bytes) -> int:
        count = 0
        for i in range(len(a)):
            if a[i] != b[i]:
                count += 1
        return count

    def hamming_numba(a: str, b: str) -> int:
        return _hamming_numba_core(a.encode(), b.encode())

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    hamming_numba = None


def benchmark_function(func, pairs: list[tuple[str, str]], name: str) -> float:
    """Benchmark a hamming distance function and return elapsed time."""
    start = time.perf_counter()
    for a, b in pairs:
        func(a, b)
    elapsed = time.perf_counter() - start
    return elapsed


def warmup_numba(pairs: list[tuple[str, str]]) -> None:
    """Warm up Numba JIT compilation."""
    if HAS_NUMBA:
        # Run a few times to trigger JIT compilation
        for a, b in pairs[:100]:
            hamming_numba(a, b)


def main() -> None:
    print("=" * 60)
    print("Hamming Distance Benchmark")
    print("=" * 60)
    print(f"UMI length: {UMI_LENGTH} bp")
    print(f"Number of unique UMIs: {NUM_UMIS:,}")
    print(f"Number of comparisons: {NUM_COMPARISONS:,}")
    print()

    print("Generating test data...")
    random.seed(42)  # Reproducible results
    pairs = generate_umi_pairs(NUM_UMIS, NUM_COMPARISONS)
    print(f"Generated {len(pairs):,} UMI pairs")
    print()

    # Warm up Numba if available
    if HAS_NUMBA:
        print("Warming up Numba JIT...")
        warmup_numba(pairs)
        print()

    # Run benchmarks
    results = []

    print("Running benchmarks...")
    print("-" * 60)

    # Original
    elapsed = benchmark_function(hamming_original, pairs, "Original")
    results.append(("Original (zip + generator)", elapsed))
    print(f"Original (zip + generator):  {elapsed:.3f}s")

    # Native map
    elapsed = benchmark_function(hamming_native_map, pairs, "Native map")
    results.append(("Native (map + ne)", elapsed))
    print(f"Native (map + ne):           {elapsed:.3f}s")

    # Numba
    if HAS_NUMBA:
        elapsed = benchmark_function(hamming_numba, pairs, "Numba")
        results.append(("Numba JIT", elapsed))
        print(f"Numba JIT:                   {elapsed:.3f}s")
    else:
        print("Numba JIT:                   (not installed)")

    print("-" * 60)
    print()

    # Summary
    print("Summary")
    print("=" * 60)
    baseline = results[0][1]
    for name, elapsed in results:
        speedup = baseline / elapsed
        ops_per_sec = NUM_COMPARISONS / elapsed
        print(f"{name:28s} {elapsed:8.3f}s  {speedup:6.2f}x  {ops_per_sec:,.0f} ops/s")

    print()
    if HAS_NUMBA:
        print("Numba is installed and active.")
    else:
        print("Install numba for faster performance: pip install 'umierrorcorrect[fast]'")


if __name__ == "__main__":
    main()
