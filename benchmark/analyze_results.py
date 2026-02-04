#!/usr/bin/env python3
"""
Analyze hyperfine benchmark results and generate comparison reports.

Usage:
    python analyze_results.py benchmark_results/
    python analyze_results.py benchmark_results/*.json
"""

import json
import sys
from pathlib import Path
from typing import Any


def load_hyperfine_json(path: Path) -> dict[str, Any]:
    """Load hyperfine JSON result file."""
    with Path.open(path) as f:
        return json.load(f)


def format_time(seconds: float) -> str:
    """Format time in human-readable format."""
    if seconds < 1:
        return f"{seconds * 1000:.1f} ms"
    elif seconds < 60:
        return f"{seconds:.2f} s"
    else:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"


def analyze_benchmark(data: dict[str, Any]) -> None:
    """Analyze and print benchmark results."""
    results = data.get("results", [])

    if not results:
        print("No results found")
        return

    print(f"\n{'Command':<25} {'Mean':>12} {'Stddev':>10} {'Min':>12} {'Max':>12}")
    print("-" * 75)

    times = []
    for r in results:
        name = r["command"][:24]
        mean = r["mean"]
        stddev = r["stddev"]
        min_t = r["min"]
        max_t = r["max"]
        times.append((r["command"], mean))

        print(
            f"{name:<25} {format_time(mean):>12} {format_time(stddev):>10} "
            f"{format_time(min_t):>12} {format_time(max_t):>12}"
        )

    # Calculate speedup if we have exactly 2 commands
    if len(times) == 2:
        print()
        name1, time1 = times[0]
        name2, time2 = times[1]

        if time1 > time2:
            speedup = time1 / time2
            faster = name2
        else:
            speedup = time2 / time1
            faster = name1

        print(f"Speedup: {faster} is {speedup:.2f}x faster")


def analyze_all_benchmarks(result_dir: Path) -> None:
    """Analyze all benchmark JSON files in a directory."""
    json_files = sorted(result_dir.glob("*.json"))

    if not json_files:
        print(f"No JSON files found in {result_dir}")
        return

    print("=" * 75)
    print("UMIErrorCorrect Benchmark Analysis")
    print("=" * 75)

    for json_file in json_files:
        print(f"\n## {json_file.stem}")
        data = load_hyperfine_json(json_file)
        analyze_benchmark(data)

    print("\n" + "=" * 75)


def generate_csv(result_dir: Path) -> None:
    """Generate CSV summary of all benchmarks."""
    csv_path = result_dir / "summary.csv"
    json_files = sorted(result_dir.glob("*.json"))

    with Path.open(csv_path, "w") as f:
        f.write("benchmark,command,mean_s,stddev_s,min_s,max_s,runs\n")

        for json_file in json_files:
            data = load_hyperfine_json(json_file)
            benchmark_name = json_file.stem

            for r in data.get("results", []):
                f.write(
                    f"{benchmark_name},{r['command']},{r['mean']:.4f},"
                    f"{r['stddev']:.4f},{r['min']:.4f},{r['max']:.4f},"
                    f"{len(r.get('times', []))}\n"
                )

    print(f"\nCSV summary saved to {csv_path}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = Path(sys.argv[1])

    if path.is_dir():
        analyze_all_benchmarks(path)
        generate_csv(path)
    elif path.suffix == ".json":
        # Analyze single JSON file
        data = load_hyperfine_json(path)
        print(f"## {path.stem}")
        analyze_benchmark(data)
    else:
        print(f"Unknown input: {path}")
        sys.exit(1)


if __name__ == "__main__":
    main()
