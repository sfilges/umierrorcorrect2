from __future__ import annotations

from pathlib import Path

from umierrorcorrect2.analysis.models import Mutation, MutationResult


class MutationTracker:
    """Extract mutation information from consensus TSV files."""

    def __init__(self, family_size: int = 3):
        self.family_size = family_size

    def track_mutations(
        self,
        cons_file: Path,
        mutations: list[Mutation],
        ml_plasma: float | None = None,
    ) -> list[MutationResult]:
        """
        Find mutation positions in consensus TSV file and calculate metrics for a specific family size.
        """
        results = []
        # Index mutations by (contig, position) for faster lookup
        # BED is 0-based start, .cons is 1-based pos
        mut_map = {(m.chromosome, m.position + 1): m for m in mutations}
        found_muts = {}

        with cons_file.open() as f:
            for line in f:
                parts = line.strip().split("	")
                if len(parts) < 17:
                    continue

                contig = parts[1]
                pos = int(parts[2])
                fsize = int(parts[13])

                if (contig, pos) in mut_map and fsize == self.family_size:
                    mutation = mut_map[(contig, pos)]

                    # Columns for alleles: A=5, C=6, G=7, T=8, I=9, D=10, N=11 (0-indexed)
                    # Mapping base to column index
                    base_map = {"A": 5, "C": 6, "G": 7, "T": 8, "I": 9, "D": 10}

                    alt = mutation.alternate
                    if alt not in base_map:
                        # Handle cases like complex indels if needed, but BED usually has simple ref/alt
                        mm_count = 0
                    else:
                        idx = base_map[alt]
                        mm_count = int(parts[idx])

                    total_reads = int(parts[12])  # tot coverage at this fsize
                    vaf = (mm_count / total_reads * 100.0) if total_reads > 0 else 0.0

                    mm_per_ml = (mm_count / ml_plasma) if (ml_plasma and ml_plasma > 0) else None

                    res = MutationResult(
                        mutation_name=mutation.name,
                        vaf=vaf,
                        mm_count=mm_count,
                        total_reads=total_reads,
                        mm_per_ml=mm_per_ml,
                    )
                    found_muts[mutation.name] = res

        # Ensure all requested mutations have a result (even if 0)
        for m in mutations:
            if m.name in found_muts:
                results.append(found_muts[m.name])
            else:
                results.append(
                    MutationResult(
                        mutation_name=m.name,
                        vaf=0.0,
                        mm_count=0,
                        total_reads=0,
                        mm_per_ml=0.0 if ml_plasma else None,
                    )
                )

        return results

    def create_simplified_cons(
        self,
        cons_file: Path,
        mutations: list[Mutation],
        output_file: Path,
    ) -> None:
        """Create a simplified consensus TSV file with only tracked mutations at the family size."""
        mut_map = {(m.chromosome, m.position + 1) for m in mutations}

        with cons_file.open() as f_in, output_file.open("w") as f_out:
            for line in f_in:
                parts = line.strip().split("	")
                if len(parts) < 17:
                    continue

                contig = parts[1]
                pos = int(parts[2])
                fsize = int(parts[13])

                if (contig, pos) in mut_map and fsize == self.family_size:
                    f_out.write(line)
