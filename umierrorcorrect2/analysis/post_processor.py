from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from umierrorcorrect2.analysis.models import load_mutations
from umierrorcorrect2.core.constants import CONSENSUS_TSV_SUFFIX
from umierrorcorrect2.core.get_regions_from_bed import get_all_annotations, read_bed
from umierrorcorrect2.core.logging_config import get_logger

logger = get_logger("post_processor")

_SUMMARY_STATS_SUFFIX = "_summary_statistics.txt"


class PostProcessor:
    """Aggregate and summarize UMI error correction results across multiple samples."""

    def __init__(self, results_dir: Path, family_size: int = 3):
        self.results_dir = results_dir
        self.family_size = family_size

    # -------------------------------------------------------------------------
    # Public methods
    # -------------------------------------------------------------------------

    def aggregate_cons(
        self,
        sample_dirs: list[Path] | None = None,
        region_bed: Path | None = None,
        mutation_bed: Path | None = None,
    ) -> pd.DataFrame:
        """Load and aggregate consensus TSV files from multiple sample directories.

        Equivalent to R's import_cons_data() in simsenR. Returns a combined DataFrame
        filtered to self.family_size. Optionally annotates with region BED intervals
        and/or mutation positions.
        """
        dirs = sample_dirs if sample_dirs is not None else self._discover_sample_dirs()
        if not dirs:
            logger.warning(f"No sample directories found under {self.results_dir}")
            return pd.DataFrame()

        frames = []
        for sample_dir in dirs:
            cons_files = list(sample_dir.glob(f"*{CONSENSUS_TSV_SUFFIX}"))
            if not cons_files:
                logger.warning(f"No {CONSENSUS_TSV_SUFFIX} found in {sample_dir}, skipping")
                continue
            try:
                df = pd.read_csv(cons_files[0], sep="\t")
                df = df[df["Consensus group size"] == self.family_size]
                frames.append(df)
            except Exception as e:
                logger.warning(f"Failed to read {cons_files[0]}: {e}")

        if not frames:
            logger.warning("No consensus data loaded from any sample directory")
            return pd.DataFrame()

        combined = pd.concat(frames, ignore_index=True)

        if region_bed is not None:
            combined = self._annotate_regions(combined, region_bed)

        if mutation_bed is not None:
            combined = self._join_mutations(combined, mutation_bed)

        return combined

    def compute_mutation_metrics(
        self,
        cons_df: pd.DataFrame,
        ml_plasma_map: dict[str, float] | None = None,
    ) -> pd.DataFrame:
        """Compute per-mutation metrics (VAF, ctDNA ppm, mm_per_ml) across samples.

        cons_df must have is_mutation, expected_alt, and alt_matches columns,
        as produced by aggregate_cons() with mutation_bed.
        """
        required = {"is_mutation", "expected_alt", "alt_matches"}
        if missing := required - set(cons_df.columns):
            raise ValueError(f"cons_df missing columns: {missing}. Call aggregate_cons(mutation_bed=...) first.")

        df = cons_df[cons_df["is_mutation"]].copy()
        if df.empty:
            return pd.DataFrame()

        # Zero out count/frequency where alt doesn't match expected
        not_detected = ~df["alt_matches"]
        df.loc[not_detected, "Max Non-ref Allele Count"] = 0
        df.loc[not_detected, "Max Non-ref Allele Frequency"] = 0.0

        coverage = df["Coverage"].replace(0, float("nan"))
        df["VAF (%)"] = (df["Max Non-ref Allele Count"] / coverage * 100.0).fillna(0.0)
        df["ctDNA ppm"] = (df["Max Non-ref Allele Count"] / coverage * 1e6).fillna(0.0)

        if ml_plasma_map:
            df["mm_per_ml"] = df.apply(
                lambda row: (
                    row["Max Non-ref Allele Count"] / ml_plasma_map[row["Sample Name"]]
                    if row["Sample Name"] in ml_plasma_map and ml_plasma_map[row["Sample Name"]] > 0
                    else None
                ),
                axis=1,
            )

        return df.reset_index(drop=True)

    def compute_on_target_fractions(
        self,
        sample_dirs: list[Path] | None = None,
        region_bed: Path | None = None,
        mutation_bed: Path | None = None,
    ) -> pd.DataFrame:
        """Compute per-sample on-target read fractions.

        Total reads denominator: fastp.json summary.after_filtering.total_reads.
        On-target numerator at Consensus group size=0:
          - region_bed provided: named-region reads from _summary_statistics.txt
          - mutation_bed provided (no region_bed): Coverage at mutation positions from _cons.tsv
          - neither: on-target not measurable
        """
        dirs = sample_dirs if sample_dirs is not None else self._discover_sample_dirs()

        # Pre-load mutation positions if needed
        mutation_positions: set[tuple[str, int]] | None = None
        if mutation_bed is not None and region_bed is None:
            mutations = load_mutations(mutation_bed)
            # BED 0-based start → 1-based cons Position
            mutation_positions = {(m.chromosome, m.position + 1) for m in mutations}

        rows = [self._compute_sample_on_target(sample_dir, region_bed, mutation_positions) for sample_dir in dirs]
        return pd.DataFrame(rows)

    # -------------------------------------------------------------------------
    # Private helpers
    # -------------------------------------------------------------------------

    def _compute_sample_on_target(
        self,
        sample_dir: Path,
        region_bed: Path | None,
        mutation_positions: set[tuple[str, int]] | None,
    ) -> dict:
        sample_name = sample_dir.name
        row: dict = {"sample_name": sample_name}

        # Total input reads from fastp.json
        fastp_path = self._find_fastp_json(sample_dir)
        if fastp_path:
            fastp_data = self._parse_fastp_json(fastp_path)
        else:
            fastp_data = None
            logger.warning(f"No fastp JSON found for {sample_name}; on-target fraction not computed")

        row["total_input_reads"] = fastp_data["total_reads_after"] if fastp_data else None
        row["q20_rate"] = fastp_data.get("q20_rate") if fastp_data else None
        row["duplication_rate"] = fastp_data.get("duplication_rate") if fastp_data else None

        # On-target reads
        if region_bed is not None:
            on_target, method = self._on_target_from_summary_stats(sample_dir)
        elif mutation_positions is not None:
            on_target = self._on_target_from_cons(sample_dir, mutation_positions)
            method = "mutation_bed"
        else:
            on_target = None
            method = "unmeasured"

        row["on_target_reads"] = on_target
        row["on_target_method"] = method

        total = row["total_input_reads"]
        row["on_target_fraction"] = (on_target / total) if (total and on_target is not None) else None

        return row

    def _discover_sample_dirs(self) -> list[Path]:
        """Return subdirectories of results_dir that contain a *_cons.tsv file."""
        dirs = [
            d for d in sorted(self.results_dir.iterdir()) if d.is_dir() and list(d.glob(f"*{CONSENSUS_TSV_SUFFIX}"))
        ]
        if not dirs:
            logger.warning(f"No sample subdirectories with {CONSENSUS_TSV_SUFFIX} found in {self.results_dir}")
        return dirs

    def _find_fastp_json(self, sample_dir: Path) -> Path | None:
        """Find a fastp JSON file in sample_dir."""
        for pattern in ("*.fastp.json", "*fastp*.json"):
            matches = list(sample_dir.glob(pattern))
            if matches:
                return matches[0]
        return None

    def _parse_fastp_json(self, path: Path) -> dict | None:
        """Extract total reads, Q20 rate, and duplication rate from fastp JSON."""
        try:
            with path.open() as f:
                data = json.load(f)
            after = data.get("summary", {}).get("after_filtering", {})
            dup = data.get("duplication", {})
            return {
                "total_reads_after": after.get("total_reads"),
                "q20_rate": after.get("q20_rate"),
                "duplication_rate": dup.get("rate"),
            }
        except Exception as e:
            logger.warning(f"Failed to parse fastp JSON {path}: {e}")
            return None

    def _parse_summary_stats(self, path: Path) -> pd.DataFrame:
        """Parse a _summary_statistics.txt file (7-column TSV, no header)."""
        cols = ["regionid", "pos", "name", "family_size", "fraction", "total_reads", "umis"]
        df = pd.read_csv(path, sep="\t", header=None, names=cols)
        # Drop header row if file has one (identified by non-numeric family_size)
        df = df[pd.to_numeric(df["family_size"], errors="coerce").notna()].copy()
        df["family_size"] = df["family_size"].astype(int)
        df["total_reads"] = df["total_reads"].astype(int)
        return df

    def _annotate_regions(self, df: pd.DataFrame, region_bed: Path) -> pd.DataFrame:
        """Update 'Name' column for rows where it is empty/NaN using BED intervals."""
        regions = read_bed(str(region_bed))
        df = df.copy()
        needs_annotation = df["Name"].isna() | (df["Name"] == "")
        if not needs_annotation.any():
            return df

        def annotate(row):
            contig = row["Contig"]
            # cons.tsv Position is 1-based; convert to 0-based for BED lookup
            pos = int(row["Position"]) - 1
            return get_all_annotations(regions.get(contig, []), pos)

        df.loc[needs_annotation, "Name"] = df.loc[needs_annotation].apply(annotate, axis=1)
        return df

    def _join_mutations(self, df: pd.DataFrame, mutation_bed: Path) -> pd.DataFrame:
        """Add mutation annotation columns by left-joining on (Contig, Position)."""
        mutations = load_mutations(mutation_bed)
        mut_df = pd.DataFrame(
            [
                {
                    "Contig": m.chromosome,
                    "Position": m.position + 1,  # BED 0-based → 1-based cons position
                    "mutation_name": m.name,
                    "expected_alt": m.alternate,
                }
                for m in mutations
            ]
        )
        df = df.merge(mut_df, on=["Contig", "Position"], how="left")
        df["is_mutation"] = df["mutation_name"].notna()
        df["alt_matches"] = df["is_mutation"] & (df["Max Non-ref Allele"] == df["expected_alt"])
        return df

    def _on_target_from_summary_stats(self, sample_dir: Path) -> tuple[int | None, str]:
        """Sum on-target reads from _summary_statistics.txt at family_size=0."""
        stat_files = list(sample_dir.glob(f"*{_SUMMARY_STATS_SUFFIX}"))
        if not stat_files:
            logger.warning(f"No {_SUMMARY_STATS_SUFFIX} found in {sample_dir}")
            return None, "summary_stats_missing"
        try:
            stats_df = self._parse_summary_stats(stat_files[0])
            # Named region rows at fsize=0; exclude the all-regions aggregate row
            on_target_df = stats_df[
                (stats_df["family_size"] == 0)
                & stats_df["name"].notna()
                & (stats_df["name"] != "")
                & (~stats_df["regionid"].isin(["all_regions", "All", "all"]))
            ]
            return int(on_target_df["total_reads"].sum()), "region_bed"
        except Exception as e:
            logger.warning(f"Failed to parse summary stats for {sample_dir.name}: {e}")
            return None, "error"

    def _on_target_from_cons(
        self,
        sample_dir: Path,
        mutation_positions: set[tuple[str, int]],
    ) -> int | None:
        """Sum Coverage at mutation positions from _cons.tsv at Consensus group size=0."""
        cons_files = list(sample_dir.glob(f"*{CONSENSUS_TSV_SUFFIX}"))
        if not cons_files:
            logger.warning(f"No {CONSENSUS_TSV_SUFFIX} found in {sample_dir}")
            return None
        try:
            df = pd.read_csv(cons_files[0], sep="\t")
            df = df[df["Consensus group size"] == 0]
            mut_df = pd.DataFrame(list(mutation_positions), columns=["Contig", "Position"])
            merged = df.merge(mut_df, on=["Contig", "Position"], how="inner")
            return int(merged["Coverage"].sum())
        except Exception as e:
            logger.warning(f"Failed to compute on-target from cons for {sample_dir.name}: {e}")
            return None
