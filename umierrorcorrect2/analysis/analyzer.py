from __future__ import annotations

import json
from pathlib import Path

from umierrorcorrect2.analysis.models import AnalysisSample, Mutation
from umierrorcorrect2.analysis.mutation_tracker import MutationTracker
from umierrorcorrect2.core.logging_config import get_logger

logger = get_logger("analysis")


class Analyzer:
    """Orchestrate the analysis of processed UMI data."""

    def __init__(self, family_size: int = 3):
        self.family_size = family_size
        self.tracker = MutationTracker(family_size=family_size)

    def analyze_sample(
        self,
        sample: AnalysisSample,
        output_dir: Path,
    ) -> None:
        """Run analysis for a single sample."""
        if not sample.mutation_bed:
            logger.warning(f"No mutation BED for sample {sample.name}, skipping mutation tracking.")
            return

        # Load mutations
        mutations = Mutation.load_from_bed(sample.mutation_bed)

        # Locate .cons file
        cons_file = output_dir / f"{sample.name}.cons"
        if not cons_file.exists():
            # Try searching for it if name differs
            cons_files = list(output_dir.glob("*.cons"))
            if cons_files:
                cons_file = cons_files[0]
            else:
                logger.error(f"Consensus file not found for {sample.name} in {output_dir}")
                return

        # Track mutations
        results = self.tracker.track_mutations(cons_file, mutations, sample.ml_plasma)
        sample.results = results

        # Create simplified cons file
        simplified_cons = output_dir / f"{sample.name}_patient_mutations.cons"
        self.tracker.create_simplified_cons(cons_file, mutations, simplified_cons)

        logger.info(f"Analysis complete for {sample.name}.")

    def save_results(self, sample: AnalysisSample, output_path: Path) -> None:
        """Save analysis results to a JSON file."""
        data = {
            "sample_name": sample.name,
            "results": [r.model_dump() for r in sample.results],
        }
        with output_path.open("w") as f:
            json.dump(data, f, indent=4)
