from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from umierrorcorrect2.models.models import Sample


class Mutation(BaseModel):
    """Patient-specific mutation to track (parsed from mutation BED file).

    The mutation BED file format is:
        chrom  start  end  name  ref  alt

    Example:
        chr3  178936091  178936091  PIK3CA_p.E545K  G  A
    """

    name: str
    chromosome: str
    position: int  # 0-based start position
    reference: str
    alternate: str

    @classmethod
    def from_bed_line(cls, line: str) -> Mutation:
        """Parse a single line from a mutation BED file."""
        fields = line.strip().split("	")
        if len(fields) < 6:
            raise ValueError(f"Mutation BED line requires 6 fields (chrom, start, end, name, ref, alt): {line}")
        return cls(
            chromosome=fields[0],
            position=int(fields[1]),
            name=fields[3],
            reference=fields[4],
            alternate=fields[5],
        )

    @classmethod
    def load_from_bed(cls, bed_path: Path) -> list[Mutation]:
        """Load mutations from a BED file."""
        mutations = []
        with bed_path.open() as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    mutations.append(cls.from_bed_line(line))
        return mutations


class MutationResult(BaseModel):
    """Analysis result for a tracked mutation in a specific sample."""

    mutation_name: str
    vaf: float | None = None
    mm_count: int | None = None  # mutant molecule count
    total_reads: int | None = None
    mm_per_ml: float | None = None

    @field_validator("vaf")
    @classmethod
    def validate_vaf(cls, v: float | None) -> float | None:
        if v is not None and not 0.0 <= v <= 100.0:  # Allow 0-100 for percentage
            raise ValueError(f"VAF must be between 0 and 100, got {v}")
        return v


class AnalysisSample(Sample):
    """A sample extended with metadata for detailed analysis."""

    patient_id: str | None = None
    mutation_bed: Path | None = None
    sample_type: Literal["validation", "ctdna"] | None = None
    ml_plasma: float | None = None
    ng_input: float | None = None
    collection_date: str | None = None
    region_bed: Path | None = None

    # Results (populated after analysis)
    results: list[MutationResult] = []

    def get_result(self, mutation_name: str) -> MutationResult | None:
        for r in self.results:
            if r.mutation_name == mutation_name:
                return r
        return None


# Required columns for extended analysis
ANALYSIS_SAMPLESHEET_REQUIRED = {"sample_name", "read1", "mutation_bed"}

ANALYSIS_SAMPLESHEET_OPTIONAL = {
    "read2": Path,
    "patient_id": str,
    "sample_type": str,
    "ml_plasma": float,
    "ng_input": float,
    "collection_date": str,
    "region_bed": Path,
}


class AnalysisSampleSheet(BaseModel):
    """Parses a sample sheet for extended analysis."""

    csv_path: Path
    samples: list[AnalysisSample] = []
    base_path: Path | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_and_parse(self) -> AnalysisSampleSheet:
        if not self.csv_path.exists():
            raise ValueError(f"Sample sheet not found: {self.csv_path}")

        if not self.samples:
            self.samples = self._parse_csv()

        return self

    def _resolve_path(self, path_str: str) -> Path:
        path = Path(path_str)
        if not path.is_absolute() and self.base_path:
            path = self.base_path / path
        return path

    def _parse_csv(self) -> list[AnalysisSample]:
        samples = []
        with self.csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError("Sample sheet is empty")

            missing = ANALYSIS_SAMPLESHEET_REQUIRED - set(reader.fieldnames)
            if missing:
                raise ValueError(f"Missing required columns for analysis: {missing}")

            for row in reader:
                # Build a data dictionary for Pydantic validation
                data: dict[str, Any] = {
                    "name": row["sample_name"].strip(),
                    "read1": self._resolve_path(row["read1"].strip()),
                    "mutation_bed": self._resolve_path(row["mutation_bed"].strip()),
                }

                for col, col_type in ANALYSIS_SAMPLESHEET_OPTIONAL.items():
                    if col in row and row[col].strip():
                        val = row[col].strip()
                        if col_type is Path:
                            data[col] = self._resolve_path(val)
                        elif col_type is float:
                            data[col] = float(val)
                        else:
                            data[col] = val

                samples.append(AnalysisSample.model_validate(data))
        return samples
