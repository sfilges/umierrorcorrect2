from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from umierrorcorrect2.core.check_args import is_tool
from umierrorcorrect2.core.utils import check_output_directory, get_sample_name


@dataclass
class FastpResult:
    """Result of fastp preprocessing."""

    filtered_read1: Path | None = None
    filtered_read2: Path | None = None
    merged_reads: Path | None = None
    fastp_json: Path | None = None
    fastp_html: Path | None = None


class FastpConfig(BaseModel):
    """Configuration for fastp preprocessing."""

    enabled: bool = True
    phred_score: int = 20
    merge_reads: bool = True
    trim_adapters: bool = True
    threads: int = 4
    # UMI extraction settings
    umi_enabled: bool = True
    umi_length: int = 0
    umi_skip: int = 0  # spacer length
    umi_loc: Literal["read1", "read2", "per_read"] = "read1"
    keep_unmerged: bool = False  # If True AND merge_reads=True, output unmerged reads


# =============================================================================
# Mutation Models
# =============================================================================


class Mutation(BaseModel):
    """Patient-specific mutation to track (parsed from mutation BED file).

    The mutation BED file format is:
        chrom  start  end  name  ref  alt   patient_id

    Example:
        chr3  178936091  178936091  KRAS_G12D  G  A   P001
    """

    name: str
    chromosome: str
    position: int  # 0-based start position
    reference: str
    alternate: str

    @classmethod
    def from_bed_line(cls, line: str) -> Mutation:
        """Parse a single line from a mutation BED file."""
        fields = line.strip().split("\t")
        if len(fields) < 7:
            raise ValueError(
                f"Mutation BED line requires 7 fields (chrom, start, end, name, ref, alt, patient_id): {line}"
            )
        return cls(
            chromosome=fields[0], position=int(fields[1]), name=fields[3], reference=fields[4], alternate=fields[5]
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

    mutation_name: str  # Links to Mutation.name
    vaf: float | None = None
    mm_count: int | None = None  # mutant molecule count
    total_reads: int | None = None
    mm_per_ml: float | None = None  # requires ml_plasma from sample

    @field_validator("vaf")
    @classmethod
    def validate_vaf(cls, v: float | None) -> float | None:
        if v is not None and not 0.0 <= v <= 1.0:
            raise ValueError(f"VAF must be between 0 and 1, got {v}")
        return v


# =============================================================================
# Sample Models
# =============================================================================


# TODO: The consensus information should be added to the sample model after the main pipeline has been run.abs
#       Together with that and the mutation be information and other metdata, the MutationResult model can be
#       calculated. For example, the VAF can be calculated as the ratio of mutant molecule count to total reads at
#       the mutation site. divinding the mutant molecule count by the ml_plasma will give the mm_per_ml.
class Sample(BaseModel):
    """A single sequencing sample with FASTQ files and metadata."""

    name: str
    read1: Path
    read2: Path | None = None

    # Sample metadata
    sample_type: Literal["validation", "ctdna"] | None = None
    replicate: int | None = None
    ng_input: float | None = None
    ml_plasma: float | None = None
    region_bed: Path | None = None
    collection_date: str | None = None

    # Analysis results (populated after processing)
    results: list[MutationResult] = []

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_files_exist(self) -> Sample:
        """Validate that input files exist."""
        if not self.read1.exists():
            raise ValueError(f"Read1 file not found for sample {self.name}: {self.read1}")
        if self.read2 is not None and not self.read2.exists():
            raise ValueError(f"Read2 file not found for sample {self.name}: {self.read2}")
        return self

    @model_validator(mode="after")
    def validate_date(self) -> Sample:
        """Validate that the collection date is in the correct format."""
        if self.collection_date is not None:
            try:
                datetime.strptime(self.collection_date, "%Y-%m-%d")
            except ValueError:
                raise ValueError(
                    f"Collection date {self.collection_date} is not in the correct format (YYYY-MM-DD)."
                ) from None
        return self

    def get_result(self, mutation_name: str) -> MutationResult | None:
        """Get the result for a specific mutation by name."""
        for result in self.results:
            if result.mutation_name == mutation_name:
                return result
        return None


# =============================================================================
# Patient Model
# =============================================================================


class Patient(BaseModel):
    """A patient with tracked mutations and associated samples."""

    patient_id: str
    mutation_bed: Path
    mutations: list[Mutation] = []
    samples: list[Sample] = []

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_and_load_mutations(self) -> Patient:
        """Validate mutation_bed exists and load mutations if not already loaded."""
        if not self.mutation_bed.exists():
            raise ValueError(f"Mutation BED file not found for patient {self.patient_id}: {self.mutation_bed}")
        if not self.mutations:
            self.mutations = Mutation.load_from_bed(self.mutation_bed)
        return self

    def get_sample(self, sample_name: str) -> Sample | None:
        """Get a sample by name."""
        for sample in self.samples:
            if sample.name == sample_name:
                return sample
        return None

    def get_samples_by_type(self, sample_type: Literal["validation", "ctdna"]) -> list[Sample]:
        """Get all samples of a specific type."""
        return [s for s in self.samples if s.sample_type == sample_type]

    def get_samples_sorted_by_date(self) -> list[Sample]:
        """Get samples sorted by collection date (earliest first)."""
        dated = [(s, s.collection_date) for s in self.samples if s.collection_date]
        dated.sort(key=lambda x: x[1])
        return [s for s, _ in dated]


# =============================================================================
# Sample Sheet Model
# =============================================================================

# Required columns in sample sheet CSV
# Minimal columns needed to run the main analysis pipeline from fastq to UMI Error Correction
SAMPLESHEET_REQUIRED_COLUMNS = {"sample_name", "read1"}

# Optional columns with their types
# These are needed for downstream analysis and reporting
SAMPLESHEET_OPTIONAL_COLUMNS = {
    "read2": Path,  # Needed for paired-end reads
    "region_bed": Path,  # Needed for annotation of regions during UMI Error Correction
    "patient_id": str,  # Needed to link samples and mutations to patients
    "mutation_bed": Path,  # Needed to link samples to patient-specific mutations
    "sample_type": str,  # Defines what types of reports to generate for the sample
    "replicate": int,  # If replicates exist they should be handled downstream to provide a single report
    "ng_input": float,  # Input DNA concentration in ng
    "ml_plasma": float,  # Plasma volume in ml, used to calculate mutant molecules per mL plasma
    "collection_date": str,  # Collection date in YYYY-MM-DD format
}


# TODO: Need to make sure that downstream analysis can handle the optional columns and that the analysis can handle
#       missing optional columns.
# TODO: This might conflict with the sample sheet parse method from batch.py. The basic sample sheet for processing
#       should follow the batch.py format (i.e. just need the sample name and fastq_1 with optional fastq_2).That is enough
#       to run the analysis. The additional data is only needed for downstream analysis and reporting.
class SampleSheet(BaseModel):
    """Parses and validates a sample sheet CSV, grouping samples by patient.

    Expected CSV format:
        patient_id,mutation_bed,sample_name,read1,read2,sample_type,replicate,ng_input,ml_plasma,region_bed,collection_date
        P001,/path/mutations.bed,P001_T0,/path/R1.fq.gz,/path/R2.fq.gz,ctdna,1,10.5,2.0,/path/regions.bed,2024-01-15

    Required columns: sample_name, read1
    Optional columns: read2, patient_id, mutation_bed, sample_type, replicate, ng_input, ml_plasma, region_bed, collection_date
    """

    csv_path: Path
    patients: list[Patient] = []
    base_path: Path | None = None  # Base path for resolving relative paths

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_and_parse(self) -> SampleSheet:
        """Validate the CSV file and parse into Patient objects."""
        if not self.csv_path.exists():
            raise ValueError(f"Sample sheet not found: {self.csv_path}")

        if not self.patients:
            self.patients = self._parse_csv()

        return self

    def _resolve_path(self, path_str: str) -> Path:
        """Resolve a path, making it relative to base_path if not absolute."""
        path = Path(path_str)
        if not path.is_absolute() and self.base_path:
            path = self.base_path / path
        return path

    def _parse_csv(self) -> list[Patient]:
        """Parse the CSV file and return a list of Patient objects."""
        with self.csv_path.open(newline="") as f:
            reader = csv.DictReader(f)

            # Validate required columns
            if reader.fieldnames is None:
                raise ValueError("Sample sheet is empty or has no header")

            missing_cols = SAMPLESHEET_REQUIRED_COLUMNS - set(reader.fieldnames)
            if missing_cols:
                raise ValueError(f"Sample sheet missing required columns: {missing_cols}")

            # Group rows by patient_id
            patient_data: dict[str, dict] = {}

            for row_num, row in enumerate(reader, start=2):  # start=2 accounts for header
                patient_id = row["patient_id"].strip()
                if not patient_id:
                    raise ValueError(f"Row {row_num}: patient_id cannot be empty")

                mutation_bed = self._resolve_path(row["mutation_bed"].strip())

                # Initialize patient if first occurrence
                if patient_id not in patient_data:
                    patient_data[patient_id] = {
                        "patient_id": patient_id,
                        "mutation_bed": mutation_bed,
                        "samples": [],
                    }
                else:
                    # Validate mutation_bed is consistent for the same patient
                    if patient_data[patient_id]["mutation_bed"] != mutation_bed:
                        raise ValueError(
                            f"Row {row_num}: Inconsistent mutation_bed for patient {patient_id}. "
                            f"Expected {patient_data[patient_id]['mutation_bed']}, got {mutation_bed}"
                        )

                # Build sample data
                sample_data = {
                    "name": row["sample_name"].strip(),
                    "read1": self._resolve_path(row["read1"].strip()),
                }

                # Add optional fields if present and non-empty
                for col, col_type in SAMPLESHEET_OPTIONAL_COLUMNS.items():
                    if col in row and row[col].strip():
                        value = row[col].strip()
                        if col_type is Path:
                            sample_data[col] = self._resolve_path(value)
                        elif col_type is int:
                            sample_data[col] = int(value)
                        elif col_type is float:
                            sample_data[col] = float(value)
                        else:
                            sample_data[col] = value

                patient_data[patient_id]["samples"].append(sample_data)

        # Create Patient objects with their samples
        patients = []
        for pdata in patient_data.values():
            samples = [Sample(**s) for s in pdata["samples"]]
            patient = Patient(
                patient_id=pdata["patient_id"],
                mutation_bed=pdata["mutation_bed"],
                samples=samples,
            )
            patients.append(patient)

        return patients

    def get_patient(self, patient_id: str) -> Patient | None:
        """Get a patient by ID."""
        for patient in self.patients:
            if patient.patient_id == patient_id:
                return patient
        return None

    def get_all_samples(self) -> list[Sample]:
        """Get all samples across all patients."""
        samples = []
        for patient in self.patients:
            samples.extend(patient.samples)
        return samples

    @property
    def patient_count(self) -> int:
        return len(self.patients)

    @property
    def sample_count(self) -> int:
        return sum(len(p.samples) for p in self.patients)


class PreprocessConfig(BaseModel):
    """Pydantic model for preprocessing."""

    read1: Path
    read2: Path | None = None
    output_path: Path
    umi_length: int
    spacer_length: int = 0
    sample_name: str | None = None
    num_threads: int = 1
    dual_index: bool = False
    reverse_index: bool = False
    adapter_trimming: bool = False
    adapter_sequence: str = "illumina"
    force: bool = False
    tmpdir: Path | None = None
    fastp_config: FastpConfig | None = None  # Optional fastp configuration

    # Derived fields
    mode: Literal["single", "paired"] = "single"
    gziptool: Literal["gzip", "pigz"] = "gzip"

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_and_configure(self) -> PreprocessConfig:
        # Check output directory
        self.output_path = Path(check_output_directory(str(self.output_path)))

        # Determine if fastp will handle adapter trimming
        fastp_handles_adapters = (
            self.fastp_config is not None and self.fastp_config.enabled and self.fastp_config.trim_adapters
        )

        # Only require cutadapt if adapter trimming is needed and fastp won't handle it
        if self.adapter_trimming and not fastp_handles_adapters and not is_tool("cutadapt"):
            raise ValueError('Cannot find program "cutadapt". Please install it and add it to the path.')

        if not is_tool("bwa"):
            raise ValueError('Cannot find program "bwa". Please install it and add it to the path.')

        if is_tool("pigz"):
            self.gziptool = "pigz"
        elif is_tool("gzip"):
            self.gziptool = "gzip"
        else:
            raise ValueError('Cannot find program "gzip" or "pigz". Install one of them and add to the path.')

        # Determine mode
        if self.read2:
            self.mode = "paired"
        else:
            self.mode = "single"

        # Derive sample name if not provided
        if not self.sample_name:
            self.sample_name = get_sample_name(str(self.read1), self.mode)

        # Logical constraints
        if self.dual_index and self.mode != "paired":
            raise ValueError("Dual index can only be used when both an R1 and R2 file are supplied.")

        if self.reverse_index and self.mode != "paired":
            raise ValueError("Reverse index can only be used when both an R1 and R2 file are supplied.")

        # Check input file existence
        if not self.read1.is_file():
            raise ValueError(f"The file specified as r1 ({self.read1}) does not exist.")

        if self.mode == "paired" and self.read2 and not self.read2.is_file():
            raise ValueError(f"The file specified as r2 ({self.read2}) does not exist.")

        # Check for existing output files
        if self.mode == "paired":
            f1file = self.output_path / f"{self.sample_name}_R1_umis_in_header.fastq.gz"
            f2file = self.output_path / f"{self.sample_name}_R2_umis_in_header.fastq.gz"
            if f1file.is_file() or f2file.is_file():
                if not self.force:
                    raise ValueError(f"The file {f1file} already exists. Overwrite by setting force=True")
                else:
                    if f1file.exists():
                        f1file.unlink()
                    if f2file.exists():
                        f2file.unlink()
        else:
            f1file = self.output_path / f"{self.sample_name}_umis_in_header.fastq.gz"
            if f1file.is_file():
                if not self.force:
                    raise ValueError(f"The file {f1file} already exists. Overwrite by setting force=True")
                else:
                    f1file.unlink()

        return self


class UMIErrorCorrectConfig(BaseModel):
    """Configuration for UMI error correction and consensus generation."""

    # Required
    reference_file: Path
    output_path: Path

    # Optional (can be auto-detected)
    bam_file: Path | None = None
    sample_name: str | None = None
    bed_file: Path | None = None

    # Clustering parameters
    position_threshold: int = 20
    edit_distance_threshold: int = 1

    # Processing
    num_threads: int | None = None  # None = use cpu_count()
    include_singletons: bool = True
    consensus_method: Literal["position", "most_common", "msa"] = "position"
    regions_from_bed: bool = False
    regions_from_tag: bool = False

    # Frequency thresholds
    indel_frequency_threshold: float = 60.0
    consensus_frequency_threshold: float = 60.0

    # Output options
    output_json: bool = False
    remove_large_files: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_and_configure(self) -> UMIErrorCorrectConfig:
        """Auto-detect bam_file and sample_name if not provided."""
        self.output_path = Path(check_output_directory(str(self.output_path)))

        # Auto-detect bam_file if not provided
        if not self.bam_file and self.sample_name:
            testname = self.output_path / f"{self.sample_name}.sorted.bam"
            if testname.is_file():
                self.bam_file = testname

        if not self.bam_file:
            bamfiles = list(self.output_path.glob("*sorted.bam"))
            if len(bamfiles) == 1:
                self.bam_file = bamfiles[0]
            elif len(bamfiles) > 1:
                raise ValueError(
                    "Too many sorted.bam files in output folder, please specify "
                    "which sample to run with sample_name or bam_file parameter."
                )

        # Auto-detect sample_name from bam_file if not provided
        if self.bam_file and not self.sample_name:
            self.sample_name = get_sample_name(str(self.bam_file))

        return self
