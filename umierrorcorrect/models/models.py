from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, model_validator

from umierrorcorrect.core.check_args import is_tool
from umierrorcorrect.core.utils import check_output_directory, get_sample_name


class PreprocessConfig(BaseModel):
    """Pydantic model for preprocessing."""
    read1: Path
    read2: Optional[Path] = None
    output_path: Path
    umi_length: int
    spacer_length: int = 0
    sample_name: Optional[str] = None
    num_threads: int = 1
    dual_index: bool = False
    reverse_index: bool = False
    adapter_trimming: bool = False
    adapter_sequence: str = "illumina"
    force: bool = False
    tmpdir: Optional[Path] = None

    # Derived fields
    mode: Literal["single", "paired"] = "single"
    gziptool: Literal["gzip", "pigz"] = "gzip"

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode='after')
    def validate_and_configure(self) -> 'PreprocessConfig':
        # Check output directory
        self.output_path = Path(check_output_directory(str(self.output_path)))

        # Check environment tools
        is_pigz = is_tool("pigz")
        is_gzip = is_tool("gzip")
        is_bwa = is_tool("bwa")

        if self.adapter_trimming:
            if not is_tool("cutadapt"):
                raise ValueError('Cannot find program "cutadapt". Please install it and add it to the path.')

        if not is_bwa:
            raise ValueError('Cannot find program "bwa". Please install it and add it to the path.')

        if is_pigz:
            self.gziptool = "pigz"
        elif is_gzip:
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
                    if f1file.exists(): f1file.unlink()
                    if f2file.exists(): f2file.unlink()
        else:
            f1file = self.output_path / f"{self.sample_name}_umis_in_header.fastq.gz"
            if f1file.is_file():
                if not self.force:
                    raise ValueError(f"The file {f1file} already exists. Overwrite by setting force=True")
                else:
                    f1file.unlink()

        return self
