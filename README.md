# UMIErrorCorrect

Pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI).

## Reference

UMIErrorCorrect has been published in Clinical Chemistry:

> Osterlund T., Filges S., Johansson G., Stahlberg A. *UMIErrorCorrect and UMIAnalyzer: Software for Consensus Read Generation, Error Correction, and Visualization Using Unique Molecular Identifiers*, Clinical Chemistry, 2022, hvac136

[Link to the paper](https://doi.org/10.1093/clinchem/hvac136)

## Installation

### Using pip

```bash
pip install umierrorcorrect
```

### Using uv (recommended)

```bash
uv pip install umierrorcorrect
```

### Using Docker

```bash
docker pull ghcr.io/sfilges/umierrorcorrect:latest
# Or build locally
docker build -t umierrorcorrect docker/
```

See the [Docker documentation](doc/docker.md) for more details.

### Verify installation

```bash
run_umierrorcorrect --help
```

## Dependencies

UMIErrorCorrect requires Python 3.8+ and the following:

**Python libraries** (installed automatically):

- pysam (≥0.8.4)
- scipy
- matplotlib

**External programs** (must be in PATH):

- `bwa` - for read mapping
- `gzip` or `pigz` - for compression

### Reference genome setup

The pipeline uses `bwa` for mapping, so you need an indexed reference genome:

```bash
bwa index -a bwtsw reference.fa
```

## Usage

### Full pipeline

```bash
run_umierrorcorrect -r1 read1.fastq.gz -r2 read2.fastq.gz \
    -ul <umi_length> -sl <spacer_length> \
    -r reference.fa -o output_directory
```

### Pipeline steps

The `run_umierrorcorrect` command performs:

1. **Preprocessing** - Extract UMI from reads and add to header
2. **Mapping** - Align reads to reference genome with BWA
3. **UMI clustering** - Group reads by UMI similarity
4. **Error correction** - Generate consensus for each UMI cluster
5. **Statistics** - Create consensus output file with per-position counts
6. **Variant calling** - Call variants from consensus data

### Running individual steps

Each step can be run independently:

```bash
preprocess --help
run_mapping --help
umi_error_correct --help
get_consensus_statistics --help
call_variants --help
filter_bam --help
filter_cons --help
```

## Documentation

- [Tutorial](https://github.com/stahlberggroup/umierrorcorrect/wiki/Tutorial)
- [UMI definition options](https://github.com/stahlberggroup/umierrorcorrect/wiki/UMI-definition-options)
