# Docker documentation

If you have Docker installed, pull the Docker image from the GitHub Container Registry.

```bash
docker pull ghcr.io/sfilges/umierrorcorrect2:latest
```

Download a reference genome fasta file and mount the reference directory and data directory (including fastq files and BED files) to the docker container. The container is configured with `umierrorcorrect` as the entrypoint.

To view the help message:

```bash
docker run --rm -it ghcr.io/sfilges/umierrorcorrect2:latest --help
```

## Running the pipeline

Since the umierrorcorrect pipeline uses `bwa` for mapping reads, a bwa-indexed reference genome is required.

To run the full pipeline, use the `batch` command. Note that all file paths must be relative to the mapped volumes inside the container (e.g., `/data/...` or `/references/...`).

Example command:

```bash
docker run --rm -v /path/to/references/:/references/ -v /path/to/data/:/data/ -it ghcr.io/sfilges/umierrorcorrect2:latest \
    batch \
    -r1 /data/read1.fastq.gz \
    -r2 /data/read2.fastq.gz \
    -ul 19 \
    -sl 16 \
    -r /references/reference.fa \
    -o /data/output_directory
```

You can also run individual steps, for example `preprocess` or `consensus` instead of `batch`.
