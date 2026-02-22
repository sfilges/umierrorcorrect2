# Specs for the analysis module

The analysis module is located in umierrorcorrect2/analysis/

The main purpose is to provide functions and extended data models
for analysis of data processed with the core pipeline.

Data models are inherited from the base models and extended with analysis-
specific fields.

The actual processing logic is based on an Analyzer class that handles
the analysis post core processing (either right after or called
separately on already analyzed data).

- **Core:** Responsible for FASTQ -> BAM -> Consensus -> VCF.
- **Analysis:** Responsible for Consensus/VCF + Metadata -> Domain Report.


## Input files

For analysis samples and their processed data must be associated with rich metadata, which the the user provides via samplesheets.

The core analysis can be run with or without a baisc samplesheet, but the analysis needs _some_ additional data.

### Patient-specific mutations

For analysis of patient-specific mutations, an additional bed file must
be provided for each input sample (for an example see mutation_bed_template.bed).

Post core processing, patient specific mutations are retained in the cons file at the chosen cons level (the original file contains multiple rows per position at different cons levels). This creates a simplified cons file
with only the patient-specific mutations at the desired level. This simplified cons file is stored next to the original cons file.

**On-target calculation:** the sum of cons0 reads for each patient-specific
mutation represents the on-target reads. Calculate on-target fractions based on:
- **Global Efficiency:** On-target reads / Total reads in the FASTQ.
- **Enrichment Efficiency:** On-target reads / Total reads within targeted regions (if a regions BED was provided).

This step can be performed for all samples with a valid mutation bed file.

#### Reporting

A per-sample level report (html) is prepared that summarizes and visualizes the analysis. Items to report are:

- Sample name
- A table with patient-specific mutations:
    - Total reads (cons0), cons reads (chosen level), VAF (%, patient specific), Mutant molecules (mutated cons reads), Mutant molecules (all except patient-specific)
- On-target statistics:
    - Global on-target fraction
    - Enrichment on-target fraction (if regions BED provided)
- A bar chart with cons reads for each mutation




