# ViralMet-Pipeline: Viral Metagenomic Analysis Pipeline

ViralMet-Pipeline is a comprehensive pipeline for analyzing viral metagenomic data. It includes quality control, host removal, assembly, taxonomic classification, and functional annotation steps.

## Table of Contents
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Output](#output)
6. [Troubleshooting](#troubleshooting)
7. [Contributing](#contributing)
8. [License](#license)

## Overview

ViralMet-Pipeline performs the following steps:
1. Quality control (FastQC)
2. Trimming (Trimmomatic)
3. Host removal (Bowtie2)
4. Digital normalization (khmer)
5. Assembly (SPAdes)
6. BLAST and DIAMOND searches
7. Read mapping and coverage assessment
8. Multiple sequence alignment (MAFFT)
9. Phylogenetic tree reconstruction (IQ-TREE)
10. Taxonomic classification (Kraken2)
11. Functional annotation (Prokka)

## Requirements

- Linux-based operating system
- Conda (Miniconda or Anaconda)
- At least 16GB RAM (32GB or more recommended)
- At least 100GB free disk space


# ViralMet Pipeline Usage Guide

## Prerequisites
- Conda installed on your system
- Git (for cloning the repository)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/sittiton/ViralMet.git
cd ViralMet
```

2. Run the installation script:
```bash
bash install.sh
```

## Running the Pipeline

1. Activate the Conda environment:
```bash
conda activate viralmet
```

2. Create a configuration file named `config.ini` in the same directory as the pipeline script:
```ini
# Input/Output
INPUT_DIR="raw_data"
OUTPUT_DIR="analysis_results"
R1_FILE="sample_R1.fastq"
R2_FILE="sample_R2.fastq"

# Resources
THREADS=4
MEMORY=16

# Reference Files
REF_HOST="/path/to/host_reference.fasta"
REF_VIRUS_NUCL="/path/to/virus_nucleotide.fasta"
REF_VIRUS_PROT="/path/to/virus_protein.fasta"
ADAPTER_FILE="/path/to/adapters.fa"
KRAKEN_DB="/path/to/kraken_db"

# Pipeline Options
CLEAN_INTERMEDIATE=false
```

3. Run the pipeline:
```bash
bash viralmet-pipeline_1.0.1.sh
```

## Troubleshooting

If you encounter issues with Bowtie2, the pipeline will automatically handle environment switching. No manual intervention is required.


## Output

The pipeline will create an output directory with the following structure:

```
output_directory/
├── fastqc_results/
├── trimmed_reads/
├── host_removed_reads/
├── assembly/
├── blast_results/
├── diamond_results/
├── coverage/
├── aligned_contigs.fasta
├── phylogenetic_tree/
├── kraken_results/
└── prokka_results/
```

For more detailed error messages, check the log files in the output directory.

## Contributing

We welcome contributions to ViralMet-Pipeline!

## License

This project is licensed under the PANDASIA project License. Other purposes using is prohibited.
