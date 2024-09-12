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

## Installation

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/ViralMet-Pipeline.git
   cd ViralMet-Pipeline
   ```

2. Create and activate the Conda environment:
   ```
   conda env create -f environment.yml
   conda activate viralmet
   ```

3. Download and set up necessary databases:
   ```
   bash setup_databases.sh
   ```

## Usage

1. Activate the Conda environment:
   ```
   conda activate viralmet
   ```

2. Run the pipeline:
   ```
   bash viralmet-pipeline.sh
   ```

   You can modify the parameters in the script or pass them as command-line arguments.

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

## Troubleshooting

If you encounter any issues, please check the following:

1. Ensure all required software is installed and in your PATH.
2. Check that you have sufficient disk space and memory.
3. Verify that your input files are in the correct format and location.

For more detailed error messages, check the log files in the output directory.

## Contributing

We welcome contributions to ViralMet-Pipeline!

## License

This project is licensed under the PANDASIA project License. Other purposes using is prohibited.
