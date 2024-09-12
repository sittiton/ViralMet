#!/bin/bash

   # Create a directory for databases
   mkdir -p databases
   cd databases

   # Download and extract Kraken2 database
   wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
   tar -xzvf k2_standard_08gb_20230314.tar.gz
   rm k2_standard_08gb_20230314.tar.gz

   # Download virus reference databases
   wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
   wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
   gunzip viral.1.1.genomic.fna.gz
   gunzip viral.1.protein.faa.gz

   # Download host reference genome (example: human)
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
   gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

   # Download Trimmomatic adapters
   wget https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE-2.fa

   cd ..
   