# Create a new conda environment for ViralMet-Pipeline
conda create -n viralmet python=3.8

# Activate the environment
conda activate viralmet

# Install required tools
conda install -c bioconda fastqc trimmomatic bowtie2 spades blast diamond samtools mafft iqtree kraken2 prokka

# Install khmer for digital normalization
conda install -c bioconda khmer

# Additional dependencies
conda install -c conda-forge pigz

# Kraken2 database setup (replace with appropriate database URL)
mkdir kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
tar -xzvf k2_standard_08gb_20230314.tar.gz -C kraken2_db

# Download and set up necessary reference databases (virus nucleotide, virus protein, host genome)
# Replace these URLs with appropriate database sources
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz -O ref_virus_nucl.fasta.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz -O ref_virus_prot.fasta.gz
gunzip ref_virus_nucl.fasta.gz ref_virus_prot.fasta.gz

# Download host reference genome (replace with appropriate host genome)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz -O reference_host.fasta.gz
gunzip reference_host.fasta.gz

# Download Trimmomatic adapters
wget https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE-2.fa
