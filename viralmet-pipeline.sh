#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status.

# ViralMet-Pipeline: Viral Metagenomic Analysis Pipeline

# Parameters (adjust as needed)
THREADS=$(nproc)  # Use all available cores
MEMORY=$(free -g | awk '/^Mem:/{print int($2*0.8)}')  # Use 80% of available memory
ADAPTER_FILE="TruSeq3-PE-2.fa"
REF_HOST="reference_host.fasta"
REF_VIRUS_NUCL="ref_virus_nucl.fasta"
REF_VIRUS_PROT="ref_virus_prot.fasta"
KRAKEN_DB="/path/to/kraken2_db"  # Update this path

# Input/Output Directories
INPUT_DIR="raw_data"
OUTPUT_DIR="analysis_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p $OUTPUT_DIR

# Function to check if a command exists
command_exists () {
    type "$1" &> /dev/null ;
}

# Check for required software
for cmd in fastqc trimmomatic bowtie2 spades blastn blastx diamond samtools mafft iqtree kraken2 prokka normalize-by-median.py
do
    if ! command_exists $cmd ; then
        echo "Error: $cmd is not installed or not in PATH" >&2
        exit 1
    fi
done

# Check for input files
if [ ! -f "${INPUT_DIR}/sample_R1.fastq" ] || [ ! -f "${INPUT_DIR}/sample_R2.fastq" ]; then
    echo "Error: Input FASTQ files not found in ${INPUT_DIR}" >&2
    exit 1
fi

echo "Starting ViralMet-Pipeline..."

# --- Quality Control (FastQC) ---
echo "Running FastQC..."
fastqc ${INPUT_DIR}/sample_R*.fastq -t $THREADS -o $OUTPUT_DIR

# --- Trimming (Trimmomatic) ---
echo "Trimming with Trimmomatic..."
trimmomatic PE \
    ${INPUT_DIR}/sample_R1.fastq ${INPUT_DIR}/sample_R2.fastq \
    ${OUTPUT_DIR}/trimmed_R1.fastq ${OUTPUT_DIR}/unpaired_R1.fastq \
    ${OUTPUT_DIR}/trimmed_R2.fastq ${OUTPUT_DIR}/unpaired_R2.fastq \
    ILLUMINACLIP:$ADAPTER_FILE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# --- Host Removal (Bowtie2) ---
echo "Removing host reads with Bowtie2..."
bowtie2-build $REF_HOST ${OUTPUT_DIR}/host_index
bowtie2 -x ${OUTPUT_DIR}/host_index -1 ${OUTPUT_DIR}/trimmed_R1.fastq -2 ${OUTPUT_DIR}/trimmed_R2.fastq \
    --very-sensitive-local --un-conc-gz ${OUTPUT_DIR}/virus_reads -S ${OUTPUT_DIR}/unmapped.sam

# --- Digital Normalization (khmer) ---
echo "Performing digital normalization..."
normalize-by-median.py -k 20 -C 20 -N 4 -x 2e9 -p \
    ${OUTPUT_DIR}/virus_reads_1.fq.gz ${OUTPUT_DIR}/virus_reads_2.fq.gz

# --- Assembly (SPAdes) ---
echo "Assembling contigs with SPAdes..."
spades.py --meta -t $THREADS -m $MEMORY \
    -1 ${OUTPUT_DIR}/virus_reads_1.fq.gz.keep -2 ${OUTPUT_DIR}/virus_reads_2.fq.gz.keep \
    -o ${OUTPUT_DIR}/assembly

# --- BLAST ---
echo "Running BLAST..."
makeblastdb -in $REF_VIRUS_NUCL -dbtype nucl -out ${OUTPUT_DIR}/ref_virus_nucl -parse_seqids
makeblastdb -in $REF_VIRUS_PROT -dbtype prot -out ${OUTPUT_DIR}/ref_virus_prot -parse_seqids

blastn -db ${OUTPUT_DIR}/ref_virus_nucl -query ${OUTPUT_DIR}/assembly/contigs.fasta -out ${OUTPUT_DIR}/blastn_results.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
blastx -db ${OUTPUT_DIR}/ref_virus_prot -query ${OUTPUT_DIR}/assembly/contigs.fasta -out ${OUTPUT_DIR}/blastx_results.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# --- DIAMOND ---
echo "Running DIAMOND..."
diamond makedb --in $REF_VIRUS_PROT -d ${OUTPUT_DIR}/ref_virus_prot
diamond blastx -d ${OUTPUT_DIR}/ref_virus_prot -q ${OUTPUT_DIR}/assembly/contigs.fasta -o ${OUTPUT_DIR}/diamond_results.txt -f 6

# --- Read Mapping and Coverage Assessment ---
echo "Mapping reads and assessing coverage..."
bowtie2-build ${OUTPUT_DIR}/assembly/contigs.fasta ${OUTPUT_DIR}/contigs_index
bowtie2 -x ${OUTPUT_DIR}/contigs_index -1 ${OUTPUT_DIR}/trimmed_R1.fastq -2 ${OUTPUT_DIR}/trimmed_R2.fastq -S ${OUTPUT_DIR}/mapped_reads.sam

samtools view -S -b ${OUTPUT_DIR}/mapped_reads.sam > ${OUTPUT_DIR}/mapped_reads.bam
samtools sort ${OUTPUT_DIR}/mapped_reads.bam -o ${OUTPUT_DIR}/mapped_reads_sorted.bam
samtools index ${OUTPUT_DIR}/mapped_reads_sorted.bam
samtools depth ${OUTPUT_DIR}/mapped_reads_sorted.bam > ${OUTPUT_DIR}/coverage.txt

# --- Multiple Sequence Alignment (MAFFT) ---
echo "Aligning contigs with MAFFT..."
mafft --auto --maxiterate 1000 --localpair ${OUTPUT_DIR}/assembly/contigs.fasta > ${OUTPUT_DIR}/aligned_contigs.fasta

# --- Phylogenetic Tree Reconstruction (IQ-TREE) ---
echo "Building phylogenetic tree with IQ-TREE..."
iqtree -s ${OUTPUT_DIR}/aligned_contigs.fasta -m MFP -nt AUTO -bb 1000 -alrt 1000 -o ${OUTPUT_DIR}/iqtree_results

# --- Taxonomic Classification (Kraken2) ---
echo "Performing taxonomic classification with Kraken2..."
if [ -d "$KRAKEN_DB" ]; then
    kraken2 --db $KRAKEN_DB --threads $THREADS \
        --paired ${OUTPUT_DIR}/trimmed_R1.fastq ${OUTPUT_DIR}/trimmed_R2.fastq \
        --output ${OUTPUT_DIR}/kraken2_output.txt --report ${OUTPUT_DIR}/kraken2_report.txt
else
    echo "Warning: Kraken2 database not found. Skipping taxonomic classification."
fi

# --- Functional Annotation (Prokka) ---
echo "Performing functional annotation with Prokka..."
prokka ${OUTPUT_DIR}/assembly/contigs.fasta --outdir ${OUTPUT_DIR}/prokka_results --prefix viral_annotation

echo "ViralMet-Pipeline completed successfully."
