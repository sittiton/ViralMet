#!/bin/bash

# ViralMet-Pipeline: Viral Metagenomic Analysis Pipeline
# Version 1.1.0

set -eo pipefail  # Exit on error, pipe failure

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load configuration
CONFIG_FILE="${SCRIPT_DIR}/config.ini"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Configuration file $CONFIG_FILE not found." >&2
    echo "Please create a config.ini file with the required parameters." >&2
    exit 1
fi

# Set default values if not specified in config
: ${THREADS:=$(nproc)}
: ${MEMORY:=$(free -g | awk '/^Mem:/{print int($2*0.8)}')}
: ${INPUT_DIR:="raw_data"}
: ${OUTPUT_DIR:="analysis_results_$(date +%Y%m%d_%H%M%S)"}
: ${LOG_FILE:="${OUTPUT_DIR}/pipeline.log"}
: ${CLEAN_INTERMEDIATE:=false}

# Verify required parameters
required_params=(R1_FILE R2_FILE REF_HOST REF_VIRUS_NUCL REF_VIRUS_PROT ADAPTER_FILE)
for param in "${required_params[@]}"; do
    if [[ -z "${!param}" ]]; then
        echo "Error: Required parameter $param is not set in config.ini" >&2
        exit 1
    fi
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Logging function
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

# Error handling function
error_exit() {
    log "ERROR: $1" >&2
    exit 1
}

# Cleanup function
cleanup() {
    if [[ "$CLEAN_INTERMEDIATE" == "true" ]]; then
        log "Cleaning up intermediate files..."
        rm -rf "${OUTPUT_DIR}"/*.sam "${OUTPUT_DIR}"/*.bam
    fi
}

# Bowtie2 wrapper function
run_bowtie2() {
    log "Running Bowtie2 command with separate environment..."
    bash "${SCRIPT_DIR}/bowtie2_wrapper.sh" "$@"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to ensure we're in the correct conda environment
ensure_environment() {
    if [[ "$CONDA_DEFAULT_ENV" != "viralmet" ]]; then
        error_exit "Please activate the viralmet conda environment before running this script"
    fi
}

# Version checking function
check_version() {
    local cmd="$1"
    local required_version="$2"
    local version
    
    if ! version=$($cmd --version 2>&1 | head -n1 | grep -oP '(\d+\.)+\d+' | head -n1); then
        error_exit "Failed to get version for $cmd"
    fi
    
    if ! awk -v v1="$version" -v v2="$required_version" 'BEGIN{if (v1 < v2) exit 1; exit 0}'; then
        error_exit "$cmd version $version is lower than required version $required_version"
    fi
}

# Main pipeline
main() {
    ensure_environment
    
    log "Starting ViralMet-Pipeline v1.1.0..."

    # Check for required software and versions
    log "Checking required software..."
    required_tools=(
        "fastqc:0.11.9"
        "trimmomatic:0.39"
        "spades.py:3.15.0"
        "blastn:2.12.0"
        "diamond:2.0.15"
        "samtools:1.15"
        "mafft:7.505"
        "iqtree:2.2.0"
        "kraken2:2.1.2"
        "prokka:1.14.6"
        "normalize-by-median.py:3.0.0"
    )
    
    for tool_info in "${required_tools[@]}"; do
        IFS=':' read -r tool version <<< "$tool_info"
        if ! command_exists "$tool"; then
            error_exit "$tool is not installed or not in PATH"
        fi
        check_version "$tool" "$version"
    done

    # Check for input files
    log "Checking input files..."
    for file in "${INPUT_DIR}/${R1_FILE}" "${INPUT_DIR}/${R2_FILE}"; do
        if [[ ! -f "$file" ]]; then
            error_exit "Input file $file not found"
        fi
    done

    # Quality Control
    log "Running FastQC..."
    fastqc "${INPUT_DIR}"/${R1_FILE} "${INPUT_DIR}"/${R2_FILE} -t "$THREADS" -o "$OUTPUT_DIR" || error_exit "FastQC failed"

    # Trimming
    log "Trimming with Trimmomatic..."
    trimmomatic PE \
        "${INPUT_DIR}/${R1_FILE}" "${INPUT_DIR}/${R2_FILE}" \
        "${OUTPUT_DIR}/trimmed_R1.fastq" "${OUTPUT_DIR}/unpaired_R1.fastq" \
        "${OUTPUT_DIR}/trimmed_R2.fastq" "${OUTPUT_DIR}/unpaired_R2.fastq" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || error_exit "Trimmomatic failed"

    # Host Removal
    log "Removing host reads with Bowtie2..."
    run_bowtie2 bowtie2-build "$REF_HOST" "${OUTPUT_DIR}/host_index" || error_exit "Bowtie2 index building failed"
    run_bowtie2 bowtie2 -x "${OUTPUT_DIR}/host_index" \
        -1 "${OUTPUT_DIR}/trimmed_R1.fastq" -2 "${OUTPUT_DIR}/trimmed_R2.fastq" \
        --very-sensitive-local --un-conc-gz "${OUTPUT_DIR}/virus_reads" \
        -S "${OUTPUT_DIR}/unmapped.sam" || error_exit "Bowtie2 alignment failed"

    # Digital Normalization
    log "Performing digital normalization..."
    normalize-by-median.py -k 20 -C 20 -N 4 -x 2e9 -p \
        "${OUTPUT_DIR}/virus_reads_1.fq.gz" "${OUTPUT_DIR}/virus_reads_2.fq.gz" || error_exit "Digital normalization failed"

    # Assembly
    log "Assembling contigs with SPAdes..."
    spades.py --meta -t "$THREADS" -m "$MEMORY" \
        -1 "${OUTPUT_DIR}/virus_reads_1.fq.gz.keep" -2 "${OUTPUT_DIR}/virus_reads_2.fq.gz.keep" \
        -o "${OUTPUT_DIR}/assembly" || error_exit "SPAdes assembly failed"

    # BLAST analysis
    log "Running BLAST analysis..."
    makeblastdb -in "$REF_VIRUS_NUCL" -dbtype nucl -out "${OUTPUT_DIR}/ref_virus_nucl" -parse_seqids || error_exit "makeblastdb for nucleotides failed"
    makeblastdb -in "$REF_VIRUS_PROT" -dbtype prot -out "${OUTPUT_DIR}/ref_virus_prot" -parse_seqids || error_exit "makeblastdb for proteins failed"

    blastn -db "${OUTPUT_DIR}/ref_virus_nucl" -query "${OUTPUT_DIR}/assembly/contigs.fasta" \
        -out "${OUTPUT_DIR}/blastn_results.txt" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" || error_exit "blastn failed"

    blastx -db "${OUTPUT_DIR}/ref_virus_prot" -query "${OUTPUT_DIR}/assembly/contigs.fasta" \
        -out "${OUTPUT_DIR}/blastx_results.txt" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" || error_exit "blastx failed"

    # DIAMOND analysis
    log "Running DIAMOND analysis..."
    diamond makedb --in "$REF_VIRUS_PROT" -d "${OUTPUT_DIR}/ref_virus_prot" || error_exit "DIAMOND makedb failed"
    diamond blastx -d "${OUTPUT_DIR}/ref_virus_prot" -q "${OUTPUT_DIR}/assembly/contigs.fasta" \
        -o "${OUTPUT_DIR}/diamond_results.txt" -f 6 || error_exit "DIAMOND blastx failed"

    # Read Mapping and Coverage Assessment
    log "Mapping reads and assessing coverage..."
    run_bowtie2 bowtie2-build "${OUTPUT_DIR}/assembly/contigs.fasta" "${OUTPUT_DIR}/contigs_index" || error_exit "Bowtie2 contig index building failed"
    run_bowtie2 bowtie2 -x "${OUTPUT_DIR}/contigs_index" \
        -1 "${OUTPUT_DIR}/trimmed_R1.fastq" -2 "${OUTPUT_DIR}/trimmed_R2.fastq" \
        -S "${OUTPUT_DIR}/mapped_reads.sam" || error_exit "Bowtie2 read mapping failed"

    samtools view -S -b "${OUTPUT_DIR}/mapped_reads.sam" > "${OUTPUT_DIR}/mapped_reads.bam" || error_exit "Samtools view failed"
    samtools sort "${OUTPUT_DIR}/mapped_reads.bam" -o "${OUTPUT_DIR}/mapped_reads_sorted.bam" || error_exit "Samtools sort failed"
    samtools index "${OUTPUT_DIR}/mapped_reads_sorted.bam" || error_exit "Samtools index failed"
    samtools depth "${OUTPUT_DIR}/mapped_reads_sorted.bam" > "${OUTPUT_DIR}/coverage.txt" || error_exit "Samtools depth failed"

    # Multiple Sequence Alignment
    log "Aligning contigs with MAFFT..."
    mafft --auto --maxiterate 1000 --localpair "${OUTPUT_DIR}/assembly/contigs.fasta" > "${OUTPUT_DIR}/aligned_contigs.fasta" || error_exit "MAFFT alignment failed"

    # Phylogenetic Tree Reconstruction
    log "Building phylogenetic tree with IQ-TREE..."
    iqtree -s "${OUTPUT_DIR}/aligned_contigs.fasta" -m MFP -nt AUTO -bb 1000 -alrt 1000 -o "${OUTPUT_DIR}/iqtree_results" || error_exit "IQ-TREE failed"

    # Taxonomic Classification
    if [[ -n "$KRAKEN_DB" && -d "$KRAKEN_DB" ]]; then
        log "Performing taxonomic classification with Kraken2..."
        kraken2 --db "$KRAKEN_DB" --threads "$THREADS" \
            --paired "${OUTPUT_DIR}/trimmed_R1.fastq" "${OUTPUT_DIR}/trimmed_R2.fastq" \
            --output "${OUTPUT_DIR}/kraken2_output.txt" --report "${OUTPUT_DIR}/kraken2_report.txt" || error_exit "Kraken2 classification failed"
    else
        log "Warning: Kraken2 database not found or not specified. Skipping taxonomic classification."
    fi

    # Functional Annotation
    log "Performing functional annotation with Prokka..."
    prokka "${OUTPUT_DIR}/assembly/contigs.fasta" --outdir "${OUTPUT_DIR}/prokka_results" --prefix viral_annotation || error_exit "Prokka annotation failed"

    log "ViralMet-Pipeline completed successfully."
}

# Set up cleanup trap
trap cleanup EXIT

# Run the pipeline
main

log "Pipeline execution completed. Results can be found in ${OUTPUT_DIR}"
