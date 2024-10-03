#!/bin/bash

# This script will be called by the main pipeline to handle Bowtie2 execution

run_bowtie2() {
    local current_env=$CONDA_DEFAULT_ENV
    
    echo "Switching to bowtie2 environment..."
    conda deactivate
    conda activate bowtie2_env || {
        echo "Creating bowtie2 environment..."
        conda create -n bowtie2_env -c bioconda bowtie2 -y
        conda activate bowtie2_env
    }
    
    # Execute Bowtie2 command
    "$@"
    local bowtie2_exit_code=$?
    
    echo "Switching back to $current_env environment..."
    conda activate "$current_env"
    
    return $bowtie2_exit_code
}

# Execute the wrapped command
run_bowtie2 "$@"
