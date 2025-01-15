#!/bin/bash

# Default values
INPUT_DIR="processing/racon"
OUTPUT_DIR="processing/amrfinder_racon"
DATABASE_PATH=""
THREADS=4
THREADS_PER_JOB=5

# Function to print usage
usage() {
    echo "Usage: $0 [-i input_dir] [-o output_dir] [-d database_path] [-t threads] [-p threads_per_job]"
    echo ""
    echo "Options:"
    echo "  -i    Input directory containing FASTA files (default: processing/racon)"
    echo "  -o    Output directory (default: processing/amrfinder_racon)"
    echo "  -d    Path to AMRFinder database (required)"
    echo "  -t    Number of parallel jobs (default: 4)"
    echo "  -p    Threads per job (default: 5)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "i:o:d:t:p:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) DATABASE_PATH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) THREADS_PER_JOB="$OPTARG" ;;
    esac
done

# Check required parameters
if [ -z "$DATABASE_PATH" ]; then
    echo "Error: Database path (-d) is required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Export variables for use in parallel command environment
export OUTPUT_DIR DATABASE_PATH THREADS_PER_JOB

# Process files in parallel using GNU Parallel
find "$INPUT_DIR" -type f -name "*.fasta" | parallel -j "$THREADS" '
    file="{}"
    output_file="$OUTPUT_DIR/$(basename "${file}")_amrfinder.txt"
    
    if [ ! -f "$output_file" ]; then
        echo "Processing ${file}..."
        amrfinder -n "${file}" \
                  -d "$DATABASE_PATH" \
                  --threads "$THREADS_PER_JOB" \
                  --plus > "$output_file"
    else
        echo "Output for ${file} already exists, skipping."
    fi
'

echo "AMRFinder Plus analysis of assemblies completed."
