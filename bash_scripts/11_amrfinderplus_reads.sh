#!/bin/bash

# Default values
INPUT_DIR="processing/fasta_output"
OUTPUT_DIR="processing/amrfinder_reads"
DATABASE_PATH=""
THREADS=20

# Function to print usage
usage() {
    echo "Usage: $0 [-i input_dir] [-o output_dir] [-d database_path] [-t threads]"
    echo ""
    echo "Options:"
    echo "  -i    Input directory containing FASTA files (default: processing/fasta_output)"
    echo "  -o    Output directory (default: processing/amrfinder_reads)"
    echo "  -d    Path to AMRFinder database (required)"
    echo "  -t    Number of threads (default: 20)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "i:o:d:t:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) DATABASE_PATH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
    esac
done

# Check required parameters
if [ -z "$DATABASE_PATH" ]; then
    echo "Error: Database path (-d) is required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run AMRFinder Plus on each FASTA file
for fasta_file in "$INPUT_DIR"/*.fasta; do
    if [ -f "$fasta_file" ]; then
        base_name=$(basename "$fasta_file" .fasta)
        output_file="${OUTPUT_DIR}/${base_name}_amrfinder.txt"
        echo "Running AMRFinder Plus on $fasta_file..."
        amrfinder --threads "$THREADS" \
                  -n "$fasta_file" \
                  -d "$DATABASE_PATH" \
                  --plus > "$output_file"
    fi
done

echo "AMRFinder Plus analysis of reads completed."
