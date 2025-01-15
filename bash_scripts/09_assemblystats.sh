#!/bin/bash

# Default values
INPUT_DIR="processing/racon"
OUTPUT_DIR="processing/assemblystats"

# Function to print usage
usage() {
    echo "Usage: $0 [-i input_dir] [-o output_dir]"
    echo ""
    echo "Options:"
    echo "  -i    Input directory containing FASTA files (default: processing/racon)"
    echo "  -o    Output directory (default: processing/assemblystats)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
    esac
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each FASTA file
for fasta_file in "$INPUT_DIR"/*.fasta; do
    if [ -f "$fasta_file" ]; then
        base_name=$(basename "$fasta_file" .fasta)
        echo "Generating assembly stats for $fasta_file..."
        assembly-stats "$fasta_file" > "${OUTPUT_DIR}/${base_name}_assemblystats.txt"
    fi
done

echo "Assembly statistics generation completed."
