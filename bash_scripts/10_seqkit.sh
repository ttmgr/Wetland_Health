#!/bin/bash

# Default values
INPUT_DIR="processing/nanofilt"
OUTPUT_DIR="processing/fasta_output"

# Function to print usage
usage() {
    echo "Usage: $0 [-i input_dir] [-o output_dir]"
    echo ""
    echo "Options:"
    echo "  -i    Input directory containing FASTQ files (default: processing/nanofilt)"
    echo "  -o    Output directory (default: processing/fasta_output)"
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

# Convert each FASTQ file to FASTA
for fastq_file in "$INPUT_DIR"/*.fastq; do
    if [ -f "$fastq_file" ]; then
        base_name=$(basename "$fastq_file" .fastq)
        echo "Converting $fastq_file to FASTA format..."
        seqkit fq2fa "$fastq_file" -o "${OUTPUT_DIR}/${base_name}.fasta"
    fi
done

echo "FASTQ to FASTA conversion completed."
