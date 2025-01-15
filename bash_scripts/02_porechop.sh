#!/bin/bash

# Default values
THREADS=16
INPUT_DIR=""
OUTPUT_DIR="processing/porechop"

# Function to print usage
usage() {
    echo "Usage: $0 [-t threads] [-i input_dir] [-o output_dir]"
    echo ""
    echo "Options:"
    echo "  -t    Number of threads (default: 16)"
    echo "  -i    Input directory containing FASTQ files (required)"
    echo "  -o    Output directory (default: processing/porechop)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "t:i:o:h" opt; do
    case $opt in
        t) THREADS="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required parameters
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory (-i) is required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each file in directory
for i in {01..24}; do
    file="${INPUT_DIR}/barcode${i}.fastq"
    if [ -f "$file" ]; then
        echo "Processing barcode${i}..."
        porechop -i "$file" \
                 -o "${OUTPUT_DIR}/trimmed_barcode${i}.fastq" \
                 -t "$THREADS"
    else
        echo "File $file does not exist. Skipping."
    fi
done

echo "Porechop processing completed."
