#!/bin/bash

# Default values
MIN_LENGTH=100
INPUT_DIR=""
OUTPUT_DIR="processing/nanofilt"

# Function to print usage
usage() {
    echo "Usage: $0 [-l min_length] [-i input_dir] [-o output_dir]"
    echo ""
    echo "Options:"
    echo "  -l    Minimum read length (default: 100)"
    echo "  -i    Input directory containing FASTQ files (required)"
    echo "  -o    Output directory (default: processing/nanofilt)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "l:i:o:h" opt; do
    case $opt in
        l) MIN_LENGTH="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
    esac
done

# Check required parameters
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory (-i) is required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each FASTQ file
for file in "$INPUT_DIR"/*.fastq; do
    if [ -f "$file" ]; then
        base_name=$(basename "$file" .fastq)
        echo "Processing $file..."
        
        cat "$file" | \
            NanoFilt -l "$MIN_LENGTH" > \
            "${OUTPUT_DIR}/filtered_${base_name}.fastq"
    fi
done

echo "NanoFilt processing completed."
