#!/bin/bash

# Default values
THREADS=16
INPUT_DIR=""
OUTPUT_DIR="processing/flye"

# Function to print usage
usage() {
    echo "Usage: $0 [-t threads] [-i input_dir] [-o output_dir]"
    echo ""
    echo "Options:"
    echo "  -t    Number of threads (default: 16)"
    echo "  -i    Input directory containing FASTQ files (required)"
    echo "  -o    Output directory (default: processing/flye)"
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

# Create base output directory
mkdir -p "$OUTPUT_DIR"

# Process each FASTQ file
for file in "$INPUT_DIR"/*.fastq; do
    if [ -f "$file" ]; then
        base_name=$(basename "$file" .fastq)
        out_dir="${OUTPUT_DIR}/${base_name}"
        echo "Processing $file..."
        
        # Create sample-specific output directory
        mkdir -p "$out_dir"
        
        flye --nano-hq "$file" \
             --threads "$THREADS" \
             --meta \
             --out-dir "$out_dir"
    fi
done

echo "Flye assembly completed."
