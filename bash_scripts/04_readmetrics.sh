#!/bin/bash

# Default values
INPUT_DIR="processing/nanofilt"
OUTPUT_DIR="processing/nanostat_metrics"
THREADS=8

# Function to print usage
usage() {
    echo "Usage: $0 [-i input_dir] [-o output_dir] [-t threads] [-h]"
    echo ""
    echo "Options:"
    echo "  -i    Input directory containing FASTQ files (default: processing/nanofilt)"
    echo "  -o    Output directory (default: processing/nanostat_metrics)"
    echo "  -t    Number of threads (default: 8)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "i:o:t:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Generate NanoStat metrics for each FastQ file
for fastq in "$INPUT_DIR"/*.fastq; do
    if [ -f "$fastq" ]; then
        base=$(basename "$fastq" .fastq)
        echo "Generating NanoStat metrics for $fastq..."
        
        NanoStat --fastq "$fastq" \
                 --outdir "$OUTPUT_DIR" \
                 --name "${base}_nanostat_report" \
                 --threads "$THREADS"
    fi
done

echo "All tasks completed successfully. NanoStat reports are in $OUTPUT_DIR."
