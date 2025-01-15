#!/bin/bash

# Default values
THREADS=16
INPUT_DIR=""
OUTPUT_DIR="processing/kraken2"
KRAKEN_DB=""

# Function to print usage
usage() {
    echo "Usage: $0 [-t threads] [-i input_dir] [-o output_dir] [-d kraken_db]"
    echo ""
    echo "Options:"
    echo "  -t    Number of threads (default: 16)"
    echo "  -i    Input directory containing FASTQ files (required)"
    echo "  -o    Output directory (default: processing/kraken2)"
    echo "  -d    Path to Kraken2 database (required)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "t:i:o:d:h" opt; do
    case $opt in
        t) THREADS="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) KRAKEN_DB="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required parameters
if [ -z "$INPUT_DIR" ] || [ -z "$KRAKEN_DB" ]; then
    echo "Error: Input directory (-i) and Kraken2 database (-d) are required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each FASTQ file
for file in "$INPUT_DIR"/*.fastq; do
    if [ -f "$file" ]; then
        base_name=$(basename "$file" .fastq)
        echo "Processing $file..."
        
        kraken2 --db "$KRAKEN_DB" \
                --threads "$THREADS" \
                --use-names \
                --output "${OUTPUT_DIR}/kraken2_${base_name}.out" \
                --report "${OUTPUT_DIR}/kraken2_${base_name}.report" \
                "$file"
    fi
done

echo "Kraken2 processing completed."
