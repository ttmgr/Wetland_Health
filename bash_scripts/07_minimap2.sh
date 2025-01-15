#!/bin/bash

# Default values
THREADS=16
READS_DIR=""
ASSEMBLY_DIR=""
OUTPUT_DIR="processing/minimap2"

# Function to print usage
usage() {
    echo "Usage: $0 [-t threads] [-r reads_dir] [-a assembly_dir] [-o output_dir]"
    echo ""
    echo "Options:"
    echo "  -t    Number of threads (default: 16)"
    echo "  -r    Directory containing reads (required)"
    echo "  -a    Directory containing assemblies (required)"
    echo "  -o    Output directory (default: processing/minimap2)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "t:r:a:o:h" opt; do
    case $opt in
        t) THREADS="$OPTARG" ;;
        r) READS_DIR="$OPTARG" ;;
        a) ASSEMBLY_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required parameters
if [ -z "$READS_DIR" ] || [ -z "$ASSEMBLY_DIR" ]; then
    echo "Error: Both reads directory (-r) and assembly directory (-a) are required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each assembly
for assembly_dir in "$ASSEMBLY_DIR"/*; do
    if [ -d "$assembly_dir" ]; then
        base_name=$(basename "$assembly_dir")
        assembly_file="${assembly_dir}/assembly.fasta"
        reads_file="${READS_DIR}/filtered_${base_name}.fastq"
        
        if [ -f "$assembly_file" ] && [ -f "$reads_file" ]; then
            echo "Processing $base_name..."
            
            minimap2 -ax map-ont \
                     -t "$THREADS" \
                     "$assembly_file" \
                     "$reads_file" > \
                     "${OUTPUT_DIR}/${base_name}.sam"
        fi
    fi
done

echo "Minimap2 alignment completed."
