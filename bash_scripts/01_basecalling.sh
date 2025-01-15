#!/bin/bash

# Default values
DORADO_PATH="dorado"
MODEL_PATH="dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
INPUT_DIR=""
OUTPUT_DIR="output"
KIT_NAME="SQK-RBK114-24"
THREADS=8

# Function to print usage
usage() {
    echo "Usage: $0 [-d dorado_path] [-m model_path] [-i input_dir] [-o output_dir] [-k kit_name] [-t threads]"
    echo ""
    echo "Options:"
    echo "  -d    Path to dorado binary (default: 'dorado' from PATH)"
    echo "  -m    Path to model file (default: dna_r10.4.1_e8.2_400bps_sup@v4.2.0)"
    echo "  -i    Input directory containing POD5/FAST5 files (required)"
    echo "  -o    Output directory (default: output)"
    echo "  -k    Kit name (default: SQK-RBK114-24)"
    echo "  -t    Number of threads (default: 8)"
    echo "  -h    Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "d:m:i:o:k:t:h" opt; do
    case $opt in
        d) DORADO_PATH="$OPTARG" ;;
        m) MODEL_PATH="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        k) KIT_NAME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check if input directory is provided
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory (-i) is required"
    usage
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Step 1: Run dorado basecaller
echo "Running dorado basecaller..."
"$DORADO_PATH" basecaller \
    "$MODEL_PATH" \
    -r "$INPUT_DIR" \
    --kit-name "$KIT_NAME" \
    --no-trim \
    --emit-fastq \
    --threads "$THREADS" \
    > "$OUTPUT_DIR/basecalled.fastq"

# Check if basecalling was successful
if [ $? -ne 0 ]; then
    echo "Error in basecalling step."
    exit 1
fi

# Step 2: Run dorado demux
echo "Running dorado demux..."
"$DORADO_PATH" demux \
    --output-dir "$OUTPUT_DIR/demux" \
    --kit-name "$KIT_NAME" \
    --threads "$THREADS" \
    "$OUTPUT_DIR/basecalled.fastq" \
    --emit-fastq

# Check if demultiplexing was successful
if [ $? -ne 0 ]; then
    echo "Error in demultiplexing step."
    exit 1
fi

echo "Process completed successfully."
