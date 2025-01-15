#!/bin/bash

# Function to display usage information
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

A metagenomics analysis pipeline for nanopore sequencing data.

Required Options:
    -i, --input-dir DIR          Directory containing input FASTQ files
    -o, --output-dir DIR         Directory for output files
    -k, --kraken-db DIR          Path to Kraken2 database
    -a, --amrfinder-db DIR       Path to AMRFinder database

Optional Parameters:
    -t, --threads NUM            Number of threads (default: 28)
    -m, --min-length NUM         Minimum read length for filtering (default: 100)
    -c, --concat PAIRS           Barcode pairs to concatenate (optional)
                                Format: barcode01,barcode02:barcode03,barcode04
                                Multiple pairs separated by colons
    --mem MEM                    Memory limit (default: 250G)
    --time TIME                  Time limit (default: 50:00:00)
    -h, --help                   Display this help message

Examples:
# Run without concatenation:
$(basename "$0") -i /path/to/fastq -o /path/to/output -k /path/to/kraken_db -a /path/to/amrfinder_db

# Run with barcode concatenation:
$(basename "$0") -i /path/to/fastq -o /path/to/output -k /path/to/kraken_db -a /path/to/amrfinder_db -c barcode01,barcode02:barcode03,barcode04
EOF
}

# Default values
THREADS=28
MIN_READ_LENGTH=100
MEMORY="250G"
TIME="50:00:00"
PAIRS=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -k|--kraken-db)
            KRAKEN2_DB="$2"
            shift 2
            ;;
        -a|--amrfinder-db)
            AMRFINDER_DB="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--min-length)
            MIN_READ_LENGTH="$2"
            shift 2
            ;;
        -c|--concat)
            PAIRS="$2"
            shift 2
            ;;
        --mem)
            MEMORY="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required parameters
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$KRAKEN2_DB" ]] || [[ -z "$AMRFINDER_DB" ]]; then
    echo "Error: Missing required parameters"
    usage
    exit 1
fi

# Validate directories and databases
for dir in "$INPUT_DIR" "$KRAKEN2_DB" "$AMRFINDER_DB"; do
    if [[ ! -d "$dir" ]]; then
        echo "Error: Directory does not exist: $dir"
        exit 1
    fi
done

# Create SLURM script
create_slurm_script() {
    local slurm_script="${OUTPUT_DIR}/run_pipeline.slurm"
    
    cat > "$slurm_script" << EOL
#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=${MEMORY}
#SBATCH -t ${TIME}
#SBATCH --nice=10000
#SBATCH --mail-type=ALL
#SBATCH -c ${THREADS}
#SBATCH --job-name=metagenomics
#SBATCH -o ${OUTPUT_DIR}/logs/%x_%j.out
#SBATCH -e ${OUTPUT_DIR}/logs/%x_%j.err

# Source the pipeline script
source "\$(dirname "\$0")/pipeline_functions.sh"

# Set parameters
export THREADS=$THREADS
export KRAKEN2_DB="$KRAKEN2_DB"
export AMRFINDER_DB="$AMRFINDER_DB"
export MIN_READ_LENGTH=$MIN_READ_LENGTH
export INPUT_DIR="$INPUT_DIR"
export OUTPUT_DIR="$OUTPUT_DIR"
export PAIRS="$PAIRS"

# Initialize conda
source "\$(conda info --base)/etc/profile.d/conda.sh"

# Run the pipeline
main
EOL
    chmod +x "$slurm_script"
}

# Create pipeline functions script
create_pipeline_script() {
    local pipeline_script="${OUTPUT_DIR}/pipeline_functions.sh"
    
    cat > "$pipeline_script" << 'EOL'
#!/bin/bash

# Exit on error
set -euo pipefail

# Logging function
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] $1" | tee -a "${OUTPUT_DIR}/logs/pipeline.log"
}

# Error handling
handle_error() {
    local line_no=$1
    local error_code=$2
    log "ERROR: Command failed at line ${line_no} with exit code ${error_code}"
    exit "${error_code}"
}

trap 'handle_error ${LINENO} $?' ERR

# Environment management
activate_env() {
    local env_name=$1
    log "Activating environment: $env_name"
    conda activate "$env_name"
}

deactivate_env() {
    log "Deactivating current environment"
    conda deactivate
}

# Create directory structure
create_directories() {
    local dirs=(
        "${OUTPUT_DIR}/logs"
        "${OUTPUT_DIR}"/{porechop,nanofilt,kraken2_{core,contigs_core},flye,minimap2,samtools,racon,mdbg,fasta_output,amrfinder_output_{nanofilt,racon}}
    )
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
    done
}

# Optional concatenation of barcode pairs
concatenate_pairs() {
    if [[ -n "$PAIRS" ]]; then
        log "Concatenating specified barcode pairs..."
        IFS=':' read -ra PAIR_LIST <<< "$PAIRS"
        for pair in "${PAIR_LIST[@]}"; do
            IFS=',' read -r barcode1 barcode2 <<< "$pair"
            
            local file1="${OUTPUT_DIR}/nanofilt/filtered_${barcode1}_passed.fastq"
            local file2="${OUTPUT_DIR}/nanofilt/filtered_${barcode2}_passed.fastq"
            local output="${OUTPUT_DIR}/nanofilt/filtered_${barcode1}_${barcode2}_combined.fastq"
            
            if [[ -f "$file1" && -f "$file2" ]]; then
                log "Concatenating ${barcode1} and ${barcode2}..."
                cat "$file1" "$file2" > "$output"
            else
                log "Warning: One or both files missing for pair ${barcode1},${barcode2}"
            fi
        done
    fi
}

# Porechop processing
run_porechop() {
    activate_env porechop_env
    log "Starting Porechop processing..."
    
    for input in "${INPUT_DIR}"/*.fastq; do
        [[ -f "$input" ]] || continue
        local base_name=$(basename "$input" .fastq)
        local output="${OUTPUT_DIR}/porechop/trimmed_${base_name}.fastq"
        
        if [[ ! -s "$output" ]]; then
            log "Processing ${base_name} with Porechop..."
            porechop -i "$input" -o "$output" -t "$THREADS"
        else
            log "Porechop output exists for ${base_name}, skipping"
        fi
    done
    
    deactivate_env
}

# NanoFilt processing
run_nanofilt() {
    activate_env nanofilt_env
    log "Starting NanoFilt processing..."
    
    for input in "${OUTPUT_DIR}"/porechop/trimmed_*.fastq; do
        [[ -f "$input" ]] || continue
        local base_name=$(basename "$input" .fastq)
        local output="${OUTPUT_DIR}/nanofilt/filtered_${base_name}.fastq"
        
        if [[ ! -s "$output" ]]; then
            log "Filtering ${base_name} with NanoFilt..."
            cat "$input" | NanoFilt -l "$MIN_READ_LENGTH" > "$output"
        else
            log "NanoFilt output exists for ${base_name}, skipping"
        fi
    done
    
    deactivate_env
}

# Kraken2 analysis
run_kraken() {
    activate_env kraken2_env
    log "Running Kraken2 analysis..."
    
    for fastq_file in "${OUTPUT_DIR}"/nanofilt/*.fastq; do
        [[ -f "$fastq_file" ]] || continue
        local base_name=$(basename -- "$fastq_file" .fastq)
        local output_file="${OUTPUT_DIR}/kraken2_core/output_${base_name}.txt"
        
        if [[ ! -s "$output_file" ]]; then
            kraken2 --db "${KRAKEN2_DB}" --use-names \
                    --report "${OUTPUT_DIR}/kraken2_core/report_${base_name}.txt" \
                    --output "$output_file" "$fastq_file" --threads "$THREADS"
        fi
    done
    
    deactivate_env
}

# Assembly with Flye
run_assembly() {
    # Flye assembly
    activate_env flye_env
    log "Running Flye assembly..."
    
    for input in "${OUTPUT_DIR}"/nanofilt/*.fastq; do
        [[ -f "$input" ]] || continue
        local base_name=$(basename "$input" .fastq)
        
        if [[ ! -d "${OUTPUT_DIR}/flye/${base_name}" ]]; then
            flye --meta --nano-hq "$input" --threads "$THREADS" \
                 -o "${OUTPUT_DIR}/flye/${base_name}"
        fi
    done
    
    deactivate_env
    
    # Minimap2 alignment
    activate_env minimap2_env
    log "Running Minimap2 alignment..."
    
    for input in "${OUTPUT_DIR}"/nanofilt/*.fastq; do
        [[ -f "$input" ]] || continue
        local base_name=$(basename "$input" .fastq)
        local assembly="${OUTPUT_DIR}/flye/${base_name}/assembly.fasta"
        
        if [[ -f "$assembly" ]]; then
            if [[ ! -f "${OUTPUT_DIR}/minimap2/aligned_${base_name}.sam" ]]; then
                minimap2 -ax map-ont -t "$THREADS" "$assembly" "$input" > \
                        "${OUTPUT_DIR}/minimap2/aligned_${base_name}.sam"
            fi
        fi
    done
    
    deactivate_env
    
    # SAM to BAM conversion
    activate_env samtools_env
    log "Converting SAM to BAM..."
    
    for sam_file in "${OUTPUT_DIR}"/minimap2/*.sam; do
        [[ -f "$sam_file" ]] || continue
        local base_name=$(basename "$sam_file" .sam)
        
        if [[ ! -f "${OUTPUT_DIR}/samtools/${base_name}.bam" ]]; then
            samtools view -b -@ "$THREADS" "$sam_file" | \
                samtools sort -@ "$THREADS" -o "${OUTPUT_DIR}/samtools/${base_name}.bam"
        fi
    done
    
    deactivate_env
    
    # Racon polishing
    activate_env racon_env
    log "Running Racon polishing..."
    
    for input in "${OUTPUT_DIR}"/nanofilt/*.fastq; do
        [[ -f "$input" ]] || continue
        local base_name=$(basename "$input" .fastq)
        local sam_file="${OUTPUT_DIR}/minimap2/aligned_${base_name}.sam"
        local assembly="${OUTPUT_DIR}/flye/${base_name}/assembly.fasta"
        
        if [[ -f "$sam_file" && -f "$assembly" ]]; then
            if [[ ! -f "${OUTPUT_DIR}/racon/polished_${base_name}.fasta" ]]; then
                racon -t "$THREADS" "$input" "$sam_file" "$assembly" > \
                      "${OUTPUT_DIR}/racon/polished_${base_name}.fasta"
            fi
        fi
    done
    
    deactivate_env
}

# AMRFinder analysis
run_amrfinder() {
    activate_env amrfinderplus_env
    log "Running AMRFinder analysis..."
    
    # Process NanoFilt outputs
    for fastq_file in "${OUTPUT_DIR}"/nanofilt/*.fastq; do
        [[ -f "$fastq_file" ]] || continue
        local base_name=$(basename -- "$fastq_file" .fastq)
        local fasta_file="${OUTPUT_DIR}/fasta_output/${base_name}.fasta"
        local amr_output="${OUTPUT_DIR}/amrfinder_output_nanofilt/${base_name}_amrfinder.txt"
        
        if [[ ! -s "$amr_output" ]]; then
            seqtk seq -A "$fastq_file" > "$fasta_file"
            amrfinder --threads "$THREADS" -n "$fasta_file" \
                     -d "$AMRFINDER_DB" --plus > "$amr_output"
        fi
    done
    
    # Process Racon outputs
    for fasta_file in "${OUTPUT_DIR}"/racon/*.fasta; do
        [[ -f "$fasta_file" ]] || continue
        local base_name=$(basename -- "$fasta_file")
        local amr_output="${OUTPUT_DIR}/amrfinder_output_racon/${base_name}_amrfinder.txt"
        
        if [[ ! -s "$amr_output" ]]; then
            amrfinder -n "$fasta_file" -d "$AMRFINDER_DB" \
                     --threads "$THREADS" --plus > "$amr_output"
        fi
    done
    
    deactivate_env
}

# Kraken2 contig analysis
run_kraken_contigs() {
    activate_env kraken2_env
    log "Running Kraken2 contig analysis..."
    
    for fasta_file in "${OUTPUT_DIR}"/racon/*.fasta; do
        [[ -f "$fasta_file" ]] || continue
        local base_name=$(basename -- "$fasta_file" .fasta)
        local output_file="${OUTPUT_DIR}/kraken2_contigs_core/output_${base_name}.txt"
        
        if [[ ! -s "$output_file" ]]; then
            kraken2 --db "${KRAKEN2_DB}" --use-names \
                    --report "${OUTPUT_DIR}/kraken2_contigs_core/report_${base_name}.txt" \
                    --output "$output_file" "$fasta_file" --threads "$THREADS"
        fi
    done
    
    deactivate_env
}

# Create deduplication script
create_dedup_script() {
    cat > remove_duplicates.py << 'EOL'
import os
import gzip
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import logging

logging.basicConfig(level=logging.INFO)

def process_fastq_chunk(chunk):
    read_ids = defaultdict(list)
    for read in chunk:
        read_id_key = read[0].split()[0]
        read_ids[read_id_key].append(read)
    return read_ids

def check_and_remove_duplicates(fastq_file):
    logging.info(f"Processing {fastq_file}")
    read_ids = defaultdict(list)
    chunk_size = 10000  # Process reads in chunks
    current_chunk = []
    
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    
    with open_func(fastq_file, 'rt') as f:
        while True:
            read_id = f.readline().strip()
            if not read_id:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            current_chunk.append((read_id, seq, plus, qual))
            
            if len(current_chunk) >= chunk_size:
                with ProcessPoolExecutor() as executor:
                    chunk_results = executor.submit(process_fastq_chunk, current_chunk).result()
                    for key, reads in chunk_results.items():
                        read_ids[key].extend(reads)
                current_chunk = []

    # Process remaining reads
    if current_chunk:
        with ProcessPoolExecutor() as executor:
            chunk_results = executor.submit(process_fastq_chunk, current_chunk).result()
            for key, reads in chunk_results.items():
                read_ids[key].extend(reads)

    # Write deduplicated reads
    temp_output_file = fastq_file + '.tmp'
    with open_func(temp_output_file, 'wt') as out_f:
        for reads in read_ids.values():
            read = reads[0]  # Take first occurrence of each read
            out_f.write(f"{read[0]}\n{read[1]}\n{read[2]}\n{read[3]}\n")

    shutil.move(temp_output_file, fastq_file)
    logging.info(f"Completed processing {fastq_file}")

def process_directory():
    target_dir = "./processing/nanofilt"
    if not os.path.exists(target_dir):
        logging.error(f"Directory {target_dir} does not exist")
        return

    for file in os.listdir(target_dir):
        if file.endswith(('.fastq', '.fastq.gz')):
            fastq_file = os.path.join(target_dir, file)
            check_and_remove_duplicates(fastq_file)

if __name__ == "__main__":
    process_directory()
EOL
}

# Run deduplication
run_deduplication() {
    activate_env python_env
    log "Running deduplication..."
    create_dedup_script
    python remove_duplicates.py
    deactivate_env
}

# Main execution
main() {
    log "Starting pipeline execution..."
    create_directories
    run_porechop
    run_nanofilt
    concatenate_pairs
    run_deduplication
    run_kraken
    run_assembly
    run_amrfinder
    run_kraken_contigs
    log "Pipeline completed successfully"
}
EOL
    chmod +x "$pipeline_script"
}

# Setup and create scripts
mkdir -p "$OUTPUT_DIR/logs"
create_slurm_script
create_pipeline_script

echo "Pipeline scripts have been created in $OUTPUT_DIR"
echo "To submit the job, run: sbatch $OUTPUT_DIR/run_pipeline.slurm"
