# DNA Shotgun Metagenomics Analysis Pipeline

This document provides detailed commands for the DNA shotgun metagenomics analysis workflow, from raw Nanopore data to AMR gene identification and taxonomic assignment.

**Prerequisites:**
* Ensure all necessary tools are installed as per the [`Installation_tutorial.md`](./Installation_tutorial.md).
* Raw POD5 data should be available.
* You should have consulted the [`dna_shotgun_data_guide.md`](./dna_shotgun_data_guide.md) for sample-specific information, barcoding details, and any pre-processing steps for raw data.

**Workflow Overview:**

1.  Basecalling and Demultiplexing
2.  Read Processing (Adapter Trimming, Quality/Length Filtering)
3.  Taxonomic Classification (Reads)
4.  Metagenome Assembly and Polishing
5.  Antimicrobial Resistance (AMR) Gene Detection
6.  Taxonomic Origin of AMR on Contigs

---

## 1. Basecalling and Demultiplexing (Dorado)

This step converts the raw Nanopore signal data (POD5 files) into FASTQ sequence files and separates reads based on barcodes. The manuscript used Dorado v5.0.0 with the `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` model and `SQK-RBK114-24` kit for DNA shotgun sequencing. [cite: 83]

**Important:**
* The commands below are generic. Replace placeholders like `/path/to/...`, `<your_input_pod5_dir>`, `<your_output_dir>`, etc., with your actual paths and desired names.
* Refer to the Dorado documentation and the [`dna_shotgun_data_guide.md`](./dna_shotgun_data_guide.md) for specific model paths and kit names if they differ.
* Ensure Dorado is installed and accessible in your environment.

```bash
# Define variables (modify these paths and names as needed)
DORADO_BIN="/path/to/dorado_executable/dorado" # e.g., /opt/dorado/bin/dorado
CONFIG_FILE="/path/to/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0" # Model dna_r10.4.1_e8.2_400bps_sup@v5.0.0 used
INPUT_POD5_DIR="/path/to/your_raw_pod5_data_recursive_search" # Dorado will search recursively
BASECALLED_FASTQ="basecalled_all_samples.fastq" # Intermediate file
DEMUX_OUTPUT_DIR="demultiplexed_fastq"
KIT_NAME="SQK-RBK114-24" # Rapid Barcoding Kit 114-24 (RBK114-24) used for DNA libraries

# Create output directory for demultiplexed files
mkdir -p ${DEMUX_OUTPUT_DIR}

# Step 1: Run Dorado Basecaller
# Basecalling performed on raw sequencing data (.pod5 format) in super-accuracy (SUP) mode.
# Sequencing runs lasted 24 hours with a 5kHz sampling frequency and a minimum read length of 20 bases.
# User script example used --no-trim.
echo "Running Dorado basecaller..."
${DORADO_BIN} basecaller \
    ${CONFIG_FILE} \
    ${INPUT_POD5_DIR} \
    --kit-name ${KIT_NAME} \
    --no-trim \
    --emit-fastq > ${BASECALLED_FASTQ}

# Check if basecalling was successful
if [ $? -ne 0 ]; then
    echo "Error in Dorado basecalling step. Please check the command and paths."
    exit 1
fi
echo "Basecalling completed."

# Step 2: Run Dorado Demux
# Dorado subsequently demultiplexed FASTQ files based on assigned barcodes.
echo "Running Dorado demux..."
${DORADO_BIN} demux \
    --output-dir ${DEMUX_OUTPUT_DIR}/ \
    --kit-name ${KIT_NAME} \
    ${BASECALLED_FASTQ} \
    --emit-fastq

# Check if demultiplexing was successful
if [ $? -ne 0 ]; then
    echo "Error in Dorado demultiplexing step. Please check the command and paths."
    exit 1
fi
echo "Demultiplexing completed. Demultiplexed FASTQ files are in ${DEMUX_OUTPUT_DIR}/"
echo "Process completed successfully."
```
**Output:** Demultiplexed FASTQ files (one per barcode/sample) in the specified output directory. Each subsequent step will typically be run on these individual demultiplexed FASTQ files.

---

## You are absolutely right. My apologies for not adhering strictly to the "one tool, one command block" rule in the previous iteration. Let's revise that to ensure every distinct tool invocation or data manipulation step is clearly separated.

Here is the corrected dna_shotgun_analysis_pipeline.md with each command or logical tool step separated:
Markdown

# DNA Shotgun Metagenomics Analysis Pipeline

This document provides detailed commands for the DNA shotgun metagenomics analysis workflow, from raw Nanopore data to AMR gene identification and taxonomic assignment.

**Prerequisites:**
* Ensure all necessary tools are installed as per the [`Installation_tutorial.md`](./Installation_tutorial.md).
* Raw POD5 data should be available.
* You should have consulted the [`dna_shotgun_data_guide.md`](./dna_shotgun_data_guide.md) for sample-specific information, barcoding details, and any pre-processing steps for raw data.

**Workflow Overview:**

1.  Basecalling
2.  Demultiplexing
3.  Adapter and Barcode Trimming
4.  Read Quality and Length Filtering
5.  Read Downsampling (for PCoA)
6.  Taxonomic Classification of Reads
7.  Metagenome Assembly (Option A: metaFlye)
8.  Assembly Polishing (Option A: Racon - 3 rounds)
9.  Assembly Polishing (Option A: Medaka on metaFlye+Racon assembly)
10. Metagenome Assembly (Option B: nanoMDBG/metaMDBG)
11. Assembly Polishing (Option B: Medaka on nanoMDBG/metaMDBG assembly)
12. Read Preparation for AMR Detection (Downsampling and Format Conversion)
13. AMR Gene Detection (on Processed Reads)
14. AMR Gene Detection (on Assembled Contigs)
15. Taxonomic Classification of AMR-carrying Contigs (DIAMOND)
16. Taxonomic Classification of AMR-carrying Contigs (Kraken2)
17. Comparative Taxonomic Assignment of AMR Contigs

---

## 1. Basecalling (Dorado)

This step converts the raw Nanopore signal data (POD5 files) into FASTQ sequence files. The manuscript used Dorado v5.0.0 with the `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` model for DNA shotgun sequencing[cite: 1]. Basecalling was performed in super-accuracy (SUP) mode[cite: 1]. Sequencing runs lasted 24 hours with a 5kHz sampling frequency and a minimum read length of 20 bases[cite: 1].

**Important:**
* Ensure Dorado is installed and accessible in your environment.
* Replace placeholders with your actual paths and model names.

```bash
# Define variables for Basecalling
DORADO_BIN="/path/to/dorado_executable/dorado"
CONFIG_FILE="/path/to/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0" # [cite: 1]
INPUT_POD5_DIR="/path/to/your_raw_pod5_data_recursive_search"
BASECALLED_FASTQ="basecalled_all_samples.fastq" # Intermediate file for all reads before demux
KIT_NAME_DNA="SQK-RBK114-24" # Rapid Barcoding Kit 114-24 (RBK114-24) used for DNA libraries [cite: 1]

# Run Dorado Basecaller
echo "Running Dorado basecaller..."
${DORADO_BIN} basecaller \
    ${CONFIG_FILE} \
    ${INPUT_POD5_DIR} \
    --kit-name ${KIT_NAME_DNA} \
    --no-trim \
    --emit-fastq > ${BASECALLED_FASTQ}

# Check if basecalling was successful
if [ $? -ne 0 ]; then
    echo "Error in Dorado basecalling step. Please check the command and paths."
    exit 1
fi
echo "Basecalling completed."

Output: A single FASTQ file (basecalled_all_samples.fastq) containing all basecalled reads.
2. Demultiplexing (Dorado)

This step separates the basecalled reads from the combined FASTQ file into individual files based on their assigned barcodes. Dorado was used for this purpose.

Important:

    This step uses the output from the basecalling step.

Bash

# Define variables for Demultiplexing
DORADO_BIN="/path/to/dorado_executable/dorado" # Same as above
BASECALLED_FASTQ="basecalled_all_samples.fastq" # Input from previous step
DEMUX_OUTPUT_DIR="demultiplexed_fastq"
KIT_NAME_DNA="SQK-RBK114-24" # [cite: 1]

# Create output directory for demultiplexed files
mkdir -p ${DEMUX_OUTPUT_DIR}

# Run Dorado Demux
echo "Running Dorado demux..."
${DORADO_BIN} demux \
    --output-dir ${DEMUX_OUTPUT_DIR}/ \
    --kit-name ${KIT_NAME_DNA} \
    ${BASECALLED_FASTQ} \
    --emit-fastq

# Check if demultiplexing was successful
if [ $? -ne 0 ]; then
    echo "Error in Dorado demultiplexing step. Please check the command and paths."
    exit 1
fi
echo "Demultiplexing completed. Demultiplexed FASTQ files are in ${DEMUX_OUTPUT_DIR}/"
```

Output: Demultiplexed FASTQ files (one per barcode/sample) in the ${DEMUX_OUTPUT_DIR} directory. Each subsequent processing step will typically be run on these individual demultiplexed FASTQ files.

---

## 2. Adapter and Barcode Trimming (Porechop)

This step removes sequencing adapters and barcodes from the demultiplexed FASTQ files using Porechop v0.2.4.

Tool: Porechop (v0.2.4) 
Input: A single demultiplexed FASTQ file (e.g., barcodeXX.fastq) from the Dorado demultiplexing output.
Output: A trimmed FASTQ file (e.g., barcodeXX.trimmed.fastq).

```bash
# Activate Mamba environment for Porechop
# mamba activate porechop_env

# Define variables for a single sample
INPUT_DEMUX_FASTQ="demultiplexed_fastq/barcodeXX.fastq" # Replace barcodeXX with actual file
TRIMMED_FASTQ_DIR="processed_reads/trimmed"
TRIMMED_FASTQ_SAMPLE="${TRIMMED_FASTQ_DIR}/barcodeXX.trimmed.fastq"
THREADS=<N> # Number of threads to use, e.g., 4

# Create output directory
mkdir -p ${TRIMMED_FASTQ_DIR}

echo "Running Porechop for adapter/barcode trimming on: ${INPUT_DEMUX_FASTQ}"
porechop \
    -i ${INPUT_DEMUX_FASTQ} \
    -o ${TRIMMED_FASTQ_SAMPLE} \
    --threads ${THREADS}

if [ $? -ne 0 ]; then
    echo "Error in Porechop step for ${INPUT_DEMUX_FASTQ}."
    # exit 1 # Decide error handling
fi
echo "Porechop completed for ${INPUT_DEMUX_FASTQ}. Output: ${TRIMMED_FASTQ_SAMPLE}"

# Deactivate environment
# mamba deactivate
```

---

## 3. Quality and Length Filtering with NanoFilt
# Reads were filtered using NanoFilt v2.8.0, keeping only reads with a minimum length of 100 base pairs[cite: 85].

```bash
# Activate Mamba environment for NanoFilt
# mamba activate nanofilt_env

# Define variables for a single sample
TRIMMED_FASTQ_SAMPLE="processed_reads/trimmed/barcodeXX.trimmed.fastq" # Input from previous step
FILTERED_FASTQ_DIR="processed_reads/filtered"
FILTERED_FASTQ_SAMPLE="${FILTERED_FASTQ_DIR}/barcodeXX.filtered.fastq"

# Create output directory
mkdir -p ${FILTERED_FASTQ_DIR}

echo "Running NanoFilt for quality/length filtering on: ${TRIMMED_FASTQ_SAMPLE}"
NanoFilt \
    --length 100 \
    < ${TRIMMED_FASTQ_SAMPLE} \
    > ${FILTERED_FASTQ_SAMPLE}
    # Quality threshold not explicitly specified in the PDF for this NanoFilt step,
    # but could be added if desired (e.g., --quality 8).

if [ $? -ne 0 ]; then
    echo "Error in NanoFilt step for ${TRIMMED_FASTQ_SAMPLE}."
    # exit 1 # Decide error handling
fi
echo "NanoFilt completed for ${TRIMMED_FASTQ_SAMPLE}. Output: ${FILTERED_FASTQ_SAMPLE}"

# Deactivate environment
# mamba deactivate
```

*Repeat these steps for each demultiplexed FASTQ file.*

---

## 4. Read Downsampling (for PCoA) (Seqkit)

For Principal Coordinate Analysis (PCoA), filtered reads were randomly downsampled to 14,000 reads per sample using Seqkit v2.10.0. Samples with fewer than 14,000 reads were excluded from PCoA (e.g., Air Sample 1 mentioned in the manuscript ).

Tool: Seqkit (v2.10.0) 
Input: A filtered FASTQ file (e.g., barcodeXX.filtered.fastq).
Output: A downsampled FASTQ file (e.g., barcodeXX.14k.fastq).

```bash
# Activate Mamba environment for Seqkit
# mamba activate seqkit_env

# Define variables
FILTERED_FASTQ_SAMPLE="processed_reads/filtered/barcodeXX.filtered.fastq" # Input from previous step
DOWNSAMPLED_FASTQ_DIR="processed_reads/downsampled_for_pcoa"
DOWNSAMPLED_FASTQ_SAMPLE="${DOWNSAMPLED_FASTQ_DIR}/barcodeXX.14k.fastq"
PCOA_DOWNSAMPLE_THRESHOLD=14000 # [cite: 1]

# Create output directory
mkdir -p ${DOWNSAMPLED_FASTQ_DIR}

echo "Downsampling reads for PCoA from: <span class="math-inline">\{FILTERED\_FASTQ\_SAMPLE\}"
\# Optional\: Check read count before downsampling if you want to handle exclusions explicitly here
\# num\_reads\=</span>(seqkit stats -T ${FILTERED_FASTQ_SAMPLE} | awk 'NR==2 {print $4}' | sed 's/,//g')
# if [ "$num_reads" -ge "$PCOA_DOWNSAMPLE_THRESHOLD" ]; then
    seqkit sample -n ${PCOA_DOWNSAMPLE_THRESHOLD} -s 11 ${FILTERED_FASTQ_SAMPLE} -o ${DOWNSAMPLED_FASTQ_SAMPLE}
    echo "Downsampling completed. Output: ${DOWNSAMPLED_FASTQ_SAMPLE}"
# else
#    echo "Sample ${FILTERED_FASTQ_SAMPLE} has fewer than ${PCOA_DOWNSAMPLE_THRESHOLD} reads. Excluding from PCoA downsampling (as per methods)."
# fi

if [ $? -ne 0 ]; then # This check applies if seqkit sample was run
    echo "Error in Seqkit downsampling step for ${FILTERED_FASTQ_SAMPLE}."
    # exit 1 # Decide error handling
fi

# Deactivate environment
# mamba deactivate
```

---

## 5. Taxonomic Classification of Reads (Kraken2)

Taxonomic classification was performed on the downsampled 14,000 read subsets (for PCoA) using Kraken2 v2.1.2 with the nt_core database (accessed May 2025). For general taxonomic profiling not limited to PCoA, this step can be run on the full set of filtered reads from Step 4.

Tool: Kraken2 (v2.1.2) 
Input: A downsampled FASTQ file (e.g., barcodeXX.14k.fastq) or a full filtered FASTQ file.
Output: Kraken2 output and report files.

```bash
# Activate Mamba environment for Kraken2
# mamba activate kraken2_env

# Define variables
INPUT_FASTQ_FOR_KRAKEN="processed_reads/downsampled_for_pcoa/barcodeXX.14k.fastq" # Or "processed_reads/filtered/barcodeXX.filtered.fastq" for all reads
KRAKEN_DB_PATH="/path/to/your_kraken2_db/nt_core_May2025" # nt_core database, accessed May 2025 [cite: 1]
KRAKEN_OUTPUT_DIR="kraken2_output/pcoa_set" # Or a different dir for all_reads_set
KRAKEN_OUTPUT_SAMPLE="<span class="math-inline">\{KRAKEN\_OUTPUT\_DIR\}/barcodeXX\.kraken\_output\.txt"
KRAKEN\_REPORT\_SAMPLE\="</span>{KRAKEN_OUTPUT_DIR}/barcodeXX.kraken_report.txt"
THREADS=<N> # e.g., 8

# Create output directory
mkdir -p ${KRAKEN_OUTPUT_DIR}

echo "Running Kraken2 taxonomic classification on: ${INPUT_FASTQ_FOR_KRAKEN}"
kraken2 \
    --db ${KRAKEN_DB_PATH} \
    --threads ${THREADS} \
    --output ${KRAKEN_OUTPUT_SAMPLE} \
    --report ${KRAKEN_REPORT_SAMPLE} \
    ${INPUT_FASTQ_FOR_KRAKEN}

if [ $? -ne 0 ]; then
    echo "Error in Kraken2 step for ${INPUT_FASTQ_FOR_KRAKEN}."
    # exit 1 # Decide error handling
fi
echo "Kraken2 classification completed. Output: ${KRAKEN_OUTPUT_SAMPLE}, Report: ${KRAKEN_REPORT_SAMPLE}"

# Deactivate environment
# mamba deactivate
```

---

## 6. Metagenome Assembly (Option A: metaFlye)

De novo assemblies were generated using metaFlye v2.9.6.

Tool: metaFlye (v2.9.6) 
Input: A filtered FASTQ file (e.g., barcodeXX.filtered.fastq from Step 4).
Output: Raw metaFlye assembly directory.

```bash
# Activate Mamba environment for metaFlye
# mamba activate metaflye_env

# Define variables
FILTERED_FASTQ_SAMPLE="processed_reads/filtered/barcodeXX.filtered.fastq" # From Step 4
ASSEMBLY_DIR_METAFYLYE_RAW="assembly_output/barcodeXX_metaflye/flye_raw_assembly"
THREADS=<N> # e.g., 16

# Create parent output directory
mkdir -p assembly_output/barcodeXX_metaflye

echo "Running metaFlye assembly for ${FILTERED_FASTQ_SAMPLE}..."
flye \
    --nano-hq ${FILTERED_FASTQ_SAMPLE} \
    --out-dir ${ASSEMBLY_DIR_METAFYLYE_RAW} \
    --meta \
    --threads ${THREADS}

if [ $? -ne 0 ]; then echo "Error in metaFlye step for ${FILTERED_FASTQ_SAMPLE}."; exit 1; fi
echo "metaFlye assembly completed. Raw assembly in ${ASSEMBLY_DIR_METAFYLYE_RAW}"
```

---

## 7. Assembly Polishing (Option A: Racon - 3 rounds with Minimap2)

metaFlye assemblies were polished using Minimap2 v2.28 and three rounds of Racon v1.5.

Tools: Minimap2 (v2.28), Racon (v1.5), SAMtools (v1.17 for sorting/indexing) 
Input: Raw metaFlye assembly (assembly.fasta) and the corresponding filtered FASTQ file.
Output: Racon-polished assembly FASTA file.

```bash
# Activate Mamba environments
# mamba activate minimap2_dna_env # Ensure this is v2.28
# mamba activate racon_env
# mamba activate samtools_env

# Define variables
FILTERED_FASTQ_SAMPLE="processed_reads/filtered/barcodeXX.filtered.fastq" # From Step 4
INITIAL_ASSEMBLY_FASTA="assembly_output/barcodeXX_metaflye/flye_raw_assembly/assembly.fasta" # From Step 7
RACCON_POLISH_DIR="assembly_output/barcodeXX_metaflye/racon_polished"
THREADS=<N> # e.g., 16

# Create output directory
mkdir -p ${RACCON_POLISH_DIR}

POLISHED_FA_RAC_ROUNDS=$INITIAL_ASSEMBLY_FASTA
echo "Starting Racon polishing (3 rounds) for metaFlye assembly..."
for i in {1..3}; do
    echo "Racon polishing round <span class="math-inline">\{i\}\.\.\."
ALIGNED\_SAM\_FOR\_RACON\="</span>{RACCON_POLISH_DIR}/racon_round${i}.aligned.sam"
    ALIGNED_BAM_FOR_RACON="<span class="math-inline">\{RACCON\_POLISH\_DIR\}/racon\_round</span>{i}.aligned.bam"
    CURRENT_POLISHED_FA="<span class="math-inline">\{RACCON\_POLISH\_DIR\}/racon\_round</span>{i}.fasta"

    echo "  Running Minimap2 (v2.28) for Racon round ${i}..."
    minimap2 -ax map-ont ${POLISHED_FA_RAC_ROUNDS} ${FILTERED_FASTQ_SAMPLE} -t ${THREADS} > ${ALIGNED_SAM_FOR_RACON}
    if [ $? -ne 0 ]; then echo "Error in Minimap2 for Racon round ${i}."; exit 1; fi

    echo "  Sorting SAM to BAM with SAMtools (v1.17) for Racon round ${i}..."
    samtools sort -@ ${THREADS} -o ${ALIGNED_BAM_FOR_RACON} ${ALIGNED_SAM_FOR_RACON}
    if [ $? -ne 0 ]; then echo "Error in samtools sort for Racon round ${i}."; exit 1; fi
    rm ${ALIGNED_SAM_FOR_RACON} # Clean up SAM

    echo "  Indexing BAM with SAMtools (v1.17) for Racon round ${i}..."
    samtools index ${ALIGNED_BAM_FOR_RACON}
    if [ $? -ne 0 ]; then echo "Error in samtools index for Racon round ${i}."; exit 1; fi

    echo "  Running Racon (v1.5) round ${i}..."
    racon \
        -t ${THREADS} \
        ${FILTERED_FASTQ_SAMPLE} \
        ${ALIGNED_BAM_FOR_RACON} \
        ${POLISHED_FA_RAC_ROUNDS} \
        > ${CURRENT_POLISHED_FA}
    if [ $? -ne 0 ]; then echo "Error in Racon round <span class="math-inline">\{i\}\."; exit 1; fi
POLISHED\_FA\_RAC\_ROUNDS\=</span>{CURRENT_POLISHED_FA}
    rm ${ALIGNED_BAM_FOR_RACON} ${ALIGNED_BAM_FOR_RACON}.bai # Clean up BAM and index
done
echo "Racon polishing completed. Final Racon polished assembly: ${POLISHED_FA_RAC_ROUNDS}"
```

---

## 8. Assembly Polishing (Option A: Medaka on metaFlye+Racon assembly)

The metaFlye assemblies, after Racon polishing, were further polished with Medaka v2.0.1.

Tool: Medaka (v2.0.1) 
Input: Racon-polished metaFlye assembly and the corresponding filtered FASTQ file.
Output: Medaka-polished assembly directory.

```bash
# Activate Mamba environment for Medaka
# mamba activate medaka_env

# Define variables
FILTERED_FASTQ_SAMPLE="processed_reads/filtered/barcodeXX.filtered.fastq" # From Step 4
RACON_POLISHED_ASSEMBLY="assembly_output/barcodeXX_metaflye/racon_polished/racon_round3.fasta" # From Step 8
MEDAKA_POLISHED_DIR_METAFYLYE="assembly_output/barcodeXX_metaflye/medaka_polished"
THREADS=<N> # e.g., 16
MEDAKA_MODEL="<model_appropriate_for_r10.4.1_data>" # e.g., r1041_e82_400bps_sup_g615

# Create output directory (Medaka creates its own subdirectory structure within this)
mkdir -p ${MEDAKA_POLISHED_DIR_METAFYLYE}

echo "Starting Medaka (v2.0.1) polishing for metaFlye+Racon assembly..."
medaka_consensus \
    -i ${FILTERED_FASTQ_SAMPLE} \
    -d ${RACON_POLISHED_ASSEMBLY} \
    -o ${MEDAKA_POLISHED_DIR_METAFYLYE} \
    -t ${THREADS} \
    -m ${MEDAKA_MODEL}

if [ $? -ne 0 ]; then echo "Error in Medaka polishing for metaFlye assembly."; exit 1; fi
echo "Medaka polishing completed. Final metaFlye polished assembly: ${MEDAKA_POLISHED_DIR_METAFYLYE}/consensus.fasta"
```

---

## 9.Metagenome Assembly (Option B: nanoMDBG/metaMDBG)

De novo assemblies were also generated using nanoMDBG v1.1  (user noted this tool might be called metaMDBG, command often used with --in-ont for Nanopore).

Tool: nanoMDBG v1.1 / metaMDBG 
Input: A filtered FASTQ file (e.g., barcodeXX.filtered.fastq from Step 4).
Output: Raw nanoMDBG/metaMDBG assembly directory.

```bash
# Activate Mamba environment for nanoMDBG/metaMDBG
# mamba activate nanomdbg_env # or metamdbg_env

# Define variables
FILTERED_FASTQ_SAMPLE="processed_reads/filtered/barcodeXX.filtered.fastq" # From Step 4
ASSEMBLY_DIR_METAMDBG_RAW="assembly_output/barcodeXX_metamdbg/metamdbg_raw_assembly" # Changed from nanomdbg to reflect user's note
THREADS=<N> # e.g., 16

# Create parent output directory
mkdir -p assembly_output/barcodeXX_metamdbg

echo "Running metaMDBG (nanoMDBG v1.1) assembly for ${FILTERED_FASTQ_SAMPLE}..."
# Consult metaMDBG/nanoMDBG documentation for exact command structure for Nanopore reads.
# The PDF refers to "nanoMDBG v1.1". User noted tool might be "metaMDBG" and used with --in-ont.
# This is a placeholder command structure.
# Example metaMDBG command (structure is hypothetical, verify with tool's documentation):
# metamdbg assemble \
#    --in-ont ${FILTERED_FASTQ_SAMPLE} \
#    --out-dir ${ASSEMBLY_DIR_METAMDBG_RAW} \
#    --prefix barcodeXX \
#    --threads ${THREADS}
# Ensure the output contig file path matches INITIAL_ASSEMBLY_FASTA_METAMDBG in the next step.

echo "Placeholder: Run metaMDBG/nanoMDBG command here. Ensure output contigs are generated in ${ASSEMBLY_DIR_METAMDBG_RAW}"
# For example, if the output is ${ASSEMBLY_DIR_METAMDBG_RAW}/graph_contigs.fasta
# if [ $? -ne 0 ]; then echo "Error in metaMDBG/nanoMDBG step."; exit 1; fi
echo "metaMDBG/nanoMDBG assembly completed (using placeholder command)."
```

---


