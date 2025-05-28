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

## 2. Read Processing (Adapter Trimming, Quality/Length Filtering)
This step cleans the demultiplexed FASTQ files by removing any remaining adapter sequences and filtering reads based on quality and length.
**Tools:** Porechop (v0.2.4)[cite: 84], NanoFilt (v2.8.0) [cite: 85]
**Input:** A single demultiplexed FASTQ file (e.g., `barcode01.fastq`) from the Dorado output.
**Output:** A processed FASTQ file (e.g., `barcode01.filtered.fastq`).

```bash
# Activate Mamba environments
# mamba activate porechop_env
# mamba activate nanofilt_env # Or run in sequence if preferred

# Define variables for a single sample
INPUT_DEMUX_FASTQ="demultiplexed_fastq/barcodeXX.fastq" # Replace barcodeXX with actual file
TRIMMED_FASTQ="processed_reads/barcodeXX.trimmed.fastq"
FILTERED_FASTQ="processed_reads/barcodeXX.filtered.fastq"
THREADS=<N> # Number of threads to use, e.g., 4

# Create output directory
mkdir -p processed_reads

echo "Processing sample: ${INPUT_DEMUX_FASTQ}"

# Step 2.1: Adapter and Barcode Trimming with Porechop
# Sequencing adapters and barcodes were removed using Porechop v0.2.4[cite: 84].
echo "Running Porechop for adapter trimming..."
porechop \
    -i ${INPUT_DEMUX_FASTQ} \
    -o ${TRIMMED_FASTQ} \
    --threads ${THREADS}

if [ $? -ne 0 ]; then
    echo "Error in Porechop step for ${INPUT_DEMUX_FASTQ}."
    # exit 1 # Decide if you want to exit or continue with other samples
fi
echo "Porechop completed for ${INPUT_DEMUX_FASTQ}."
```

# Step 2.2: Quality and Length Filtering with NanoFilt
# Reads were filtered using NanoFilt v2.8.0, keeping only reads with a minimum length of 100 base pairs[cite: 85].
echo "Running NanoFilt for quality and length filtering..."
NanoFilt \
    --length 100 \
    < ${TRIMMED_FASTQ} \
    > ${FILTERED_FASTQ}
    # Quality threshold not specified in PDF for NanoFilt for this step,
    # but a common default or Q score like 8 might be used if desired.
    # For example, adding --quality 8 if desired:
    # NanoFilt --length 100 --quality 8 < ${TRIMMED_FASTQ} > ${FILTERED_FASTQ}


if [ $? -ne 0 ]; then
    echo "Error in NanoFilt step for ${TRIMMED_FASTQ}."
    # exit 1
fi
echo "NanoFilt completed. Processed file: ${FILTERED_FASTQ}"

# Deactivate environments if you activated them for single tool usage
# mamba deactivate
```

*Repeat these steps for each demultiplexed FASTQ file.*

---

## 3. Taxonomic Classification (Reads)

Assign taxonomy to the processed reads. For Principal Coordinate Analysis (PCoA), reads were downsampled to 14,000 reads per sample. [cite: 86]

**Tools:** Kraken2 (v2.1.2)[cite: 88], Seqkit (v2.10.0 for downsampling) [cite: 86]

**Input:** A processed FASTQ file (e.g., `barcodeXX.filtered.fastq`).
**Output:** Kraken2 output and report files.

*<Bash commands for Seqkit (downsampling) and Kraken2 will be added here in the next step.>*

*Repeat for each sample. The PCoA itself using Python libraries (scikit-bio v0.6.3[cite: 90], Matplotlib v3.10.0[cite: 90], Pandas v2.2.3[cite: 90], NumPy v1.26.4 [cite: 90]) would be a separate script using these Kraken2 outputs.*

---

## 4. Metagenome Assembly and Polishing

Assemble reads into contigs and then polish these assemblies to improve accuracy. De novo assemblies were generated using metaFlye v2.9.6 [cite: 91] and nanoMDBG v1.1[cite: 92]. metaFlye assemblies were polished using Minimap2 v2.28 [cite: 91] and three rounds of Racon v1.5[cite: 91]. Assemblies from metaFlye and nanoMDBG were then polished with Medaka v2.0.1. [cite: 92]

**Tools:** metaFlye (v2.9.6)[cite: 91], nanoMDBG (v1.1)[cite: 92], Minimap2 (v2.28)[cite: 91], Racon (v1.5)[cite: 91], Medaka (v2.0.1) [cite: 92]

**Input:** A processed FASTQ file (e.g., `barcodeXX.filtered.fastq`).
**Output:** Polished assembly FASTA file(s).

### 4.1 Assembly Option A: metaFlye + Racon + Medaka

*<Bash commands for metaFlye, Minimap2, Racon, and Medaka (for metaFlye assembly) will be added here in the next step.>*

### 4.2 Assembly Option B: nanoMDBG + Medaka

(nanoMDBG assemblies were used for further resistance gene detection and species identification [cite: 139])
*<Bash commands for nanoMDBG and Medaka (for nanoMDBG assembly) will be added here in the next step.>*

*Repeat for each sample. The nanoMDBG polished assembly is used for downstream AMR analysis on contigs.*

---

## 5. Antimicrobial Resistance (AMR) Gene Detection

Detect AMR genes from both processed reads (downsampled) and assembled contigs (from nanoMDBG+Medaka). AMRFinderPlus v3.12.8 was used. [cite: 95]

**Tools:** AMRFinderPlus (v3.12.8)[cite: 95], Seqkit (v2.10.0 for downsampling reads) [cite: 96, 97]

**Input:**
* Processed FASTQ file (e.g., `barcodeXX.filtered.fastq`).
* Polished nanoMDBG assembly (e.g., `barcodeXX_nanomdbg/medaka_polished/consensus.fasta`).
**Output:** AMRFinderPlus report files for reads and contigs.

*<Bash commands for Seqkit (downsampling for AMR) and AMRFinderPlus will be added here in the next step.>*

*Repeat for each sample.*

---

## 6. Taxonomic Origin of AMR on Contigs

For AMR genes identified on contigs from nanoMDBG assemblies, determine their likely host species using a dual-approach involving DIAMOND (NCBI nr database, accessed May 2025) [cite: 101] and Kraken2 v2.1.2 (NCBI nt_core database, accessed May 2025)[cite: 101]. An AMR gene was assigned to a specific species only if both DIAMOND and Kraken2 showed the same species-level classification for that contig. [cite: 102] Mismatched classifications resulted in exclusion. [cite: 103]

**Tools:** DIAMOND[cite: 101], Kraken2 (v2.1.2) [cite: 101]

**Input:**
* Contigs identified by AMRFinderPlus as carrying AMR genes (from `amrfinder_contigs_output.txt`). You'll need to extract these contig sequences into individual FASTA files or a multi-FASTA file.
**Output:** Taxonomic assignment for each AMR-carrying contig.

*<Bash commands for DIAMOND and Kraken2 (for AMR contigs) will be added here in the next step.>*

*This section would typically be implemented with a script to automate the extraction of AMR-positive contigs and iterate the DIAMOND/Kraken2 analysis over them, followed by parsing and comparison.*

---

This pipeline provides a comprehensive workflow for analyzing your DNA shotgun metagenomic data. Remember to adapt paths, resource allocations (threads), and specific parameters to your computational environment and research questions.
