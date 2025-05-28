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

*<Bash commands for Dorado basecalling and demux will be added here in the next step.>*

**Output:** Demultiplexed FASTQ files (one per barcode/sample) in the specified output directory. Each subsequent step will typically be run on these individual demultiplexed FASTQ files.

---

## 2. Read Processing (Adapter Trimming, Quality/Length Filtering)
This step cleans the demultiplexed FASTQ files by removing any remaining adapter sequences and filtering reads based on quality and length.
**Tools:** Porechop (v0.2.4)[cite: 84], NanoFilt (v2.8.0) [cite: 85]
**Input:** A single demultiplexed FASTQ file (e.g., `barcode01.fastq`) from the Dorado output.
**Output:** A processed FASTQ file (e.g., `barcode01.filtered.fastq`).

*<Bash commands for Porechop and NanoFilt will be added here in the next step.>*

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
