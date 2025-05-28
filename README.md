# Wetland Metagenomics & Viromics by Nanopore Sequencing (PRJEBXXXXX)

## Project Overview

This project utilizes metagenomic (DNA) and viromic (RNA) analysis of environmental water and air samples to characterize microbial communities, antimicrobial resistance (AMR) genes, and Avian Influenza Viruses (AIV), employing Oxford Nanopore Technologies (ONT) sequencing. This repository contains workflows designed to process raw sequencing data for these analyses.

## ENA Files: [https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX](https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX) ---

**❗ Important First Step ❗**

Before using the general pipelines described below, you **must** consult the specific data guides. These guides contain critical information about:

* **Data Access:** Where to find the raw (POD5) and/or processed (FASTQ) data.
* **Sample Mapping & Barcoding:** Which barcode corresponds to which sample for different sequencing runs (DNA shotgun[cite: 78], AIV [cite: 108]). Details on library prep (e.g., RBK114-24 for DNA[cite: 78], SQK-RBK114.24 for AIV [cite: 108]) and pooling strategies[cite: 79, 80, 81, 82].
* **Basecalling/Demultiplexing:** Specific Dorado commands and configurations used for the initial conversion of raw POD5 data to FASTQ files[cite: 83, 109], including required demultiplexing based on barcodes[cite: 83].

**Find the essential guides here:**

* **DNA Shotgun Data Guide:** [`dna_shotgun_data_guide.md`](./dna_shotgun_data_guide.md)
    * Describes samples for DNA metagenomic analysis (active water, passive water, air).
    * Uses **POD5** data format and **Dorado** (v5.0.0, model dna_r10.4.1_e8.2_400bps_sup@v5.0.0) for basecalling/demultiplexing[cite: 83].
* **AIV (RNA) Data Guide:** [`aiv_rna_data_guide.md`](./aiv_rna_data_guide.md)
    * Describes samples positive for AIV and processed for RNA sequencing.
    * Uses **POD5** data format and **Dorado** for basecalling/demultiplexing (primers and adapters removed by Dorado)[cite: 109].

**The pipelines described below assume you have already completed the necessary steps from the relevant guide and have demultiplexed FASTQ files ready for analysis.**

---

## Analysis Pipeline Overview

This repository outlines two main analysis pipelines:

1.  **DNA Shotgun Metagenomics:** Processes FASTQ files from environmental DNA samples.
    * Read Processing: Adapter trimming and quality/length filtering[cite: 84, 85].
    * Taxonomic Classification: Assigning taxonomy to reads[cite: 88].
    * Metagenome Assembly: Assembling reads into contigs using two different assemblers (metaFlye, nanoMDBG)[cite: 91, 92].
    * Assembly Polishing: Improving assembly accuracy using Racon and Medaka[cite: 91, 92].
    * AMR Gene Detection: Identifying antimicrobial resistance genes from reads and contigs[cite: 95, 98].
    * Taxonomic Origin of AMR: Determining the host of AMR genes on contigs[cite: 99, 101, 102].
2.  **AIV (RNA) Analysis:** Processes FASTQ files from samples amplified for AIV.
    * Read Processing: Quality and length filtering[cite: 110].
    * Alignment: Aligning reads to AIV reference genomes[cite: 111, 112].
    * Consensus Sequence Generation: Creating a consensus AIV genome for each sample[cite: 113, 114, 115, 116].

## Tools Used

This project integrates the following key bioinformatics tools:

* **Dorado:** Basecaller for ONT data (v5.0.0 for DNA[cite: 83], specified version for AIV [cite: 109]).
* **Porechop:** Adapter and barcode trimming (v0.2.4)[cite: 84].
* **NanoFilt:** Quality and length filtering for ONT reads (v2.8.0)[cite: 85].
* **Filtlong:** Quality and length filtering for AIV reads[cite: 110].
* **Kraken2:** K-mer based taxonomic classification (v2.1.2)[cite: 88].
* **metaFlye:** Long-read assembler for metagenomes (v2.9.6)[cite: 91].
* **nanoMDBG:** Long-read assembler for metagenomes (v1.1)[cite: 92].
* **Minimap2:** Long-read alignment (v2.28 for polishing DNA assemblies[cite: 91], v2.26 for AIV alignment [cite: 111]).
* **Racon:** Consensus correction/polishing for assemblies (v1.5)[cite: 91].
* **Medaka:** Consensus correction/polishing for assemblies (v2.0.1)[cite: 92].
* **AMRFinderPlus:** Detection of AMR genes (v3.12.8)[cite: 95].
* **DIAMOND:** Protein sequence alignment for taxonomic assignment of AMR-carrying contigs[cite: 101].
* **Seqkit:** Toolkit for FASTA/Q sequence manipulation (v2.10.0)[cite: 86, 96, 97].
* **SAMtools:** Utilities for SAM/BAM alignment files (v1.17)[cite: 113].
* **BCFtools:** Utilities for variant calling and consensus generation (v1.17)[cite: 116].
* **Python Libraries for PCoA:** scikit-bio v0.6.3, Matplotlib v3.10.0, Pandas v2.2.3, NumPy v1.26.4[cite: 90].

## Repository Structure

* `dna_shotgun_data_guide.md`: Essential guide for accessing raw data, sample mapping, barcoding, and basecalling/demultiplexing for DNA shotgun data.
* `aiv_rna_data_guide.md`: Essential guide for accessing raw data, sample mapping, barcoding, and basecalling/demultiplexing for AIV (RNA) data.
* `dna_shotgun_analysis_pipeline.md`: Detailed commands and scripts for the DNA shotgun metagenomics analysis workflow.
* `aiv_rna_analysis_pipeline.md`: Detailed commands and scripts for the AIV (RNA) analysis workflow.
* `Installation_tutorial.md`: Step-by-step guide for installing all required tools and databases.
* `/scripts` (example directory): May contain supplementary helper scripts.

## Installation

For comprehensive installation instructions for all pipeline tools and required databases, please refer to the [`Installation_tutorial.md`](./Installation_tutorial.md) file.

## Usage Workflow

This section outlines the main steps for processing samples. Detailed commands and scripts can be found in the linked pipeline documents.

### DNA Shotgun Metagenomics Workflow

1.  **Read Processing:** Sequencing adapters and barcodes are removed[cite: 84], and reads are filtered by length (minimum 100 bp)[cite: 85].
2.  **Taxonomic Classification:** Filtered reads are classified using Kraken2[cite: 88]. For Principal Coordinate Analysis (PCoA), reads are first downsampled (e.g., to 14,000 reads per sample)[cite: 86].
3.  **Metagenome Assembly & Polishing:** Assemblies are generated using metaFlye (v2.9.6) [cite: 91] (polished with Minimap2 v2.28 and three rounds of Racon v1.5 [cite: 91]) and nanoMDBG (v1.1)[cite: 92]. Assemblies from both are further polished with Medaka (v2.0.1)[cite: 92].
4.  **AMR Gene Detection:** AMRFinderPlus (v3.12.8) is used on reads (downsampled to a fixed threshold per sample type) [cite: 96, 97] and on nanoMDBG contigs[cite: 98].
5.  **Taxonomic Origin of AMR on Contigs:** Contigs carrying AMR genes are taxonomically classified using a dual approach with DIAMOND (NCBI nr database) [cite: 101] and Kraken2 (NCBI nt_core database)[cite: 101]. A species assignment requires agreement between both tools[cite: 102].

    *For detailed commands and scripts, see [`dna_shotgun_analysis_pipeline.md`](./dna_shotgun_analysis_pipeline.md).*

### AIV (RNA) Analysis Workflow

1.  **Read Processing:** Following basecalling and demultiplexing with Dorado (which also removed primers and adapters)[cite: 109], reads are filtered for quality (minimum Phred score >8) and length (>150 bp) using Filtlong[cite: 110].
2.  **Alignment to Reference Genomes:** Filtered reads are aligned to a comprehensive AIV reference database (European sequences from NCBI Influenza Virus Database as of 04/03/2023) using Minimap2 (v2.26) with the `-ax map-ont` setting[cite: 111, 112].
3.  **Consensus Sequence Generation:** Alignment files are processed using SAMtools (v1.17)[cite: 113]. The best reference for each of the eight AIV segments is determined by identifying the reference to which most reads mapped[cite: 114, 115]. Consensus sequences for each segment are then generated using BCFtools (v1.17)[cite: 116].

    *For detailed commands and scripts, see [`aiv_rna_analysis_pipeline.md`](./aiv_rna_analysis_pipeline.md).*

**Notes:**

* Replace placeholders in linked files (e.g., paths, sample names, database locations, Medaka models) with your actual information.
* Adjust thread counts based on your system's resources.
* This workflow provides a template based on the provided PDF[cite: 1]. Refer to individual tool documentation for detailed usage and optimization.
