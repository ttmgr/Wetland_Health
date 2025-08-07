# Wetland Metagenomics, Viromics, and Vertebrate DNA by Nanopore Sequencing (PRJEBXXXXX)

## Project Overview

This project utilizes metagenomic (DNA), viral metagenomic (viromic), vertebrate 12S rRNA gene, and Avian Influenza Virus (AIV) RNA analysis of environmental water and air samples. The goal is to characterize microbial and viral communities, identify vertebrate species, and detect antimicrobial resistance (AMR) genes and AIV, employing Oxford Nanopore Technologies (ONT) sequencing. This repository contains workflows designed to process raw sequencing data for these analyses.

## ENA Files: [https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX](https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX) ---

**❗ Important First Step ❗**

Before using the general pipelines described below, you **must** consult the specific data guides. These guides contain critical information about:

* **Data Access:** Where to find the raw (POD5) and/or processed (FASTQ) data.
* **Sample Mapping & Barcoding:** Which barcode corresponds to which sample for different sequencing runs. Details on library prep (e.g., RBK114-24 for DNA, SQK-RBK114.24 for AIV) and pooling strategies.
* **Basecalling/Demultiplexing:** Specific Dorado commands and configurations used for the initial conversion of raw POD5 data to FASTQ files, including required demultiplexing based on barcodes.

**Find the essential guides here:**

* **DNA Shotgun Data Guide:** [`dna_shotgun_data_guide.md`](./dna_shotgun_data_guide.md)
    * Describes samples for DNA metagenomic analysis (active water, passive water, air).
    * Uses **POD5** data format and **Dorado** (v5.0.0, model dna_r10.4.1_e8.2_400bps_sup@v5.0.0) for basecalling/demultiplexing.
* **AIV (RNA) Data Guide:** [`aiv_rna_data_guide.md`](./aiv_rna_data_guide.md)
    * Describes samples positive for AIV and processed for RNA sequencing.
    * Uses **POD5** data format and **Dorado** for basecalling/demultiplexing (primers and adapters removed by Dorado).

**The pipelines described below assume you have already completed the necessary steps from the relevant guide and have demultiplexed FASTQ files ready for analysis.**

---

## Analysis Pipeline Overview

This repository outlines four main analysis pipelines:

1.  **DNA Shotgun Metagenomics:** Processes FASTQ files from environmental DNA samples.
    * Read Processing: Adapter trimming and quality/length filtering.
    * Taxonomic Classification: Assigning taxonomy to reads.
    * Metagenome Assembly: Assembling reads into contigs using two different assemblers (metaFlye, nanoMDBG).
    * Assembly Polishing: Improving assembly accuracy using Racon and Medaka.
    * AMR Gene Detection: Identifying antimicrobial resistance genes from reads and contigs.
    * Taxonomic Origin of AMR: Determining the host of AMR genes on contigs.
2.  **AIV (RNA) Analysis:** Processes FASTQ files from samples amplified for AIV.
    * Read Processing: Quality and length filtering.
    * Alignment: Aligning reads to AIV reference genomes.
    * Consensus Sequence Generation: Creating a consensus AIV genome for each sample.
3.  **Viral Metagenomics (Viromics):** Processes FASTQ files from total RNA converted to cDNA.
    * Read Processing: Adapter trimming and quality filtering.
    * Metagenome Assembly: Assembling reads into contigs.
    * Viral Contig Identification: Identifying and classifying viral sequences from the assembly.
4.  **12S Vertebrate Genomics:** Processes FASTQ files from 12S rRNA gene amplicons.
    * Read Processing: Primer trimming and quality/length filtering.
    * Clustering: Grouping similar sequences into Operational Taxonomic Units (OTUs).
    * Taxonomic Classification: Assigning taxonomy to identify vertebrate species.

## Tools Used

This project integrates the following key bioinformatics tools:

* **AMRFinderPlus:** Detection of AMR genes (v3.12.8).
* **BCFtools:** Utilities for variant calling and consensus generation (v1.17).
* **DIAMOND:** Protein sequence alignment for taxonomic assignment of AMR-carrying contigs.
* **Dorado:** Basecaller for ONT data (v5.0.0 for DNA, specified version for AIV).
* **Filtlong:** Quality and length filtering for AIV reads.
* **Kraken2:** K-mer based taxonomic classification (v2.1.2).
* **Medaka:** Consensus correction/polishing for assemblies (v2.0.1).
* **metaFlye:** Long-read assembler for metagenomes (v2.9.6).
* **Minimap2:** Long-read alignment (v2.28 for polishing DNA assemblies, v2.26 for AIV alignment).
* **MMseqs2:** High-sensitivity sequence searching for 12S classification.
* **NanoFilt:** Quality and length filtering for ONT reads (v2.8.0).
* **nanoMDBG:** Long-read assembler for metagenomes (v1.1).
* **Porechop:** Adapter and barcode trimming (v0.2.4).
* **Python Libraries for PCoA:** scikit-bio v0.6.3, Matplotlib v3.10.0, Pandas v2.2.3, NumPy v1.26.4.
* **Racon:** Consensus correction/polishing for assemblies (v1.5).
* **SAMtools:** Utilities for SAM/BAM alignment files (v1.17).
* **Seqkit:** Toolkit for FASTA/Q sequence manipulation (v2.10.0).
* **VirSorter2:** Identification of viral contigs from assemblies.

## Repository Structure

* `dna_shotgun_data_guide.md`: Essential guide for accessing raw data, sample mapping, barcoding, and basecalling/demultiplexing for DNA shotgun data.
* `aiv_rna_data_guide.md`: Essential guide for accessing raw data, sample mapping, barcoding, and basecalling/demultiplexing for AIV (RNA) data.
* `dna_shotgun_analysis_pipeline.md`: Detailed commands and scripts for the DNA shotgun metagenomics analysis workflow.
* `aiv_rna_analysis_pipeline.md`: Detailed commands and scripts for the AIV (RNA) analysis workflow.
* `virome_analysis_pipeline.md`: Detailed commands and scripts for the viral metagenomics workflow.
* `12s_vertebrate_analysis_pipeline.md`: Detailed commands and scripts for the 12S vertebrate analysis workflow.
* `Installation_tutorial.md`: Step-by-step guide for installing all required tools and databases.
* `/scripts` (example directory): May contain supplementary helper scripts.

## Installation

For comprehensive installation instructions for all pipeline tools and required databases, please refer to the [`Installation_tutorial.md`](./Installation_tutorial.md) file.

## Usage Workflow

This section outlines the main steps for processing samples. Detailed commands and scripts can be found in the linked pipeline documents.

### DNA Shotgun Metagenomics Workflow

1.  **Read Processing:** Sequencing adapters and barcodes are removed, and reads are filtered by length (minimum 100 bp).
2.  **Taxonomic Classification:** Filtered reads are classified using Kraken2. For Principal Coordinate Analysis (PCoA), reads are first downsampled (e.g., to 14,000 reads per sample).
3.  **Metagenome Assembly & Polishing:** Assemblies are generated using metaFlye (v2.9.6) (polished with Minimap2 v2.28 and three rounds of Racon v1.5) and nanoMDBG (v1.1). Assemblies from both are further polished with Medaka (v2.0.1).
4.  **AMR Gene Detection:** AMRFinderPlus (v3.12.8) is used on reads (downsampled to a fixed threshold per sample type) and on nanoMDBG contigs.
5.  **Taxonomic Origin of AMR on Contigs:** Contigs carrying AMR genes are taxonomically classified using a dual approach with DIAMOND (NCBI nr database) and Kraken2 (NCBI nt_core database). A species assignment requires agreement between both tools.

    *For detailed commands and scripts, see [`dna_shotgun_analysis_pipeline.md`](./dna_shotgun_analysis_pipeline.md).*

### AIV (RNA) Analysis Workflow

1.  **Read Processing:** Following basecalling and demultiplexing with Dorado (which also removed primers and adapters), reads are filtered for quality (minimum Phred score >8) and length (>150 bp) using Filtlong.
2.  **Alignment to Reference Genomes:** Filtered reads are aligned to a comprehensive AIV reference database (European sequences from NCBI Influenza Virus Database as of 04/03/2023) using Minimap2 (v2.26) with the `-ax map-ont` setting.
3.  **Consensus Sequence Generation:** Alignment files are processed using SAMtools (v1.17). The best reference for each of the eight AIV segments is determined by identifying the reference to which most reads mapped. Consensus sequences for each segment are then generated using BCFtools (v1.17).

    *For detailed commands and scripts, see [`aiv_rna_analysis_pipeline.md`](./aiv_rna_analysis_pipeline.md).*

### Viral Metagenomics (Viromics) Workflow

1.  **Library Prep & Read Processing:** Total RNA is reverse-transcribed to cDNA using SMART-9N primers. After sequencing, reads are trimmed for adapters and filtered for quality.
2.  **Metagenome Assembly:** Quality-filtered reads are assembled into contigs using an assembler like metaFlye.
3.  **Viral Contig Identification:** Assembled contigs are processed with VirSorter2 to identify and extract viral sequences.
4.  **Taxonomic Classification:** The identified viral contigs are then taxonomically classified to determine the composition of the viral community.

    *For detailed commands and scripts, see [`virome_analysis_pipeline.md`](./virome_analysis_pipeline.md).*

### 12S Vertebrate Genomics Workflow

1.  **Read Processing:** Reads are filtered for quality. Primers specific to the 12S rRNA gene are trimmed from both ends of the reads.
2.  **Dereplication & Clustering:** Reads are dereplicated (grouping identical sequences) and then clustered into Operational Taxonomic Units (OTUs) based on a similarity threshold (e.g., 98%).
3.  **Taxonomic Classification:** Representative sequences for each OTU are aligned against a curated 12S vertebrate reference database using MMseqs2 or BLAST to assign a taxonomic identity (e.g., species, family).

    *For detailed commands and scripts, see [`12s_vertebrate_analysis_pipeline.md`](./12s_vertebrate_analysis_pipeline.md).*

**Notes:**

* Replace placeholders in linked files (e.g., paths, sample names, database locations, Medaka models) with your actual information.
* Adjust thread counts based on your system's resources.
* This workflow provides a template based on the provided PDF. Refer to individual tool documentation for detailed usage and optimization.
