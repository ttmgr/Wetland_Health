# Wetland Metagenomic and Viromic Data Analysis Workflow

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Authors:** Tim Reska<sup>+</sup>, Albert Perlas Puente<sup>+</sup>, xxx, yyy, zzz, Lara Urban<sup>*</sup>

<sup>+</sup> contributed equally
<sup>*</sup> corresponding author

---

## Analysis Overview

This repository provides the bioinformatic workflows used for the analysis of DNA (shotgun metagenomics) and RNA (viromics, with a focus on Avian Influenza Virus - AIV) sequence data from environmental wetland samples. The primary analyses include taxonomic classification, de novo assembly, antimicrobial resistance (AMR) gene identification, and AIV characterization.

### DNA Metagenomic Analysis
* **Preprocessing:** Sequencing adapters and barcodes were removed using Porechop v0.2.4[cite: 84]. Reads were then filtered using NanoFilt v2.8.0, retaining only those with a minimum length of 100 base pairs[cite: 85].
* **Taxonomic Classification:** Kraken2 v2.1.2 (nt_core database, accessed May 2025) was used for taxonomic classification of reads[cite: 88].
* **Community Analysis:** For Principal Coordinate Analysis (PCoA), filtered reads were randomly downsampled to 14,000 reads per sample using Seqkit v2.10.0[cite: 86]. The PCoA was performed using scikit-bio v0.6.3 and visualized with Python libraries (Matplotlib v3.10.0, Pandas v2.2.3, NumPy v1.26.4)[cite: 90].
* **Assembly:** De novo assemblies were generated using metaFlye v2.9.6 (polished with Minimap2 v2.28 and three rounds of Racon v1.5) [cite: 91] and nanoMDBG v1.1[cite: 92]. Both assembler outputs were further polished with Medaka v2.0.1[cite: 92].
* **Gene Annotation:** Prokka v1.14.5 was used for gene annotation on assembled contigs[cite: 94].
* **AMR Gene Detection:** AMRFinderPlus v3.12.8 was applied to reads (downsampled per sample type: 87k for Passive Water, 93k for Active Water, 14k for Air [cite: 96]) and to nanoMDBG contigs[cite: 95, 98].
* **Taxonomic Origin of AMR on Contigs:** A dual approach using DIAMOND (NCBI nr database, accessed May 2025) and Kraken2 v2.1.2 (NCBI nt_core database, accessed May 2025) was employed. An AMR gene was assigned to a species if both tools provided the same species-level classification for the contig[cite: 101, 102].

### AIV (RNA Virus) Analysis
* **Sample Processing:** AIV sequencing of positive samples followed cDNA conversion and multi-segment amplification using M-RTPCR targeting conserved regions across all AIV segments[cite: 104]. Primers MBTuni-12 (5'-ACGCGTGATCAGCAAAAGCAGG) and MBTuni-13 (5'-ACGCGTGATCAGTAGAAACAAGG) were used[cite: 107].
* **Sequencing & Basecalling:** Sequencing was performed on a MinION Mk1C using R10.4.1 chemistry and rapid barcoding (SQK-RBK114.24)[cite: 108]. Raw POD5 files were basecalled using Dorado, which also removed primers and adapters[cite: 109].
* **Quality/Length Filtering:** Reads were filtered using Filtlong (minimum Phred score >8, length >150 bp)[cite: 110].
* **Alignment:** Filtered FASTQ files were aligned to reference genomes using Minimap2 (v2.26) with the -ax map-ont setting[cite: 111]. The reference database was generated for each segment from the NCBI Influenza Virus Database (all AIV nucleotide sequences from Europe as of 04/03/2023)[cite: 112].
* **Consensus Generation:** SAM files were converted to BAM, sorted, and indexed using SAMtools (v1.17)[cite: 113]. The segment reference to which most reads mapped was selected using `samtools idxstats`[cite: 114]. Consensus sequences were obtained with BCFtools (v1.17)[cite: 116].

---

## Workflow Details and Scripts

For more detailed information on specific parts of our workflow, please see the following documents in this repository:

* **Barcoding and Sample Scheme for Sequencing:** [`barcoding_sample_scheme.md`](./barcoding_sample_scheme.md)
    * This document outlines the DNA library preparation using the Rapid Barcoding Kit 114-24 (RBK114-24), barcoding strategy for different sample types (active water, passive water, air), and how samples were pooled for MinION R10.4.1 sequencing runs[cite: 78, 79, 80, 81, 82]. It also covers AIV library preparation (SQK-RBK114.24) for R10.4.1 chemistry[cite: 108].
* **Tool Installation Guide:** [`tool_installation.md`](./tool_installation.md)
    * Instructions and recommendations for installing the bioinformatic tools used in this project (e.g., Dorado v5.0.0[cite: 83], Porechop v0.2.4[cite: 84], NanoFilt v2.8.0[cite: 85], Kraken2 v2.1.2[cite: 88], metaFlye v2.9.6[cite: 91], nanoMDBG v1.1[cite: 92], AMRFinderPlus v3.12.8[cite: 95], Minimap2 v2.26/v2.28[cite: 91, 111], SAMtools v1.17[cite: 113], BCFtools v1.17[cite: 116], etc.).
* **Analysis Scripts and Commands:** [`analysis_pipelines.md`](./analysis_pipelines.md)
    * This file contains the specific bash commands and script snippets used for the DNA metagenomic and AIV analysis pipelines, from raw data processing to final output generation.

---

## Data Availability

All raw sequencing data has been made publicly available via the European Nucleotide Archive (ENA).
* **Study Accession Number:** PRJEBXXXXX ---

## Code Availability

Core bioinformatic processing scripts and pipeline workflows as described in `analysis_pipelines.md` are available in this GitHub repository:
* **GitHub Repository:** [https://github.com/ttmgr/wetland_health](https://github.com/ttmgr/wetland_health) [cite: 164]

Please refer to the specific directories within the repository and the linked markdown files for detailed instructions and dependencies. Scripts for final figure generation or specific statistical analyses presented in the manuscript are not included in this repository.

---

## How to Cite

If you use data or code from this project in your research, please cite our publication:

Reska, T., Perlas Puente, A., [xxx, yyy, zzz], & Urban, L. (Year). Title of your Manuscript. *Journal Name*. DOI/Link
---

## Acknowledgements

---

## License

The code in this repository is licensed under the MIT License. See the `LICENSE` file for more details.
