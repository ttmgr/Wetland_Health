# Wetland Health Surveillance: Microbial and Viral Dynamics in Response to Anthropogenic Pressures

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Authors:** Tim Reska<sup>+</sup>, Albert Perlas Puente<sup>+</sup>, xxx, yyy, zzz, Lara Urban<sup>*</sup>

<sup>+</sup> contributed equally
<sup>*</sup> corresponding author

---

## Abstract

This project investigates the hypothesis that anthropogenic pressures on wetland environments alter their microbial and viral landscapes, potentially impacting wildlife health and increasing the environmental prevalence or diversity of pathogenic agents and antimicrobial resistance (AMR) determinants[cite: 35]. We aimed to develop and apply non-invasive surveillance methods coupled with nanopore technology to characterize the RNA virome (with a specific focus on Avian Influenza Virus - AIV) and to identify AMR genes and potential bacterial pathogens present in environmental samples from wetland ecosystems subject to varying degrees of anthropogenic impact[cite: 36]. Our approach integrated AI-based land usage classification from environmental remote sensing data with the analysis of nucleic acids extracted from non-invasive environmental samples (water and air)[cite: 37]. This methodology employed nanopore sequencing technology for data generation, enabling subsequent interpretation for viral, bacterial, and AMR profiling from the sampled environment[cite: 38].

---

## Introduction

Wetlands cover a significant portion of the Earth's surface and are indispensable for global biodiversity and providing critical ecosystem services[cite: 9, 10]. Despite their importance, wetlands face rapid degradation[cite: 11]. The health of wildlife populations in these areas, particularly migratory birds, is crucial from a One Health perspective, given their potential to spread zoonotic diseases and antimicrobial resistances[cite: 13, 14]. This study focuses on understanding how varying anthropogenic stressors impact the microbial and viral dynamics within wetland environments[cite: 34, 35].

---

## Research Goals

The primary goals of this study were to:
* Develop and apply non-invasive, surveillance methods and nanopore technology to characterize the RNA virome (with a specific focus on AIV)[cite: 36].
* Identify AMR genes and potential bacterial pathogens present in environmental samples from wetland ecosystems subject to varying degrees of anthropogenic impact[cite: 36].
* Assess how different degrees of anthropogenic impact influence the microbial and viral profiles of wetland ecosystems by integrating Al-based land usage classification from environmental remote sensing data with the analysis of nucleic acids extracted from non-invasive environmental samples (water and air)[cite: 37].

---

## Methodology Overview

### Study Sites and Sampling
* **Locations:** 12 distinct sites across Germany, France, and Spain, chosen to represent a range of environmental conditions and assigned unique identifiers[cite: 40, 41]. Sites were categorized by country (G, F, S) and as anthropogenic ('H') or non-anthropogenic ('N')[cite: 43, 44].
    * Germany (Greifswald area): GH1 ("Greifswald-landfill"), GH2 ("Greifswald Wackerow"), GH3 ("Greifswald City center"), GN1 ("Greifswald Karrendorfer Wiesen"), GN2 ("Greifswald Karrendorfer Wiesen Sea"), GN3 ("Greifswald Mannhagener Moor")[cite: 46, 47, 48, 49, 50].
    * France (Toulouse region): FH1 ("Toulouse Farm Dupoux Roucolle"), FH2 ("Toulouse Farm SCEA"), FN1 ("Toulouse Natural area")[cite: 51, 52].
    * Spain (Ciudad Real province): SH1 ("Ciudad Real Veguilla"), SN1 ("Ciudad Real Peralbillo"), SN2 ("Ciudad Real Daimiel")[cite: 53, 54, 55].
* **Sampling Protocols:**
    * **Air:** Active air sampling using the Coriolis® μ liquid impinger (Bertin Technologies, France)[cite: 60]. The sampler was operated at 1.5 m above ground, for 3-hour events, using 10 mL of 0.005% (v/v) Triton X in ultrapure water as collection liquid, at an air flow rate of 300 L/min[cite: 61, 62].
    * **Water:**
        * Active: Two 300 mL water samples collected per site, filtered through $0.45~\mu m$ pore size mixed cellulose ester filter discs[cite: 58, 59].
        * Passive: Torpedo-shaped sampling devices deployed (three per site in Germany, two per site in France/Spain), each with three $0.45~\mu m$ pore size mixed cellulose ester filter discs, retrieved three days post-installation[cite: 63, 64, 65].

### Laboratory Processing & Sequencing
* **Nucleic Acid Extraction:** Co-extraction of DNA and RNA from all samples using the AllPrep PowerFecal Pro DNA/RNA Kit (QIAGEN) with sample-specific pre-processing[cite: 66]. DNA was eluted in 25 µL and RNA in 50 µL of nuclease-free water[cite: 75].
* **Quantification:** DNA and RNA concentrations determined using a Qubit 4.0 fluorometer with Qubit dsDNA HS and RNA HS Assay Kits[cite: 76].
* **DNA Shotgun Sequencing:**
    * Library Prep: ONT Rapid Barcoding Kit 114-24 (RBK114-24)[cite: 78].
    * Sequencing: MinION device (Oxford Nanopore Technologies) using R10.4.1 MinION flow cells, 24-hour runs, 5kHz sampling frequency, minimum read length of 20 bases[cite: 77, 78].
    * Basecalling: Dorado v5.0.0 in super-accuracy (SUP) mode (model dna_r10.4.1_e8.2_400bps_sup@v5.0.0)[cite: 83].
* **RNA Processing & AIV Analysis:**
    * qPCR for influenza[cite: 1].
    * Targeted sequencing (SMART9n or other)[cite: 1].
    * AIV multi-segment M-RTPCR targeting conserved regions across all AIV segments[cite: 104]. Primers MBTuni-12 and MBTuni-13 were used[cite: 107].
    * Sequencing: MinION Mk1C device, FLO-MIN114 R10.4.1 flow cell, rapid barcoding library preparation (SQK-RBK114.24)[cite: 108].
    * Basecalling: Dorado, which also removed primers and adapters[cite: 109].

### Downstream Analysis
* **DNA Metagenomics:**
    * Preprocessing: Porechop v0.2.4 for adapter/barcode removal[cite: 84], NanoFilt v2.8.0 for read filtering (min length 100 bp)[cite: 85].
    * Taxonomic Classification: Kraken2 v2.1.2 (nt_core database)[cite: 88].
    * Community Analysis: PCoA on reads downsampled to 14,000 per sample (using Seqkit v2.10.0)[cite: 86], visualized using Python libraries[cite: 90].
    * Assembly: metaFlye v2.9.6 and nanoMDBG v1.1[cite: 91, 92]. Assemblies polished with Minimap2 v2.28, Racon v1.5 (three rounds for metaFlye), and Medaka v2.0.1 (for both)[cite: 91, 92].
    * Gene Annotation: Prokka v1.14.5[cite: 94].
    * AMR Gene Detection: AMRFinderPlus v3.12.8 on reads (downsampled per sample type) and nanoMDBG contigs[cite: 95, 96, 98].
    * Taxonomic Origin of AMR on Contigs: Dual approach using DIAMOND (NCBI nr database) and Kraken2 v2.1.2 (NCBI nt_core database); species assigned if both tools agreed[cite: 99, 101, 102].
* **AIV Analysis:**
    * Quality/Length Filtering: Filtlong (minimum Phred score >8, length >150 bp)[cite: 110].
    * Alignment: Minimap2 (v2.26) with -ax map-ont setting to a reference database from NCBI Influenza Virus Database (Europe sequences as of 04/03/2023)[cite: 111, 112].
    * Consensus Generation: SAMtools (v1.17) for BAM processing, selection of best reference per segment via samtools idxstats, and BCFtools (v1.17) for consensus sequence[cite: 113, 114, 115, 116].

---

## Data Availability

All raw sequencing data has been made publicly available via the European Nucleotide Archive (ENA).
* **Study Accession Number:** PRJEBXXXXX ---

## Code Availability

Relevant scripts for the core bioinformatic processing and analysis pipelines described in the manuscript are publicly available on this GitHub repository:
* **GitHub Repository:** [https://github.com/ttmgr/wetland_health](https://github.com/ttmgr/wetland_health) [cite: 164]

This includes scripts for:
* Downstream metagenomic analysis (e.g., basecalling, demultiplexing, filtering, taxonomic classification, assembly, AMR gene finding).
* AIV analysis pipelines.

Please refer to the specific directories within the repository for detailed instructions and dependencies. Scripts for final figure generation or specific statistical analyses presented in the manuscript are not included in this repository.

---

## How to Cite

If you use data or code from this project in your research, please cite our publication:

Reska, T., Perlas Puente, A., [xxx, yyy, zzz], & Urban, L. (Year). Title of your Manuscript. *Journal Name*. DOI/Link
---

## Acknowledgements

---

## License

The code in this repository is licensed under the MIT License. See the `LICENSE` file for more details.
