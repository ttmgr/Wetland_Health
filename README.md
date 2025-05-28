# Wetland Health Surveillance: Microbial and Viral Dynamics in Response to Anthropogenic Pressures

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Authors:** Tim Reska<sup>+</sup>, Albert Perlas Puente<sup>+</sup>, xxx, yyy, zzz, Lara Urban<sup>*</sup>

<sup>+</sup> contributed equally
<sup>*</sup> corresponding author

---

## Abstract

This project investigates the hypothesis that anthropogenic pressures on wetland environments alter their microbial and viral landscapes, potentially impacting wildlife health and increasing the environmental prevalence or diversity of pathogenic agents and antimicrobial resistance (AMR) determinants. We aimed to develop and apply non-invasive surveillance methods coupled with nanopore technology to characterize the RNA virome (with a specific focus on Avian Influenza Virus - AIV) and to identify AMR genes and potential bacterial pathogens present in environmental samples from wetland ecosystems subject to varying degrees of anthropogenic impact[cite: 35, 36]. Our approach integrated AI-based land usage classification from environmental remote sensing data with the analysis of nucleic acids extracted from non-invasive environmental samples (water and air)[cite: 37]. This methodology employed nanopore sequencing technology for data generation, enabling subsequent interpretation for viral, bacterial, and AMR profiling from the sampled environment[cite: 38].

---

## Introduction

Wetlands cover a significant portion of the Earth's surface (seven percent) and are indispensable for global biodiversity, supporting 40% of the world's biodiversity, and providing critical ecosystem services like climate regulation and water filtration[cite: 9, 10]. Despite their importance, wetlands face rapid degradation[cite: 11]. The health of wildlife populations in these areas, particularly migratory birds, is crucial from a One Health perspective, given their potential to spread zoonotic diseases and antimicrobial resistances[cite: 13, 14]. This study focuses on understanding how varying anthropogenic stressors impact the microbial and viral dynamics within wetland environments[cite: 34].

---

## Research Goals

The primary goals of this study were to:
* Develop and apply non-invasive, surveillance methods and nanopore technology to characterize the RNA virome (with a specific focus on AIV)[cite: 36].
* Identify AMR genes and potential bacterial pathogens present in environmental samples from wetland ecosystems subject to varying degrees of anthropogenic impact[cite: 36].
* Assess how different degrees of anthropogenic impact influence the microbial and viral profiles of wetland ecosystems by integrating Al-based land usage classification from environmental remote sensing data with the analysis of nucleic acids extracted from non-invasive environmental samples (water and air)[cite: 37, 38, 39].

---

## Methodology Overview

### Study Sites and Sampling
* **Locations:** 12 distinct sites across Germany, France, and Spain, chosen to represent a range of environmental conditions and assigned unique identifiers[cite: 40, 41]. Sites were categorized by country (G, F, S) and as anthropogenic ('H') or non-anthropogenic ('N')[cite: 42, 43, 44].
    * Germany (Greifswald area): GH1 ("Greifswald-landfill"), GH2 ("Greifswald Wackerow"), GH3 ("Greifswald City center"), GN1 ("Greifswald Karrendorfer Wiesen"), GN2 ("Greifswald Karrendorfer Wiesen Sea"), GN3 ("Greifswald Mannhagener Moor")[cite: 46, 47, 48, 49, 50].
    * France (Toulouse region): FH1 ("Toulouse Farm Dupoux Roucolle"), FH2 ("Toulouse Farm SCEA"), FN1 ("Toulouse Natural area")[cite: 51, 52].
    * Spain (Ciudad Real province): SH1 ("Ciudad Real Veguilla"), SN1 ("Ciudad Real Peralbillo"), SN2 ("Ciudad Real Daimiel")[cite: 53, 54, 55].
* **Sampling Protocols:** Environmental surveillance was conducted using both active and passive modalities for water and air samples[cite: 56].
    * **Air:** Active air sampling used the Coriolis® μ liquid impinger (Bertin Technologies, France)[cite: 60]. The sampler was operated at 1.5 m above ground, for 3-hour events, using 10 mL of 0.005% (v/v) Triton X in ultrapure water as collection liquid, at an air flow rate of 300 L/min[cite: 60, 61, 62].
    * **Water:**
        * Active: Two 300 mL water samples collected per site, filtered through $0.45~\mu m$ pore size mixed cellulose ester filter discs (Whatman® filter, Merck, USA)[cite: 58, 59].
        * Passive: Torpedo-shaped sampling devices deployed (three per site in Germany, two per site in France/Spain), each with three $0.45~\mu m$ pore size mixed cellulose ester filter discs, retrieved three days post-installation[cite: 63, 64, 65].

*(For details on laboratory processing, sequencing, and specific barcoding schemes, please refer to the linked documents in the "Workflow Details and Scripts" section below and the full manuscript.)*

### Downstream Analysis
* **DNA Metagenomics:**
    * Preprocessing: Porechop v0.2.4 for adapter/barcode removal, NanoFilt v2.8.0 for read filtering (min length 100 bp)[cite: 84, 85].
    * Taxonomic Classification: Kraken2 v2.1.2 (nt_core database, accessed May 2025)[cite: 88].
    * Community Analysis: PCoA on reads downsampled to 14,000 per sample (using Seqkit v2.10.0), visualized using Python libraries (Matplotlib v3.10.0, Pandas v2.2.3, NumPy v1.26.4, scikit-bio v0.6.3)[cite: 86, 89, 90].
    * Assembly: metaFlye v2.9.6 and nanoMDBG v1.1[cite: 91, 92]. Assemblies polished with Minimap2 v2.28, Racon v1.5 (three rounds for metaFlye), and Medaka v2.0.1 (for both)[cite: 91, 92].
    * Gene Annotation: Prokka v1.14.5[cite: 94].
    * AMR Gene Detection: AMRFinderPlus v3.12.8 on reads (downsampled per sample type: 87k for Passive Water, 93k for Active Water, 14k for Air) and nanoMDBG contigs[cite: 95, 96, 98].
    * Taxonomic Origin of AMR on Contigs: Dual approach using DIAMOND (NCBI nr database, accessed May 2025) and Kraken2 v2.1.2 (NCBI nt_core database, accessed May 2025); species assigned if both tools agreed[cite: 100, 101, 102].
* **AIV Analysis:**
    * AIV sequencing followed cDNA conversion and multi-segment amplification using M-RTPCR[cite: 104]. Primers MBTuni-12 and MBTuni-13 were used[cite: 107].
    * Raw POD5 files were basecalled using Dorado, which also removed primers and adapters[cite: 109].
    * Quality/Length Filtering: Filtlong (minimum Phred score >8, length >150 bp)[cite: 110].
    * Alignment: Minimap2 (v2.26) with -ax map-ont setting to a reference database generated for each segment from the NCBI Influenza Virus Database (Europe sequences as of 04/03/2023)[cite: 111, 112].
    * Consensus Generation: SAMtools (v1.17) for BAM processing, selection of best reference per segment via samtools idxstats, and BCFtools (v1.17) for consensus sequence[cite: 113, 114, 116].

---

## Workflow Details and Scripts

For more detailed information on specific parts of our workflow, please see the following documents in this repository:

* **Barcoding and Sample Scheme:** [`barcoding_sample_scheme.md`](./barcoding_sample_scheme.md)
    * This document outlines the DNA library preparation, barcoding strategy for different sample types (active water, passive water, air), and how samples were pooled for MinION sequencing runs[cite: 78, 79, 80, 81, 82].
* **Tool Installation Guide:** [`tool_installation.md`](./tool_installation.md)
    * Instructions and recommendations for installing the bioinformatic tools used in this project (e.g., Dorado, Porechop, NanoFilt, Kraken2, metaFlye, nanoMDBG, AMRFinderPlus, Minimap2, SAMtools, BCFtools, etc.).
* **Analysis Scripts and Commands:** [`analysis_pipelines.md`](./analysis_pipelines.md)
    * This file contains the specific bash commands and script snippets used for the DNA metagenomic and AIV analysis pipelines, from raw data processing to final output generation.

---

## Data Availability

All raw sequencing data has been made publicly available via the European Nucleotide Archive (ENA)[cite: 163].
* **Study Accession Number:** PRJEBXXXXX ---

## Code Availability

Core bioinformatic processing scripts and pipeline workflows as described in `analysis_pipelines.md` are available in this GitHub repository[cite: 164]:
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
