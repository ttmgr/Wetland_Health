# Wetland Health Surveillance: Microbial and Viral Dynamics in Response to Anthropogenic Pressures

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) **Authors:** Tim Reska<sup>+</sup>, Albert Perlas Puente<sup>+</sup>, xxx, yyy, zzz, Lara Urban<sup>*</sup>

<sup>+</sup> contributed equally
<sup>*</sup> corresponding author

---

## Abstract

This project investigates the hypothesis that anthropogenic pressures on wetland environments alter their microbial and viral landscapes, potentially impacting wildlife health and increasing the environmental prevalence or diversity of pathogenic agents and antimicrobial resistance (AMR) determinants[cite: 35]. We aimed to develop and apply non-invasive surveillance methods coupled with nanopore technology to characterize the RNA virome (with a specific focus on Avian Influenza Virus - AIV) and to identify AMR genes and potential bacterial pathogens[cite: 36]. Our approach integrated AI-based land usage classification from environmental remote sensing data with the analysis of nucleic acids extracted from environmental water and air samples[cite: 37]. This methodology allows for the characterization of key elements of the microbial landscape within these vital ecosystems[cite: 39].

---

## Introduction

Wetlands are critical ecosystems supporting immense biodiversity and providing essential services such as climate regulation and water filtration[cite: 9, 10]. However, they face rapid degradation[cite: 11]. The health of their wildlife, particularly migratory birds, is crucial from a One Health perspective due to their potential to spread zoonotic diseases and AMR[cite: 13, 14]. This study focuses on understanding how varying anthropogenic stressors impact the microbial and viral dynamics within wetland environments, employing advanced molecular and computational techniques[cite: 34, 35].

---

## Research Goals

The primary goals of this study were to:
* Develop and apply non-invasive environmental surveillance methods (air and water sampling) for wetlands[cite: 36].
* Utilize nanopore sequencing technology to characterize the RNA virome, with a specific focus on Avian Influenza Viruses (AIV)[cite: 36].
* Identify and quantify antimicrobial resistance (AMR) genes and potential bacterial pathogens in collected samples[cite: 36].
* Assess how different degrees of anthropogenic impact influence the microbial and viral profiles of wetland ecosystems[cite: 35].
* Integrate AI-based land usage classification with molecular data for a comprehensive ecological assessment[cite: 37].

---

## Methodology Overview

### Study Sites and Sampling
* **Locations:** 12 distinct sites across Germany, France, and Spain, categorized by country and as anthropogenic ('H') or non-anthropogenic ('N')[cite: 40, 41].
    * Germany (Greifswald area): GH1 (landfill), GH2 (village edge), GH3 (city center), GN1 (meadow plain), GN2 (coastal marsh), GN3 (peat bog)[cite: 46, 47, 48, 49, 50].
    * France (Toulouse region): FH1 & FH2 (ponds near duck farms), FN1 (natural wetland area)[cite: 51, 52].
    * Spain (Ciudad Real province): SH1 (agricultural fertile plain), SN1 (natural area Peralbillo), SN2 (wetland periphery)[cite: 53, 54, 55].
* **Sampling Protocols:**
    * **Air:** Active sampling using Coriolis® μ liquid impinger[cite: 3, 60].
    * **Water:**
        * Active: 300 mL water samples, filtered ($0.45~\mu m$ pore size)[cite: 58, 59].
        * Passive: Torpedo-shaped samplers with filter discs, deployed for 3 days[cite: 63, 64, 65].

### Laboratory Processing & Sequencing
* **Nucleic Acid Extraction:** Co-extraction of DNA and RNA using Qiagen AllPrep PowerFecal Pro DNA/RNA Kit[cite: 66].
* **Quantification:** Qubit 4.0 fluorometer[cite: 76].
* **DNA Shotgun Sequencing:**
    * Library Prep: ONT Rapid Barcoding Kit (RBK114-24)[cite: 78].
    * Sequencing: MinION R10.4.1 flow cells, 24-hour runs[cite: 77, 78].
    * Basecalling: Dorado SUP mode[cite: 83].
* **RNA Processing & AIV Analysis:**
    * qPCR for influenza[cite: 5].
    * Targeted sequencing (SMART9n or other)[cite: 5].
    * AIV multi-segment M-RTPCR[cite: 104].
    * Sequencing: MinION Mk1C, R10.4.1 flow cell, Rapid Barcoding[cite: 108].
    * Basecalling: Dorado[cite: 109].

### Downstream Analysis
* **DNA Metagenomics:**
    * Preprocessing: Porechop for adapter/barcode removal, NanoFilt for read filtering[cite: 84, 85].
    * Taxonomic Classification: Kraken2 (nt_core database)[cite: 88].
    * Community Analysis: PCoA on downsampled reads[cite: 86, 90].
    * Assembly: metaFlye and nanoMDBG; polished with Minimap2, Racon, Medaka[cite: 91, 92].
    * Gene Annotation: Prokka[cite: 94].
    * AMR Gene Detection: AMRFinderPlus on downsampled reads and nanoMDBG contigs[cite: 95, 96, 98].
    * Taxonomic Origin of AMR on Contigs: DIAMOND (NCBI nr) and Kraken2 (NCBI nt_core)[cite: 100, 101].
* **AIV Analysis:**
    * Quality/Length Filtering: Filtlong[cite: 110].
    * Alignment: Minimap2 to reference genomes (NCBI Influenza Virus Database - Europe)[cite: 111, 112].
    * Consensus Generation: SAMtools, BCFtools[cite: 113, 116].

---

## Key Findings (Summary from Manuscript)
* DNA yields varied by sample type: active water (495.00–5400.00 ng), passive water (22.68–917.50 ng), and air (15.00–917.50 ng)[cite: 118].
* Nanopore shotgun sequencing generated 93k to 381k high-quality reads for active water, 28k to 577k for passive water, and 376 to 381k for air samples[cite: 119].
* PCoA revealed distinct microbial compositions for active water, passive water, and air samples, with the first two coordinates explaining 10.81% and 18.05% of the variance, respectively[cite: 122, 123, 124].
* Taxonomic profiles differed between natural and human-impacted sites. For example, human-impacted active water often featured *Flavobacterium* and *Limnohabitans*[cite: 127], while natural active waters showed *Woronichinia* or *Microcystis* at some sites[cite: 125]. Human-impacted passive water sites like FH1 (duck farm) were dominated by *Aeromonas*, and SH1 (landfill-adjacent) showed high levels of *Aeromonas* and *Vibrio*[cite: 129].
* Comparison of assemblers showed nanoMDBG generated more contigs (median: 1,965 vs. 254) and total assembled bases (median: 6.8 Mb vs. 2.2 Mb) than metaFlye[cite: 135]. MetaFlye produced more contiguous assemblies (median contig length: 5.0 kb vs. 2.6 kb)[cite: 136]. nanoMDBG was chosen for AMR analysis due to identifying more unique species (median: 223 vs. 93) and annotating more genes (median: 6,736 vs. 154)[cite: 137, 139].
* AMR gene profiles varied by location and environmental type. Passive water samples from human-impacted sites like FH1 (duck farm) showed high AMR gene counts (e.g., *blaOXA*, *cphA*) and diversity (thirteen distinct genes in reads)[cite: 141, 144, 145]. Landfill-adjacent site SH1 also showed notable diversity with eight distinct AMR genes in reads[cite: 146].
* Contig-level analysis confirmed human-impacted sites as key locations for assembled AMR genes. For example, at FH1, *arsB* and *arsC* were assigned to *Acidiphilium multivorum*, while *blaOXA* and *cphA* were assigned to *Aeromonas veronii*[cite: 153, 154, 155].

*(More detailed findings related to specific taxa, AMR genes, and AIV can be found in the full publication.)*

---

## Data Availability

All raw sequencing data has been made publicly available via the European Nucleotide Archive (ENA)[cite: 163].
* **Study Accession Number:** `PRJEBXXXXX` (Please replace `PRJEBXXXXX` with your actual ENA study accession number once available.)

---

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

*(Please update this section with the full citation once your manuscript is published.)*

---

## Acknowledgements

*(You can add any specific acknowledgements here if you wish, mirroring your manuscript or adding any specific to the GitHub repository maintenance).*

---

## License

The code in this repository is licensed under the MIT License. See the `LICENSE` file for more details.
*(You should add a LICENSE file to your repository. MIT is a common permissive license, but you can choose another like GPL, Apache 2.0, etc., depending on your preferences.)*
