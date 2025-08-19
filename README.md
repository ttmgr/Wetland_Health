# Wetland Metagenomics, Viromics, and Vertebrate eDNA from European Wetlands (PRJEBXXXXX)

## Project Overview

This project investigates the microbial and viral communities, vertebrate biodiversity, and prevalence of antimicrobial resistance (AMR) genes in 12 European wetlands across Germany, France, and Spain. Samples were collected from sites categorized by land use (anthropogenic vs. natural) using torpedo-shaped passive water samplers.

From dual DNA/RNA extractions, four parallel analyses were conducted using Oxford Nanopore Technologies sequencing (MinION Mk1D/Mk1C, R10.4.1 flow cells):

1.  **Shotgun Metagenomics (eDNA):** Characterization of microbial communities and AMR genes.
2.  **Vertebrate Metabarcoding (eDNA):** Identification of local vertebrate fauna using 12S rRNA amplicons.
3.  **Avian Influenza Virus (AIV) Analysis (eRNA):** Targeted whole-genome sequencing of AIV-positive samples.
4.  **RNA Viromics (eRNA):** Untargeted analysis of RNA virus communities.

This repository contains the bioinformatic workflows used to process the sequencing data.

## ENA Files: [https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX](https://www.ebi.ac.uk/ena/browser/view/PRJEBXXXXX)

---

**❗ Bioinformatics - First Steps: Basecalling & QC ❗**

The pipelines below assume initial data processing has been completed. The following steps are universal or pipeline-specific first actions.

* **1. Basecalling:** Raw nanopore signal data (POD5 format) from all sequencing runs were basecalled using **Dorado v5.0.0** with the super-accuracy model (`dna_r10.4.1_e8.2_400bps_sup@v5.0.0`).

* **2. Demultiplexing & Trimming:** Barcodes and adapters were removed from the basecalled FASTQ files using pipeline-specific tools:
    * **Shotgun & Virome data:** `Porechop` (v0.2.4)
    * **12S Vertebrate data:** `OBITools4` (v1.3.1)
    * **AIV data:** Primers/adapters removed by `Dorado` during basecalling.

* **3. Initial Filtering:** Reads shorter than 100 bp were discarded using `NanoFilt` (v2.8.0) for shotgun and virome datasets. Further specific filtering is detailed in each pipeline.

**The pipelines described below assume you have completed these initial processing steps and have demultiplexed, trimmed FASTQ files ready for analysis.**

---

## Analysis Pipeline Overview

This repository outlines four main analysis pipelines:

1.  **DNA Shotgun Metagenomics:** Processes FASTQ files from environmental DNA for community analysis, assembly, and AMR/pathogen detection.
2.  **AIV (RNA) Analysis:** Aligns reads from AIV-positive samples to generate consensus genomes.
3.  **Viral Metagenomics (Viromics):** Processes cDNA reads from total environmental RNA to characterize the RNA virome.
4.  **12S Vertebrate Genomics:** Processes 12S rRNA gene amplicons to identify vertebrate species.

## Tools Used

This project integrates the following key bioinformatics tools:

* **Basecalling & QC:**
    * **Dorado:** Basecaller for ONT data (v5.0.0).
    * **Porechop:** Adapter and barcode trimming (v0.2.4).
    * **NanoFilt:** Quality and length filtering for ONT reads (v2.8.0).
    * **Seqkit:** Toolkit for FASTA/Q sequence manipulation (v2.3.0).
* **Metagenomics & Assembly:**
    * **Kraken2:** K-mer based taxonomic classification (v2.1.2).
    * **metaFlye:** Long-read assembler for metagenomes (v2.9.6).
    * **nanoMDBG:** Long-read assembler for metagenomes (v1.1).
    * **Minimap2:** Long-read alignment (v2.24 for polishing, v2.28 for AIV).
    * **Racon:** Consensus correction/polishing for assemblies (v1.5.0).
    * **Medaka:** Consensus correction/polishing for assemblies (v1.7.2).
    * **Prokka:** Prokaryotic genome annotation (v1.14.5).
* **Pathogen & AMR Analysis:**
    * **AMRFinderPlus:** Detection of AMR genes (v4.0.23).
    * **Prodigal:** ORF prediction (v2.6.3).
    * **DIAMOND:** Protein sequence alignment (v2.1.13).
    * **MEGAN-CE:** Interactive metagenomic analysis (v6.21.1).
    * **PlasmidFinder:** Plasmid detection from assemblies (v2.1.6).
* **Vertebrate & Virome Analysis:**
    * **OBITools4:** Barcode demultiplexing for amplicons (v1.3.1).
    * **Cutadapt:** Primer trimming (v4.2).
    * **VSEARCH:** Sequence analysis toolkit (clustering, chimera removal) (v2.21).
* **AIV & General Utilities:**
    * **SAMtools:** Utilities for SAM/BAM alignment files (v1.17).
    * **BCFtools:** Utilities for variant calling and consensus generation (v1.17).

## Repository Structure

* `dna_shotgun_analysis_pipeline.md`: Detailed commands for the DNA shotgun metagenomics workflow.
* `aiv_rna_analysis_pipeline.md`: Detailed commands for the AIV (RNA) analysis workflow.
* `virome_analysis_pipeline.md`: Detailed commands for the viral metagenomics workflow.
* `12s_vertebrate_analysis_pipeline.md`: Detailed commands for the 12S vertebrate analysis workflow.
* `Installation_tutorial.md`: Guide for installing all required tools and databases.

## Usage Workflow

### DNA Shotgun Metagenomics Workflow

1.  **Read Processing:** Reads were demultiplexed and trimmed with Porechop, then filtered to remove reads < 100 bp with NanoFilt.
2.  **Taxonomic Classification:** Reads were classified using Kraken2 against the NCBI nt_core database. A dual strategy was used:
    * **For community analysis (beta-diversity):** Samples were rarefied to 87,000 reads using Seqkit.
    * **For hazard detection (pathogens/AMR):** The complete, non-rarefied dataset was used to maximize sensitivity.
3.  **Metagenome Assembly & Annotation:** *De novo* assembly was performed with nanoMDBG, chosen over metaFlye because its lack of a minimum read length requirement prevented the >50% data loss that metaFlye would have caused. Assemblies were functionally annotated with Prokka.
4.  **Pathogen & AMR Gene Detection:**
    * **Pathogen ID:** Reads and contigs were aligned to nt_core with Minimap2 and assigned taxonomy with MEGAN-CE using a conservative Lowest Common Ancestor (LCA) approach.
    * **Virulence Factors:** For specific samples (e.g., high *Vibrio cholerae*), reads were screened against the Virulence Factor Database (VFDB) using DIAMOND BLASTx.
    * **AMR Genes:** AMRFinderPlus (in `--plus` mode with Prodigal for ORF prediction) was run on both reads and contigs.
    * **AMR Mobility:** Contigs were screened for plasmids using PlasmidFinder.

    *For detailed commands, see [`dna_shotgun_analysis_pipeline.md`](./dna_shotgun_analysis_pipeline.md).*

### AIV (RNA) Analysis Workflow

1.  **Read Processing:** After basecalling and demultiplexing with Dorado (which removes primers), reads were quality filtered.
2.  **Alignment to Reference Genomes:** Filtered reads were aligned to a custom European AIV reference database using Minimap2 (v2.28) with the `-ax map-ont` setting.
3.  **Consensus Sequence Generation:** SAMtools was used to process alignments and identify the best-matching reference for each of the eight AIV segments. A consensus sequence for each segment was then generated using BCFtools.

    *For detailed commands, see [`aiv_rna_analysis_pipeline.md`](./aiv_rna_analysis_pipeline.md).*

### Viral Metagenomics (Viromics) Workflow

1.  **Library Prep & Read Processing:** RNA was DNase-treated and then converted to cDNA using the **Rapid SMART-9N protocol**, which employs random priming and template switching. Barcoded amplicons were generated and sequenced. Reads were trimmed (Porechop) and filtered (NanoFilt).
2.  **Taxonomic Classification:** Reads were taxonomically classified by translated alignment using **DIAMOND BLASTx** against the NCBI non-redundant (NR) protein database. Viral hits were tallied to profile the RNA virome of each sample.

    *For detailed commands, see [`virome_analysis_pipeline.md`](./virome_analysis_pipeline.md).*

### 12S Vertebrate Genomics Workflow

1.  **Library Prep:** A ~97 bp fragment of the 12S rRNA gene was amplified using 12SV05 primers with 9 bp tags. A human-blocking oligonucleotide was included. The final library was prepared with the Ligation Sequencing Kit (SQK-LSK114).
2.  **Read Processing:** Reads were demultiplexed by their 9 bp tags using **OBITools4**. Primers were trimmed with **Cutadapt**.
3.  **OTU Clustering & Classification:** The VSEARCH pipeline was used to:
    * Filter reads by expected error (maxEE 1.0).
    * Dereplicate sequences and remove singletons.
    * Remove chimeras.
    * Cluster sequences into Operational Taxonomic Units (OTUs) at 97% similarity.
4.  **Taxonomic Assignment:** OTU representative sequences were identified via global alignment against the **MIDORI2** reference database.

    *For detailed commands, see [`12s_vertebrate_analysis_pipeline.md`](./12s_vertebrate_analysis_pipeline.md).*
