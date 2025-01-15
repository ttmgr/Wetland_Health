## Overview
This repository contains a collection of scripts for processing and analyzing metagenomic data from environmental samples, with a focus on Oxford Nanopore sequencing data.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Components](#pipeline-components)
- [Directory Structure](#directory-structure)
- [Scripts](#scripts)

## Prerequisites

### Required Software
- Dorado (for basecalling)
- Porechop (for adapter trimming)
- NanoFilt (for read filtering)
- NanoStat (for read metrics)
- Kraken2 (for taxonomic classification)
- Flye (for assembly)
- Minimap2 (for read mapping)
- Racon (for assembly polishing)
- assembly-stats (for assembly metrics)
- AMRFinder Plus (for antimicrobial resistance gene detection)
- seqtk (for FASTQ to FASTA conversion)
- GNU Parallel (for parallel processing)

### Required Databases
- Kraken2 database

## Installation

### Using Conda
```bash
# Create conda environment
conda create -n metagenome_pipeline python=3.9
conda activate metagenome_pipeline

# Install required tools
conda install -c bioconda -c conda-forge dorado
conda install -c bioconda -c conda-forge porechop
conda install -c bioconda -c conda-forge nanofilt
conda install -c bioconda -c conda-forge nanostat
conda install -c bioconda -c conda-forge kraken2
conda install -c bioconda -c conda-forge flye
conda install -c bioconda -c conda-forge minimap2
conda install -c bioconda -c conda-forge racon
conda install -c bioconda -c conda-forge assembly-stats
conda install -c bioconda -c conda-forge ncbi-amrfinderplus
conda install -c bioconda -c conda-forge seqtk
```

### Database Setup
```bash
# Download and setup Kraken2 database
kraken2-build --standard --db kraken2_db
```

## Usage

### Basic Usage
1. Clone this repository
```bash
git clone https://github.com/yourusername/metagenomic-pipeline.git
cd metagenomic-pipeline
```

2. Make scripts executable
```bash
chmod +x scripts/*.sh
```

3. Run the complete pipeline
```bash
./run_pipeline.sh -i input_dir -o output_dir -t threads
```

### Individual Steps

1. Basecalling
```bash
./scripts/01_basecalling.sh -i pod5_dir -o output_dir
```

2. Adapter Trimming
```bash
./scripts/02_porechop.sh -i fastq_dir -o trimmed_dir
```

3. Length Filtering
```bash
./scripts/03_nanofilt.sh -i trimmed_dir -o filtered_dir
```

4. Read Quality Metrics
```bash
./scripts/09_nanostat.sh -i filtered_dir -o metrics_dir
```

5. Taxonomic Classification
```bash
./scripts/04_kraken2.sh -i filtered_dir -d kraken2_db -o classified_dir
```

6. Assembly
```bash
./scripts/05_flye.sh -i filtered_dir -o assembly_dir
```

7. Read Mapping
```bash
./scripts/06_minimap2.sh -r filtered_dir -a assembly_dir -o mapped_dir
```

8. Assembly Polishing
```bash
./scripts/07_racon.sh -r filtered_dir -s mapped_dir -a assembly_dir -o polished_dir
```

9. Assembly Statistics
```bash
./scripts/10_assembly_stats.sh -i polished_dir -o stats_dir
```

## Pipeline Components

### 1. Quality Control
- Adapter removal using Porechop
- Length filtering using NanoFilt
- Read statistics using NanoStat

### 2. Taxonomic Classification
- Classification of raw reads using Kraken2
- Generation of taxonomic reports

### 3. Assembly and Polishing
- Metagenome assembly using Flye
- Read mapping using Minimap2
- Assembly polishing using Racon
- Assembly statistics using assembly-stats

## Directory Structure
```
project/
├── scripts/
│   ├── 01_basecalling.sh
│   ├── 02_porechop.sh
│   ├── 03_nanofilt.sh
│   ├── 04_kraken2.sh
│   ├── 05_flye.sh
│   ├── 06_minimap2.sh
│   ├── 07_racon.sh
│   └── 08_metamdbg.sh
├── processing/
│   ├── porechop/      # Adapter-trimmed reads
│   ├── nanofilt/      # Length-filtered reads
│   ├── kraken2/       # Taxonomic classification
│   ├── flye/          # Flye assemblies
│   ├── minimap2/      # Read alignments
│   ├── racon/         # Polished assemblies
│   └── mdbg/          # MetaMDBG assemblies
└── README.md
```
