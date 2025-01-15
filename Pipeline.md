# Metagenomic Analysis Pipeline for Environmental Samples

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
- Kraken2 (for taxonomic classification)
- Flye (for assembly)
- Minimap2 (for read mapping)
- Racon (for assembly polishing)
- MetaMDBG (for alternative assembly)

### Required Databases
- Kraken2 database

## Installation

### Using Conda
```bash
# Create conda environment
conda create -n metagenome_pipeline python=3.9
conda activate metagenome_pipeline

# Install required tools
conda install -c bioconda dorado
conda install -c bioconda porechop
conda install -c bioconda nanofilt
conda install -c bioconda kraken2
conda install -c bioconda flye
conda install -c bioconda minimap2
conda install -c bioconda racon
conda install -c bioconda metamdbg
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

4. Taxonomic Classification
```bash
./scripts/04_kraken2.sh -i filtered_dir -d kraken2_db -o classified_dir
```

5. Assembly
```bash
./scripts/05_flye.sh -i filtered_dir -o assembly_dir
```

6. Read Mapping
```bash
./scripts/06_minimap2.sh -r filtered_dir -a assembly_dir -o mapped_dir
```

7. Assembly Polishing
```bash
./scripts/07_racon.sh -r filtered_dir -s mapped_dir -a assembly_dir -o polished_dir
```

8. Alternative Assembly
```bash
./scripts/08_metamdbg.sh -i filtered_dir -o mdbg_dir
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
- Primary assembly using Flye
- Alternative assembly using MetaMDBG
- Read mapping using Minimap2
- Assembly polishing using Racon

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

## Parameters

### Default Settings
- Minimum read length: 100 bp
- Number of threads: 16
- Flye options: `--meta --nano-hq`
- Minimap2 options: `-ax map-ont`

These parameters can be adjusted using command-line options for each script.

## Output Files

### Quality Control
- Trimmed FASTQ files
- Filtered FASTQ files
- Quality reports

### Taxonomic Classification
- Kraken2 output files
- Taxonomic classification reports

### Assembly
- Assembled contigs (FASTA)
- Assembly graphs
- Assembly statistics

### Polishing
- Mapped reads (SAM/BAM)
- Polished assemblies (FASTA)

## Contributing
Feel free to submit issues and enhancement requests.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
[Your Name] - [your.email@example.com]
