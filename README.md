# Wetland Environment Zoonotic and AMR Assessment Pipeline

## Overview
This repository contains the bioinformatics pipeline for analyzing environmental samples from wetlands along the East Atlantic Flyway. The pipeline processes Oxford Nanopore sequencing data to assess zoonotic potential and antimicrobial resistance (AMR) in wetland environments.

## Project Scope
- Analysis of environmental samples from wetland sites
- Processing of water, air, and avian fecal samples
- Integration of metagenomic and AMR analysis
- Focus on wetlands with varying anthropogenic pressure

## Pipeline Components
- Read processing and quality control
- Taxonomic classification
- Metagenomic assembly and analysis
- AMR gene detection
- Quality metrics generation

## Prerequisites
### Hardware Requirements
- High-performance computing environment recommended
- Sufficient storage for sequencing data and databases (>1TB recommended)
- Minimum 32GB RAM recommended

### Software Requirements
All tools can be installed via conda:
```bash
# Create and activate conda environment
conda create -n wetland_pipeline python=3.9
conda activate wetland_pipeline

# Install required tools
conda install -c bioconda dorado porechop nanofilt nanostat kraken2 flye \
    minimap2 racon assembly-stats ncbi-amrfinderplus seqtk
conda install -c conda-forge parallel
```

## Quick Start
1. Clone this repository
```bash
git clone https://github.com/yourusername/wetland_health.git
cd wetland_health
```

2. Make scripts executable
```bash
chmod +x bash_scripts/*.sh
```

3. See `ANALYSIS.md` for detailed workflow instructions

## Repository Structure
```
wetland_health/
├── bash_scripts/   # Analysis pipeline scripts
├── processing/     # Output directories
└── docs/          # Documentation
```

## Documentation
- `README.md`: Project overview and setup instructions
- `ANALYSIS.md`: Detailed analysis workflow and usage guide
- `METHODOLOGY.md`: Comprehensive project methodology

## Data Management
- Raw sequencing data will be deposited in ENA/SRA
- Analysis results follow FAIR principles
- Pipeline and scripts available in this repository

## Contributing
Contributions are welcome! Please feel free to submit issues and enhancement requests.
