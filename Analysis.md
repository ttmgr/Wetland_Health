# Nanopore Metagenomics Analysis Pipeline

## System Requirements
- RAM: ≥16GB (32GB recommended for large datasets)
- Storage: Sufficient space for raw data + ~3x for processing
- CPU: Modern multi-core processor (≥8 cores recommended)
- OS: Linux (tested on Ubuntu 20.04 LTS)
- Conda for environment management

## Quick Start
```bash
# Clone the repository
git clone https://github.com/yourusername/repo-name.git
cd repo-name

# Create and activate conda environment
conda env create -f metagenomics.yml
conda activate metagenomics

# Download required databases
./setup/download_databases.sh
```

## Analysis Workflow

⚠️ **Important**: This pipeline must be run in the exact order presented below to ensure reproducibility. Each step depends on the output of previous steps.

## Pipeline Overview
1. Read Processing and Quality Control (Steps 1-4)
2. Taxonomic Classification (Step 5)
3. Assembly and Analysis (Steps 6-9)
4. Data Format Conversion (Step 10)
5. AMR Detection (Steps 11-12)

## Pipeline Components and Expected Runtimes
All runtimes are estimated based on a 10GB raw nanopore dataset.

### 1. Read Processing and Quality Control
- Basecalling with Dorado v0.4.0 ([01_basecalling.sh](bash_scripts/01_basecalling.sh))
  - Runtime: ~4-6 hours with GPU
- Adapter removal with Porechop v0.2.4 ([02_porechop.sh](bash_scripts/02_porechop.sh))
  - Runtime: ~30-45 minutes
- Length filtering with NanoFilt v2.8.0 ([03_nanofilt.sh](bash_scripts/03_nanofilt.sh))
  - Runtime: ~15-20 minutes
- Quality metrics with NanoStat v1.6.0 ([04_readmetrics.sh](bash_scripts/04_readmetrics.sh))
  - Runtime: ~10-15 minutes

### 2. Taxonomic Classification
- Read-based classification using Kraken2 v2.1.2 ([05_kraken2_reads.sh](bash_scripts/05_kraken2_reads.sh))
  - Runtime: ~2-3 hours

### 3. Assembly and Analysis
- Assembly with Flye v2.9.2 ([06_flye.sh](bash_scripts/06_flye.sh))
  - Runtime: ~4-6 hours
- Read mapping with Minimap2 v2.24 ([07_minimap2.sh](bash_scripts/07_minimap2.sh))
  - Runtime: ~1-2 hours
- Assembly polishing with Racon v1.5.0 ([08_racon.sh](bash_scripts/08_racon.sh))
  - Runtime: ~2-3 hours
- Assembly statistics generation ([09_assemblystats.sh](bash_scripts/09_assemblystats.sh))
  - Runtime: ~5-10 minutes

### 4. Data Format Conversion
- FASTQ to FASTA conversion using SeqKit v2.3.0 ([10_seqkit.sh](bash_scripts/10_seqkit.sh))
  - Runtime: ~10-15 minutes

### 5. AMR Detection
- AMR gene detection in reads using AMRFinderPlus v3.11.18 ([11_amrfinderplus_reads.sh](bash_scripts/11_amrfinderplus_reads.sh))
  - Runtime: ~2-3 hours
- AMR gene detection in assemblies ([12_amrfinderplus_contigs.sh](bash_scripts/12_amrfinderplus_contigs.sh))
  - Runtime: ~1-2 hours

Total expected runtime: ~18-25 hours for a 10GB dataset

## Detailed Usage

### 1. Basecalling
```bash
./bash_scripts/01_basecalling.sh -i pod5_dir -o output_dir
```
Converts raw ONT signals to FASTQ format and performs demultiplexing.
Required GPU: NVIDIA with ≥8GB VRAM

### 2. Adapter Trimming
```bash
./bash_scripts/02_porechop.sh -i fastq_dir -o trimmed_dir
```
Removes adapter sequences from reads.

[... rest of the usage instructions remain the same ...]

## Output Structure
```
processing/
├── porechop/      # Adapter-trimmed reads
├── nanofilt/      # Length-filtered reads
├── readmetrics/   # Read quality metrics
├── kraken2/       # Taxonomic classification
├── flye/          # Assemblies
├── minimap2/      # Read alignments
├── racon/         # Polished assemblies
├── assemblystats/ # Assembly statistics
├── fasta/         # Converted FASTA files
└── amrfinder/     # AMR detection results
```

## Required Databases and Setup

### Kraken2 Database
```bash
# Download and setup standard Kraken2 database
kraken2-build --standard --db kraken2_db
# Required storage: ~150GB
# Download time: ~2-4 hours depending on connection
```

### AMRFinder Database
The AMRFinder database is automatically downloaded during installation but can be updated:
```bash
amrfinder_update --force_update
# Required storage: ~20GB
# Download time: ~15-30 minutes
```

## Troubleshooting

### Common Issues
1. Insufficient memory during assembly
   - Solution: Reduce thread count or increase swap space
2. GPU errors during basecalling
   - Solution: Ensure CUDA drivers are up to date
3. Database download failures
   - Solution: Use provided mirror links in setup/download_databases.sh
```
