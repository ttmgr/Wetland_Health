# Nanopore Metagenomics Pipeline

A comprehensive pipeline for analyzing nanopore sequencing data from metagenomic samples. This pipeline includes quality control, assembly, taxonomic classification, and antimicrobial resistance gene detection.

## Features

- Adapter trimming with Porechop
- Read filtering with NanoFilt
- Optional barcode concatenation
- Read deduplication
- Taxonomic classification with Kraken2
- Metagenomic assembly with Flye
- Assembly polishing with Racon
- Antimicrobial resistance gene detection with AMRFinder
- SLURM job scheduling support
- Comprehensive logging
- Automatic environment management

## Prerequisites

- Conda/Mamba for environment management
- SLURM workload manager (for cluster deployment)
- Sufficient storage space for intermediate files

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/nanopore-metagenomics-pipeline.git
cd nanopore-metagenomics-pipeline
```

2. Create the required conda environments:

```bash
# Porechop environment
mamba create -n porechop_env -c bioconda porechop

# NanoFilt environment
mamba create -n nanofilt_env -c bioconda nanofilt

# Kraken2 environment
mamba create -n kraken2_env -c bioconda kraken2

# Flye environment
mamba create -n flye_env -c bioconda flye

# Minimap2 environment
mamba create -n minimap2_env -c bioconda minimap2

# Samtools environment
mamba create -n samtools_env -c bioconda samtools

# Racon environment
mamba create -n racon_env -c bioconda racon

# AMRFinder environment
mamba create -n amrfinderplus_env -c bioconda amrfinderplus

# Python environment for custom scripts
mamba create -n python_env python=3.9
```

3. Download required databases:
```bash
# Kraken2 database
kraken2-build --standard --threads 28 --db /path/to/kraken2_db

# AMRFinder database
amrfinder_update --force_update
```

## Usage

### Basic Usage

```bash
./metagenomics_pipeline.sh -i /path/to/fastq -o /path/to/output -k /path/to/kraken_db -a /path/to/amrfinder_db
```

### Required Parameters

- `-i, --input-dir`: Directory containing input FASTQ files
- `-o, --output-dir`: Directory for output files
- `-k, --kraken-db`: Path to Kraken2 database
- `-a, --amrfinder-db`: Path to AMRFinder database

### Optional Parameters

- `-t, --threads`: Number of threads (default: 28)
- `-m, --min-length`: Minimum read length for filtering (default: 100)
- `-c, --concat`: Barcode pairs to concatenate (optional)
- `--mem`: Memory limit (default: 250G)
- `--time`: Time limit (default: 50:00:00)

### Barcode Concatenation

To concatenate specific barcode pairs, use the `-c` option with comma-separated pairs:

```bash
./metagenomics_pipeline.sh -i /path/to/fastq -o /path/to/output -k /path/to/kraken_db -a /path/to/amrfinder_db -c barcode01,barcode02:barcode03,barcode04
```

This will concatenate barcode01 with barcode02, and barcode03 with barcode04.

## Pipeline Stages

1. **Quality Control**
   - Adapter trimming (Porechop)
   - Length filtering (NanoFilt)
   - Read deduplication (custom Python script)

2. **Taxonomic Classification**
   - Initial read classification (Kraken2)
   - Post-assembly contig classification (Kraken2)

3. **Assembly and Polishing**
   - Metagenomic assembly (Flye)
   - Read mapping (Minimap2)
   - Assembly polishing (Racon)

4. **Antimicrobial Resistance Analysis**
   - AMR gene detection on reads (AMRFinder)
   - AMR gene detection on assembled contigs (AMRFinder)

## Output Structure

```
output_directory/
├── logs/
├── porechop/
├── nanofilt/
├── kraken2_core/
├── kraken2_contigs_core/
├── flye/
├── minimap2/
├── samtools/
├── racon/
├── fasta_output/
├── amrfinder_output_nanofilt/
└── amrfinder_output_racon/
```

## Logging

The pipeline generates detailed logs in the `logs` directory:
- `pipeline.log`: Main pipeline progress
- SLURM output and error logs for each job

## Error Handling

The pipeline includes comprehensive error handling:
- Checks for required input files and directories
- Validates tool outputs
- Provides detailed error messages
- Maintains execution logs

## Resource Requirements

- RAM: Default 250GB (adjustable with --mem)
- Runtime: Default 50 hours (adjustable with --time)
- Storage: Varies with input data size (typically 3-5x input size)
- Threads: Default 28 (adjustable with -t)

## Contributing

Feel free to submit issues, fork the repository, and create pull requests for any improvements.

## License

[Your chosen license]

## Acknowledgments

This pipeline incorporates several open-source tools:
- [Porechop](https://github.com/rrwick/Porechop)
- [NanoFilt](https://github.com/wdecoster/nanofilt)
- [Kraken2](https://github.com/DerrickWood/kraken2)
- [Flye](https://github.com/fenderglass/Flye)
- [Minimap2](https://github.com/lh3/minimap2)
- [Racon](https://github.com/isovic/racon)
- [AMRFinder](https://github.com/ncbi/amr)
