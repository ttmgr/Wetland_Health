# Analysis Workflow

## Pipeline Components

### 1. Read Processing and Quality Control
- Basecalling with Dorado ([01_basecalling.sh](bash_scripts/01_basecalling.sh))
- Adapter removal with Porechop ([02_porechop.sh](bash_scripts/02_porechop.sh))
- Length filtering with NanoFilt ([03_nanofilt.sh](bash_scripts/03_nanofilt.sh))
- Quality metrics with NanoStat ([04_readmetrics.sh](bash_scripts/04_readmetrics.sh))

### 2. Taxonomic Classification
- Read-based classification using Kraken2 ([05_kraken2_reads.sh](bash_scripts/05_kraken2_reads.sh))

### 3. Assembly and Analysis
- Assembly with Flye ([06_flye.sh](bash_scripts/06_flye.sh))
- Read mapping with Minimap2 ([07_minimap2.sh](bash_scripts/07_minimap2.sh))
- Assembly polishing with Racon ([08_racon.sh](bash_scripts/08_racon.sh))
- Assembly statistics generation ([09_assemblystats.sh](bash_scripts/09_assemblystats.sh))

### 4. Data Format Conversion
- FASTQ to FASTA conversion ([10_seqkit.sh](bash_scripts/10_seqkit.sh))

### 5. AMR Detection
- AMR gene detection in reads ([11_amrfinderplus_reads.sh](bash_scripts/11_amrfinderplus_reads.sh))
- AMR gene detection in assemblies ([12_amrfinderplus_contigs.sh](bash_scripts/12_amrfinderplus_contigs.sh))

## Detailed Usage

### 1. Basecalling
```bash
./bash_scripts/01_basecalling.sh -i pod5_dir -o output_dir
```
Converts raw ONT signals to FASTQ format and performs demultiplexing.

### 2. Adapter Trimming
```bash
./bash_scripts/02_porechop.sh -i fastq_dir -o trimmed_dir
```
Removes adapter sequences from reads.

### 3. Length Filtering
```bash
./bash_scripts/03_nanofilt.sh -i trimmed_dir -o filtered_dir
```
Filters reads based on length criteria.

### 4. Read Quality Metrics
```bash
./bash_scripts/04_readmetrics.sh -i filtered_dir -o metrics_dir
```
Generates quality metrics for filtered reads.

### 5. Taxonomic Classification
```bash
./bash_scripts/05_kraken2_reads.sh -i filtered_dir -d kraken2_db -o classified_dir
```
Performs taxonomic classification of reads.

### 6. Assembly
```bash
./bash_scripts/06_flye.sh -i filtered_dir -o assembly_dir
```
Assembles reads into contigs.

### 7. Read Mapping
```bash
./bash_scripts/07_minimap2.sh -r filtered_dir -a assembly_dir -o mapped_dir
```
Maps reads back to assemblies.

### 8. Assembly Polishing
```bash
./bash_scripts/08_racon.sh -r filtered_dir -s mapped_dir -a assembly_dir -o polished_dir
```
Polishes assemblies using mapped reads.

### 9. Assembly Statistics
```bash
./bash_scripts/09_assemblystats.sh -i polished_dir -o stats_dir
```
Generates assembly quality metrics.

### 10. FASTQ to FASTA Conversion
```bash
./bash_scripts/10_seqkit.sh -i filtered_dir -o fasta_dir
```
Converts FASTQ files to FASTA format.

### 11. AMR Detection (Reads)
```bash
./bash_scripts/11_amrfinderplus_reads.sh -i fasta_dir -d amrfinder_db -o amr_reads_dir
```
Detects AMR genes in filtered reads.

### 12. AMR Detection (Assembly)
```bash
./bash_scripts/12_amrfinderplus_contigs.sh -i polished_dir -d amrfinder_db -o amr_assembly_dir
```
Detects AMR genes in polished assemblies.

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

## Required Databases

### Kraken2 Database
```bash
# Download and setup standard Kraken2 database
kraken2-build --standard --db kraken2_db
```

### AMRFinder Database
The AMRFinder database is automatically downloaded during installation but can be updated:
```bash
amrfinder_update --force_update
```
