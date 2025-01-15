# Wetland Environment Zoonotic and AMR Assessment Methodology

## Project Overview
This pilot project implements an integrated approach combining artificial intelligence-based remote sensing, non-invasive environmental sampling, and advanced genomic analyses to assess zoonotic potential and antimicrobial resistance (AMR) in wetland environments along the East Atlantic Flyway.

## Methods

### 1. AI-based Remote Sensing Analysis

#### Site Selection Criteria
- Focus on wetlands along the East Atlantic Flyway
- Prioritization based on:
  - Migratory bird importance
  - Varying levels of anthropogenic impact
  - Accessibility for sampling
  - Historical monitoring data availability

#### Satellite Data Collection
- **Primary Data Sources:**
  - Sentinel-1 (SAR imaging)
  - Sentinel-2 (multispectral imaging)
- **Temporal Resolution:** 5-day revisit time
- **Spatial Resolution:** 10m - 60m depending on bands

### 2. Non-invasive Environmental Sampling

#### Study Sites
- **Geographic Coverage:**
  - Germany (Northern region)
  - France (Central region)
  - Spain (Southern region)
- **Site Classification:**
  - 2 low anthropogenic pressure wetlands per country
  - 2 high anthropogenic pressure wetlands per country
  - Total: 12 study sites

#### Sample Collection Methods

##### Active Sampling
1. **Water Samples**
   - Volume: 1L per sample
   - Frequency: Daily collection
   - Filtration: 0.22μm sterilized filters

2. **Air Samples**
   - Technology: Liquid impingement
   - Collection medium: Triton solution
   - Duration: 3 hours per sample
   - Flow rate: 300 L/min

##### Passive Sampling
1. **Water Samples**
   - Device: Torpedo-shaped collectors
   - Membrane: Electronegative charged
   - Exposure time: 4 days
   - Depth: 0.5m below surface

2. **Air Samples**
   - Equipment: Modified Wilson and Cook (MWAC) samplers
   - Height: 1.5m above ground
   - Duration: 4 days
   - Wind direction monitoring

##### Avian Fecal Sampling
- Quantity: 30 samples per wetland
- Collection: Sterile swabs
- GPS coordinates recorded for each sampling site

### 3. Data Generation and Analysis

#### Nucleic Acid Extraction
- **Protocol:** Optimized for environmental matrices
- **Quality Control:**
  - Qubit 4.0 quantification
  - Internal extraction controls

#### Avian Influenza Virus (AIV) Analysis
1. **Initial Screening**
   - Method: Real-time RT-PCR
   - Target: M-segment
   - Controls: Positive and negative controls included
   - Ct value cutoff: ≤35 for positive samples

2. **Shotgun Sequencing**
   - Platform: Oxford Nanopore
   - Chemistry: R10.4.1
   - Library prep: RBK114.24

#### Virome Analysis
1. **Viral Enrichment**
   - Method: SISPA (Sequence-independent, single-primer amplification)
   - Primer design: Random hexamers with known tag sequence

2. **Sequencing**
   - Platform: Oxford Nanopore
   - Chemistry: R10 Direct RNA

3. **Bioinformatics**
   - TBD

#### Metagenomic Analysis
1. **Sequencing Approach**
   - Platform: Oxford Nanopore
   - Chemistry: R10.4.1
   - Kit: RBK114-24 Rapid Barcoding

3. **Taxonomic Classification**
   - Primary tools:
     - Kraken2
   - Database: NCBI nt-core
  
4. **Functional Analysis**
   - Databases:
     - KEGG
     - COG
     - eggNOG
     - NR
     - AMRFinderPlus
     - Plasmidfinder
   - Gene family analysis
   - Resistance mechanism classification
   - Mobile genetic element detection
   - Novel AMR gene prediction

## Data Integration and Analysis
- Statistical analysis of correlations between:
  - Land use patterns
  - Microbial diversity
  - AMR prevalence
  - Zoonotic potential
- Geospatial mapping of results

## Quality Control Measures
- Field blanks and technical replicates
- Negative controls through all processing steps
- Cross-validation between different analysis methods

## Expected Outcomes
1. Novel surveillance methodology framework
2. Baseline data for wetland microbiome composition
3. AMR prevalence maps
4. Zoonotic potential assessment tools
5. Recommendations for wetland management

## Data Management and Sharing
- Raw data deposited in appropriate databases (ENA/SRA)
- Analysis pipelines shared via GitHub
- Results published in open-access journals
- Data management plan following FAIR principles
