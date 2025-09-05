# rRNA Vertebrate Analysis

This pipeline describes how avian host DNA is recovered from environmental 12S rRNA amplicon sequencing data.

## Tools and Reference Databases

- **OBITools4 – `obimultiplex` (v1.3.1)**: demultiplexes reads using a tab‑delimited tag file that defines the 9 bp primer tags for each sample.
- **Cutadapt (v4.2)**: removes forward and reverse primer sequences from reads.
- **VSEARCH (v2.21)**: open‑source sequence analysis toolkit used for quality filtering, dereplication, chimera removal and OTU clustering.
- **MIDORI2 mitochondrial database – “Unique 266” release**: curated set of vertebrate mitochondrial sequences compiled from GenBank and BOLD; used for taxonomic identification.
- **eBird**: global bird‑occurrence database used to confirm the geographic plausibility of species assignments.

## Pipeline

1. **Demultiplex basecalled reads**
   - `obimultiplex` assigns reads to samples based on their 9 bp tags, allowing up to two mismatches.

2. **Trim primer sequences**
   - `cutadapt` removes anchored forward and reverse primers from each read.

3. **Filter and cluster sequences with VSEARCH**
   - **Quality filter** by expected error (`--fastq_maxee 1.0`).
   - **Dereplicate** sequences and discard singletons (`--minuniquesize 2`).
   - **Remove chimeras** de novo (`--uchime_denovo`).
   - **Cluster** remaining reads into operational taxonomic units (OTUs) at 97 % identity.

4. **Assign taxonomy**
   - OTU centroids are globally aligned against the MIDORI2 “Unique 266” database. Hits must cover ≥80 % of the query and have ≥98 % identity. OTUs with <5 reads are discarded.

5. **Validate and tally species**
   - Species not reported in the sampling region (checked with eBird) are collapsed to higher taxonomic ranks or removed. Read counts from multiple OTUs mapping to the same retained taxon are summed per sample.

## Example command sequence

```bash
# 1. Demultiplex by 9 bp tags
obimultiplex --tags tagfile.txt --primerfile primers.fasta --allow-mismatch 2 \
    --input basecalled_reads.fastq --output demultiplexed.fastq

# 2. Remove primers
cutadapt -g ^PRIMER_FWD -G ^PRIMER_REV -o trimmed.fastq demultiplexed.fastq

# 3. Filter, dereplicate and cluster
vsearch --fastq_filter trimmed.fastq --fastq_maxee 1.0 --fastaout filtered.fasta
vsearch --derep_fulllength filtered.fasta --sizeout --minuniquesize 2 --output derep.fasta
vsearch --uchime_denovo derep.fasta --nonchimeras nochimeras.fasta
vsearch --cluster_size nochimeras.fasta --id 0.97 --centroids otus.fasta

# 4. Taxonomic assignment against MIDORI2
vsearch --usearch_global otus.fasta --db MIDORI2_unique_266.fasta --id 0.98 \
    --blast6out assignments.b6

# 5. Downstream processing
#   Summarize assignments, check against eBird, and aggregate read counts per taxon
```

