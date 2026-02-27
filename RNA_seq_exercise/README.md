# Bulk RNA Sequencing Exercise

> Visualizing differential expression of bulk-RNA-seq transcript levels in WT and YAP/TAzdKO Schwann cells. 
> https://nbisweden.github.io/workshop-ngsintro/2511/topics/rnaseq/lab_rnaseq.html#pca-plot
> 


## Summary of Contents:

1. 
2. 
3. 
4. 
5. 
6. 

## 1. Nextflow nf-core.

> https://nf-co.re/rnaseq/3.22.2/
> https://github.com/nf-core/rnaseq

```
1. Merge re-sequenced FastQ files (cat)
2. Auto-infer strandedness by subsampling and pseudoalignment (fq, Salmon)
3. Read QC (FastQC)
4. UMI extraction (UMI-tools)
5. Adapter and quality trimming (Trim Galore!)
6. Removal of genome contaminants (BBSplit)
7. Removal of ribosomal RNA (SortMeRNA)
8. Choice of multiple alignment and quantification routes (For STAR the sentieon implementation can be chosen):
- STAR -> Salmon
- STAR -> RSEM
- HiSAT2 -> NO QUANTIFICATION
9. Sort and index alignments (SAMtools)
10. UMI-based deduplication (UMI-tools)
11. Duplicate read marking (picard MarkDuplicates)
12. Transcript assembly and quantification (StringTie)
13. Create bigWig coverage files (BEDTools, bedGraphToBigWig)
14. Extensive quality control:
- RSeQC
- Qualimap
- dupRadar
- Preseq
- DESeq2
- Kraken2 -> Bracken on unaligned sequences; optional
15. Pseudoalignment and quantification (Salmon or ‘Kallisto’; optional)
16. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks (MultiQC, R)

```

## 2. Installation and Environment Setup:

```bash
#!/usr/bin/bash

# installing nextflow in conda env
conda create -n nextflow bioconda::nextflow
conda activate nextflow

# exporting conda details
conda env export > nextflow_environment.yml

# nextflow info
nextflow info

#   Version: 25.10.4 build 11173     
#   Created: 10-02-2026 15:17 UTC (16:17 CEST)
#   System: Linux 6.14.5-100.fc40.x86_64        
#   Runtime: Groovy 4.0.28 on OpenJDK 64-Bit Server VM 23.0.2-internal-adhoc.conda.src
#   Encoding: UTF-8 (UTF-8)

# initiating rnaseq workflow
nextflow run nf-core/rnaseq --help

```


## 3. Working directory setup:

```
.
├── envs
├── reference
├── results
└── scripts

```

## 4. Setting up necessary input files.

> make samplesheet.csv (no spaces and comma separated)

| sample | fastq_1 | fastq_2 | strandedness |
|--------|---------|---------|--------------|
| SRR3222409_KO | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222409-19_1.fq.gz | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222409-19_2.fq.gz | auto |
| SRR3222410_KO | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222410-19_1.fq.gz | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222410-19_2.fq.gz | auto |
| SRR3222411_KO | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222411-19_1.fq.gz | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222411-19_2.fq.gz | auto |
| SRR3222412_WT | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222412-19_1.fq.gz | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222412-19_2.fq.gz | auto |
| SRR3222413_WT | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222413-19_1.fq.gz | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222413-19_2.fq.gz | auto |
| SRR3222414_WT | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222414-19_1.fq.gz | /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222414-19_2.fq.gz | auto |

## 5. Downloading reference genome.

```bash
#!/usr/bin/bash

# placed into reference directory
wget ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz

wget ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

```

## 6. Running Nextflow.

```bash
#!/usr/bin/bash

# from project root 

# had to skip fastq validation step (with skip_linting) because there are duplicate headers in some fastq files
nextflow run nf-core/rnaseq --input samplesheet.csv --outdir ./results --fasta ./reference/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz --gtf ./reference/Mus_musculus.GRCm38.99.gtf.gz -profile singularity --skip_linting

```
 
 ## 7. Using R to visualize count data.

> scripts used:\
> dge_02.R\
> This script will do:\
1. Run DESeq analysis on "/results/star_salmon/salmon.merged.gene_counts.tsv"
2. Resort by padj and produce a most differential genes tables with the 0.05 padj and log2FC of 1 cutoff.
3. Perform a VST tranformation.
4. Produce a Volcano plot.
5. Produce a Heatmap. 



> ***Other scripts are included as inspiration for the final script and were not used for the analysis.***