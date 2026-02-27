# Bulk RNA Sequencing Exercise

>


## Summary of Contents:



## 1. Installation and Environment Setup:

```bash
#!/usr/bin/bash

# installing nextflow in conda env
conda create -n nextflow bioconda::nextflow
conda activate nextflow

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


## 2. Working directory setup:

```



```

## 3. Nextflow nf-core.

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

> make samplesheet.csv

| sample | fastq_1 | fastq_2 | strandedness |
|--------|---------|---------|--------------|
| SRR3222409_KO | SRR3222409-19_1.fq.gz | SRR3222409-19_2.fq.gz | auto |
| SRR3222410_KO | SRR3222410-19_1.fq.gz | SRR3222410-19_2.fq.gz | auto |
| SRR3222411_KO | SRR3222411-19_1.fq.gz | SRR3222411-19_2.fq.gz | auto |
| SRR3222412_WT | SRR3222412-19_1.fq.gz | SRR3222412-19_2.fq.gz | auto |
| SRR3222413_WT | SRR3222413-19_1.fq.gz | SRR3222413-19_2.fq.gz | auto |
| SRR3222414_WT | SRR3222414-19_1.fq.gz | SRR3222414-19_2.fq.gz | auto |

## 4. Downloading reference genome.

```bash
#!/usr/bin/bash

wget ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz

wget ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz


```

## 5. Running Nextflow.

```bash
#!/usr/bin/bash

nextflow run nf-core/rnaseq --input samplesheet.csv --outdir ./results --fasta ./reference/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz --gtf ./reference/Mus_musculus.GRCm38.99.gtf.gz -profile singularity --skip_linting

```
## 6. Resolving Duplicate header for read.
ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (SRR3222410_KO)'

Caused by:
  Process `NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (SRR3222410_KO)` terminated with an error exit status (1)


Command executed:

  fq lint \
      --disable-validator P001 \
      SRR3222410-19_1.fq.gz SRR3222410-19_2.fq.gz > SRR3222410_KO.fq_lint.txt
  
  cat <<-END_VERSIONS > versions.yml
  "NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT":
      fq: $(echo $(fq lint --version | sed 's/fq-lint //g'))
  END_VERSIONS

Command exit status:
  1

Command output:
  (empty)

Command error:
  INFO:    Environment variable SINGULARITYENV_NXF_TASK_WORKDIR is set, but APPTAINERENV_NXF_TASK_WORKDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  SRR3222410-19_1.fq.gz:333:1: [S007] DuplicateNameValidator: duplicate name: '@SRR3222410.10880'

Work dir:
  /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/work/4a/7576f78e9d0ca437f04c4fe7f76edf

Container:
  /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/work/singularity/depot.galaxyproject.org-singularity-fq-0.12.0--h9ee0642_0.img

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

 -- Check '.nextflow.log' file for details
ERROR ~ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting

 -- Check '.nextflow.log' file for details

```bash
#!/usr/bin/bash

# check if file is corrupted
gzip -t /home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/resources/SRR3222410-19_1.fq.gz

#
zgrep -n "@SRR3222410.10880" SRR3222410-19_1.fq.gz

#
zcat SRR3222410-19_1.fq.gz | sed -n '329,336p'

```
 

 ## R.

 ```
 colnames
[1] "gene_name"     "SRR3222409_KO" "SRR3222410_KO" "SRR3222411_KO" "SRR3222412_WT"
[6] "SRR3222413_WT" "SRR3222414_WT"

head

                   gene_name SRR3222409_KO SRR3222410_KO SRR3222411_KO SRR3222412_WT
ENSMUSG00000001750 "Tcirg1"  " 151.000"    " 129.001"    " 136.000"    "  242.000"  
ENSMUSG00000003053 "Cyp2c29" "   0.000"    "   0.000"    "   0.000"    "    0.000"  
ENSMUSG00000003228 "Grk5"    " 164.000"    " 156.000"    " 186.000"    "  242.000"  
ENSMUSG00000003555 "Cyp17a1" "   1.000"    "   0.000"    "   0.000"    "    0.000"  
ENSMUSG00000003559 "As3mt"   "  73.000"    "  54.000"    "  92.001"    "  113.000"  
ENSMUSG00000003680 "Taf6l"   "  66.000"    "  37.000"    "  50.000"    "   87.000"  
                   SRR3222413_WT SRR3222414_WT
ENSMUSG00000001750 "  175.000"   "  165.000"  
ENSMUSG00000003053 "    0.000"   "    0.000"  
ENSMUSG00000003228 "  225.000"   "  239.000"  
ENSMUSG00000003555 "    0.000"   "    0.000"  
ENSMUSG00000003559 "   78.999"   "  105.000"  
ENSMUSG00000003680 "   64.000"   "   53.000"  

dim

[1] 1394    7

class

[1] "matrix" "array" 



After changing count matrix:

[1] "SRR3222409_KO" "SRR3222410_KO" "SRR3222411_KO" "SRR3222412_WT" "SRR3222413_WT" "SRR3222414_WT"
     SRR3222409_KO SRR3222410_KO SRR3222411_KO SRR3222412_WT SRR3222413_WT SRR3222414_WT
[1,]           151           129           136           242           175           165
[2,]             0             0             0             0             0             0
[3,]           164           156           186           242           225           239
[4,]             1             0             0             0             0             0
[5,]            73            54            92           113            79           105
[6,]            66            37            50            87            64            53
[1] 1394    6
[1] "matrix" "array" 

[1] "SRR3222409_KO" "SRR3222410_KO" "SRR3222411_KO" "SRR3222412_WT" "SRR3222413_WT" "SRR3222414_WT"
                   SRR3222409_KO SRR3222410_KO SRR3222411_KO SRR3222412_WT SRR3222413_WT SRR3222414_WT
ENSMUSG00000001750           151           129           136           242           175           165
ENSMUSG00000003053             0             0             0             0             0             0
ENSMUSG00000003228           164           156           186           242           225           239
ENSMUSG00000003555             1             0             0             0             0             0
ENSMUSG00000003559            73            54            92           113            79           105
ENSMUSG00000003680            66            37            50            87            64            53
[1] 1394    6
[1] "matrix" "array" 
 ```