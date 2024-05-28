# Analysis of ribosome profiling (monosome/disome) or 5'P-seq data in U34 tRNA modification mutants. 


This repository includes scripts for analyzing the data in doi:10.1101/2024.02.27.582385. 


## Prerequisites
```
* FASTX-Toolkit
* STAR (>2.7.10)
* Bowtie (1.3.1)

```

## Prepare the genome/annotation index
### Generate ncRNA (rRNA, tRNA, snRNA...)

```
STAR --runThreadN 36 --runMode genomeGenerate --genomeDir Saccharomyces_cerevisiae.R64-1-1.86_changed_non_protein_coding_tRNACCA_STAR --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.86_changed_non_protein_coding_tRNACCA.fasta  --genomeSAindexNbases 4

```
