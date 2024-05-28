# Analysis of ribosome profiling (monosome/disome) or 5'P-seq data in U34 tRNA modification mutants. 


This repository includes scripts for analyzing the data in doi:10.1101/2024.02.27.582385. 


## Prerequisites

* FASTX-Toolkit
* STAR (>2.7.10)
* Bowtie (1.3.1)


## Prepare the genome/annotation index
### Generate ncRNA (rRNA, tRNA, snRNA...)

```
STAR --runThreadN 18 --runMode genomeGenerate --genomeDir ncRNA_STAR --genomeFastaFiles SacCer3_ncRNA.fa  --genomeSAindexNbases 4
```
