# Analysis of ribosome profiling (monosome/disome) or 5'P-seq data in U34 tRNA modification mutants. 


This repository includes scripts for analyzing the data in doi:10.1101/2024.02.27.582385. 


## Prerequisites

* FASTX-Toolkit
* STAR (>2.7.10)
* Bowtie (1.3.1)


## Prepare the genome/annotation index
### Generate STAR index for ncRNA (rRNA, tRNA, snRNA...)

```
STAR --runThreadN 18 --runMode genomeGenerate --genomeDir ncRNA_STAR --genomeFastaFiles SacCer3_ncRNA.fa  --genomeSAindexNbases 4
```
### Generate bowtie index for CDS with 21 nt extension

```
bowtie-build SacCer3_CDS_ex21.fa SacCer3_CDS_ex21
```
## Preprocessing the data

* Copy the scripts inside "Ribo-seq" folder to your ribo-seq data folder
* In "Read_preprecessing.sh", Change the adapter sequence and nucleotides that need to be clipped or trimmed from your sequences, change the parameters of STAR as you need and run: 

```
./Read_preprocessing.sh 
```
After running, the clean read fastq files are in the folder `ncRNA`
