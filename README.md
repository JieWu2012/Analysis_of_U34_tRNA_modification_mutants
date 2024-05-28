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
* In `Read_preprecessing.sh`, change the adapter sequence and nucleotides that need to be clipped or trimmed from your sequences; change the parameters of STAR as you need and run: 

```
./Read_preprocessing.sh 
```
After running, the clean read fastq files are in the folder `ncRNA_rm`. 

## Mapping

Run the following command to map clean reads to CDS 
```
./Mapping.sh
```
After running, the periodicity frame plots for each sample will be generated. Select the read lengths with high abundance; select the dominant frames for each read length and corresponding offsets for each frame (if initiating ribosome footprints of frame 0 is at -12 in frame plots: `'0' -> '16', '1' -> '15', '2' -> '17'`). Generate the frame selection file for next step, see demo `Demo_FrameSelection.txt`. The second column is the sample name. The last column is either "wt" or "mutant". 

## Downstream analysis

```
./Downstream_analysis.sh SacCer3_CDS_ex21.fa Output_folder Demo_FrameSelection.txt
```
After running, it generates `gene.txt` for differential expression analysis; `all_sample_Rel_occ_filter0site_median.txt` for calculating vulnerbility score between mutant and wt; `dicodon.txt`for 2D di-codon plot. 
