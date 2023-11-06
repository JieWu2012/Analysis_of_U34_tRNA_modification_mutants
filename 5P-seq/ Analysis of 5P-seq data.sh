# Analysis of 5P-seq data.



# Remove adapter sequence. Note: this step does not discard untrimmed reads. 
for i in *R1_*.fastq;do echo $i;a=${i/.fastq/};cutadapt -j 15 -m 10  -a GATCGGAAGAGCACACGTCTGAACTCCAGTC  $i > $a"_cl_keepuntrimmed.fastq" 2> $a"_report_keepuntrimmed.txt";done;

# Extract 8nt umi at the 5' end of reads and put it into read name. 
for i in *keepuntrimmed.fastq; do echo $i;a=${i/\.fastq/};umi_tools extract --extract-method=string --bc-pattern=NNNNNNNN -L $a".log"  --stdin=$i --stdout $a"_extracted.fastq";done


# Remove ncRNA. 
for i in *_extracted.fastq ;do echo $i;a=${i/\.fastq/};b=${a/*\//};STAR --runThreadN 30 --genomeDir ~/annotation/Saccharomyces_cerevisiae.R64-1-1.86_changed_non_protein_coding_tRNACCA_STAR  --readFilesIn $i --outFilterMultimapNmax 10 --outSAMmultNmax 1 --outFilterMismatchNoverReadLmax 0.05  --alignEndsType EndToEnd --outFileNamePrefix rm_rRNA/$b"_2_rRNA_"  --outReadsUnmapped Fastx --outSAMtype BAM Unsorted;done;


# Map to transcriptome rp_analysis script. Make sure there are bam file in the output. -X: Read length (in my case 60). -I: 20. -E: 21 (extention in yor annotaiton.)
for i in rm_rRNA/*.fastq;do echo $i;./rp_analysis_for_short_reads_+8910_frag60.sh -M Mapping -Q $i  -T /home/jwu/annotation/extend_transcriptome/Saccharomyces_cerevisiae.R64-1-1.86_changed_CDS_non_overlappwithohters_protein_coding_extended_annotated_nooverlapwithTy  -O mapping_keepuntrimmed  -X 60 -I 10 -E 21;done;

# For generating the frame information at the 3' end of genes. 
for i in *.bam ;do echo $i;a=${i/\.bam/};b=${i/_R1*/};samtools view $i | awk 'BEGIN{OFS="\t"}NR==FNR{if(NR%2==1){sub(/>/,"",$1)
;hash[NR]=$1}else{len[hash[NR-1]]=length($0)}}NR!=FNR && $2==0 {sub(/M/,"",$6);if($0~/MD:Z:0[A-Z]/){if($6-1>=10 && $6-1<=60)print $4+1
-(len[$3]-21+1),$6-1}else{if($6>=10 && $6<=60)print $4-(len[$3]-21+1),$6}}' ~/annotation/extend_transcriptome/Saccharomyces_cerevisiae
.R64-1-1.86_changed_CDS_non_overlappwithohters_protein_coding_extended_annotated_nooverlapwithTy.fa  - | awk 'BEGIN{OFS="\t"}{hash[$0]
++}END{for(i in hash){print i,hash[i],"'$b'"}}' >> 5P_secondbatch_transcriptomebased_all_5end_stopcodon_RLchosed.txt;done;

# Use R code "5p_frame_information_3end.R" to generate frame information. 5P-seq has quite bad periodicity at 52 nt (=60-8)in my data.

cd rm_rRNA

# Sort bam and make index for bam.
for i in *.bam;do a=${i/\.bam/};echo $i;samtools sort $i > $a"_sorted.bam";samtools index $a"_sorted.bam";done

# make sure you install umi_tools
conda activate umi_tools

# Remove PCR duplication
for i in *_sorted.bam;do a=${i/\.bam/};echo $i;umi_tools dedup -I $i  --output-stats=$a"_deduplicated" -S $a"_dedup.bam";done

# Check the duplication level (read counts before and after deduplication)
for i in *_mapped.bam;do echo $i;a=${i/\.bam/};samtools view $i | wc -l; samtools view $a"_sorted_dedup.bam" | wc -l;done;

# Generate bed file for downstream analysis. 
for i in *_dedup.bam;do echo $i; a=${i/\.bam/};samtools view $i | awk  'BEGIN{FS="\t";OFS="\t"} $2==0{sub(/M/,"",$6);if($0~/MD:
Z:0[A-Z]/){p=substr($10,2,length($10));print $3,$4-1+1,$4-1+$6,"@"$1"|"p,"0","+"}else{print $3,$4-1,$4-1+$6,"@"$1"|"$10,"0","+"}}' > $
a".bed6";done;

# Generate raw counts for each sample. 
for i in *.bed6;do echo $i;a=${i/\.bed6/};awk '{hash[$1]++}END{for(i in hash){print i"\t"hash[i]}}' $i > $a"_count.txt";done;

# Generate raw count matrix. 
echo "Gene" > all_count_rmdup.txt

cut -f 1 *_count.txt | sort | uniq >> all_count_rmdup.txt

for i in *_R1*_count.txt;do echo $i;a=${i/_R1*/};awk 'NR==FNR{hash[$1]=$2}NR!=FNR{if(FNR==1){print $0"\t'$a'"}else if(hash[$1])
{print $0"\t"hash[$1]}else{print $0"\t0"}}' $i all_count_rmdup.txt > all_count_rmdup.txt_temp;mv all_count_rmdup.txt_temp all_count_rmdup.txt;done; 

# Use all_count_rmdup.txt for DE analysis with R code "5p_DE_analysis.R". 
