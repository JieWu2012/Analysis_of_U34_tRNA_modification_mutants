# Analysis of monosome profiling data


# Remove adapter sequence at the 3' end using fastx toolkit. (For wild type and ncs2Delp6D)

for i in *.fastq

do 

	echo $i

	a=${i/\.fastq/}

	fastx_clipper -a CTGTAGGCACCATCAAT -l 15 -c -Q 33  -i $i   -o $a"_clipped.fastq"

done

# Trim addtional 4 nt at the 3' end. 
# If there are additional nt at the 5' end, the reads need to be further trimmed using "-f". See manual of fastx toolkit. 

for i in *_clipped.fastq

do 

	echo $i

	a=${i/\.fastq/}

	fastx_trimmer -m 15 -t 4 -i $i -o $a"_trimmer_mt15.fastq"

	rm $i

done


# Fastqc

for i in *_trimmer_mt15.fastq

do 
	echo $i

	fastqc $i

done

# Remove non protein-coding reads

mkdir rRNA_rm

ncRNA_STAR=~/annotation/Saccharomyces_cerevisiae.R64-1-1.86_changed_non_protein_coding_tRNACCA_nopseudogene_STAR/

for i in data/*.fastq 

do 

	echo $i

	a=${i/\.fastq/}

	b=${a/*\//}

	STAR 	--runThreadN 20 \
			--genomeDir $ncRNA_STAR \
			--readFilesIn  \
			--outFilterMultimapNmax 30 \
			--outSAMmultNmax 1 \
			--outFilterMismatchNoverReadLmax 0.1  \
			--alignEndsType EndToEnd \
			--outFileNamePrefix rRNA_rm/$b"_toncRNA_"  \
			--outReadsUnmapped Fastx \
			--outSAMtype BAM Unsorted 
done

# Change file name

for i in *mate1;do echo $i; a=${i/\mate1/fastq}; mv $i $a; done;


# Map reads to CDS and generate periodicity plot. 

CDSbtIndex=~/annotation/extend_transcriptome/Saccharomyces_cerevisiae.R64-1-1.86_changed_CDS_non_overlappwithohters_protein_coding_extended_annotated_nooverlapwithTy

for i in rRNA_rm/*.fastq

do 
	echo $i
	./rp_analysis.sh -M Mapping \
					 -Q $i  \
					 -T $CDSbtIndex   \
					 -O mapping/  \
					 -X 34 \
					 -I 16 \
					 -E 21
done

# Generate codon plot (compare codon speed between KO and WT)
# This step will generate codon plot w/o log. 
# Frame selection file is made according to periodicity plot. Run "rp_analysis.sh -h" or check example to see how to make index plot. 

CDS_fasta=Saccharomyces_cerevisiae.R64-1-1.86_changed_CDS_non_overlappwithohters_protein_coding_extended_annotated_nooverlapwithTy.fa

Output=YourOutputFolder

FrameSelection=YourFrameSelectionFile

./rp_analysis.sh -M CodonPlot \
				 -A $CDS_fasta \
				 -D $FrameSelection \
				 -O $Output  \
				 -E 21 \
				 -F 15 \
				 -B 15 

# Generate raw count matrix for DE analysis. 

for i in *.sam

do 
	echo $i

	a=${i/\.sam/}

	awk '$1!~/^@/ && $3!="*"{

		if(!hash[$1]){

			count[$3]++
			hash[$1]=1
		}

	}END{

		for(i in count){

			print i"\t"count[i]

		}
	}' $i > $a"_count.txt"

 done


echo "Gene" > gene.txt

cut -f 1 *_count.txt | sort | uniq >> gene.txt 

for i in *_count.txt

do 

	echo $i

	a=${i/_count*.txt/}

	awk 'NR==FNR{

		hash[$1]=$2

	}NR!=FNR{

		if(FNR==1){

			print $0"\t'$a'"

		}else if(hash[$1]){

			print $0"\t"hash[$1]

		}else{

			print $0"\t0"
		}
	}' $i gene.txt  > gene.txt_temp

	mv gene.txt_temp gene.txt

done


# Calculate vulnerability score
# This step is performed in the result folder of codon plot analysis. 
 

echo -e  "Gene\tLoc\tCodon" > all_sample.txt

cat *excluded.txt >> all_sample.txt

for i in *trimmer*_offset.txt

do 

	echo $i

	a=${i/_S*/}

	awk 'NR==FNR{

		hash[$1"\t"$2/3+1]=$3

	}NR!=FNR{

		if(FNR==1){

			$0=$0"\t'$a'"

		}else if(hash[$1"\t"$2]){

			$0=$0"\t"hash[$1"\t"$2]

		}else{

			$0=$0"\t0"

		}

		print

	}' $i  all_sample.txt > all_sample.txt_temp

	mv all_sample.txt_temp all_sample.txt

done


awk 'NR==FNR && NR>1 && $2>1{

	hash1[$1]=hash1[$1]+$4+1
	hash2[$1]=hash2[$1]+$5+1
	hash3[$1]=hash3[$1]+$6+1
	hash4[$1]=hash4[$1]+$7+1
	hash5[$1]=hash5[$1]+$8+1
	hash6[$1]=hash6[$1]+$9+1
	count[$1]++

}NR!=FNR && FNR>1 && $2!=1{

	print $0"\t"($4+1)*count[$1]/hash1[$1]"\t"($5+1)*count[$1]/hash2[$1]"\t"($6+1)*count[$1]/hash3[$1]"\t"($7+1)*count[$1]/hash4[$1]"\t"($8+1)*count[$1]/hash5[$1]"\t"($9+1)*count[$1]/hash6[$1]


}' all_sample.txt all_sample.txt | awk '{print $0"\t"($10+$11+$12)/($13+$14+$15)}' | sort -k16,16rn > all_sample_Rel_occ.txt




awk '{

	if($10>=$11){

		if($12<=$11){

			a=$11
		
		}else if($12>=$11 && $12<=$10){

			a=$12

		}else{

			a=$10
		}

	}else{

		if($12<=$10){

			a=$10

		}else if($12>=$10 && $12<=$11){

			a=$12

		}else{

			a=$11
		}

	}

	if($13>=$14){

		if($15<=$14){

			b=$14

		}else if($15>=$14 && $15<=$13){

			b=$15

		}else{

			b=$13
		}

	}else{

		if($15<=$13){

			b=$13

		}else if($15>=$13 && $15<=$14){

			b=$15

		}else{

			b=$14
		}

	}

	print $0"\t"a/b

}' all_sample_Rel_occ.txt  | sort -k17,17rn    > all_sample_Rel_occ_filter0site_median.txt


# Generate dicodon plot

ExcludedFile=ExcludedFileInCodonPlotResutls

for i in  *R1*_Asite_offset.txt

do 

	echo $i

	a=${i/_*/}

	echo $a

	awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{

		len[$1]=$2-15

		hash[$1"\t"$2]=1

		codon[$1"\t"$2]=$3

	}ARGIND==2{

		a=$2/3+1

		if(hash[$1"\t"a]){

			hash1[$1]=hash1[$1]+$3

			hash2[$1"\t"a]=$3

		}

	}ARGIND==3{

		if(hash1[$1]!="" && hash1[$1]/len[$1]>=0.1){

			print  $0"\t"hash2[$1"\t"$2]*len[$1]/hash1[$1]

		}
	}' $ExcludedFile $i $ExcludedFile | awk 'BEGIN{OFS="\t";FS="\t"}NR==FNR{

		hash1[$1"\t"$2]=$3

		if(hash1[$1"\t"$2-1]){

			class[$1"\t"$2]=hash1[$1"\t"$2-1]"\t"$3

		}

	}NR!=FNR{

		if(class[$1"\t"$2]){

			value[class[$1"\t"$2]]=value[class[$1"\t"$2]]+$4

			count[class[$1"\t"$2]]++

		}}END{

			for(i in value){

				print i"\t"value[i]/count[i]"\t"value[i]"\t"count[i]"\t'$a'"

			}
			
		}' $ExcludedFile   - >> dicodon_paper.txt

done






