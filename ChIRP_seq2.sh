#QC:
ls ../01raw_data/*gz|xargs fastqc -t 10 -o  ./

#trim_galore for pair end:
for f1 in *R1.fastq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && trim_galore --illumina --paired --trim1 --no_report_file $f1 $f2; done
for file in *val_1*; do mv $file ${file/001_val_1/trimmed}; done
for file in *val_2*; do mv $file ${file/001_val_2/trimmed}; done

# Align to TAIR10 in PE mode:
for f1 in *R1_trimmed.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && nohup bowtie2 -p 10 -x TAIR10_index -1 $f1 -2 $f2 -S ${f1/R1_trimmed.fq.gz/.sam} & done
for file in *sam; do echo $file && samtools view -hq 5 -f 2 -b -@10 -u $file | samtools sort - -o ${file/.sam/_sorted.bam}; done


#picard remove duplicates
ls  *.bam  |xargs -i samtools index {} 
for f1 in *sorted.bam; do echo $f1 && nohup picard MarkDuplicates --INPUT $f1  --OUTPUT ${f1/_sorted.bam/_rmdup.bam} --METRICS_FILE ${f1/_sorted.bam/_rmdup}.log --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true & done

# Make bw files for visualization:
ls  *rmdup.bam  |xargs -i samtools index {} 
ls *rmdup.bam |while read id;do
nohup bamCoverage  --normalizeUsing RPKM -b $id -o ${id%%.*}.rmdup.bw & done

# Convert BAM files to BED:
for file in *rmdup.bam; do bedtools bamtobed -i $file > ${file/bam/bed}; done

# Pool replicates:
cat Luc1*bed Luc2*bed| sort -k1,1 -k2,2n > Luc1_Luc2_merged.bed
cat Luc1*bed Luc3*bed| sort -k1,1 -k2,2n > Luc1_Luc3_merged.bed
cat Luc2*bed Luc3*bed| sort -k1,1 -k2,2n > Luc2_Luc3_merged.bed
cat Even*bed Odd*bed |sort -k1,1 -k2,2n > Even_Odd_merged.bed

# Run CCAT for peakcalling:
CCAT Even_Odd_merged.bed Luc1_Luc2_merged.bed TAIR10_chrom_sizes.txt ./example/config_histone.txt Even_Odd_mergedvsluc1_luc2_merged
CCAT Even_Odd_merged.bed Luc1_Luc3_merged.bed TAIR10_chrom_sizes.txt ./example/config_histone.txt Even_Odd_mergedvsluc1_luc3_merged
CCAT Even_Odd_merged.bed Luc2_Luc3_merged.bed TAIR10_chrom_sizes.txt ./example/config_histone.txt Even_Odd_mergedvsluc2_luc3_merged

# Convert output files to BED:
for file in *significant.peak; do echo $file && cut -f 1,3,4,5,6,7,8 $file | sort -k1,1 -k2,2n | sed 's/ChrM/Mt/;s/ChrC/Pt/;s/Chr//'> ${file}.bed; done
