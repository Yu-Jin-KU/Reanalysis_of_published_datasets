#QC:
ls ../01raw_data/*gz|xargs fastqc -t 10 -o  ./

#trim_galore for pair end:
for f1 in *R1.fastq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && trim_galore --illumina --paired --trim1 --no_report_file $f1 $f2; done
for file in *val_1*; do mv $file ${file/001_val_1/trimmed}; done
for file in *val_2*; do mv $file ${file/001_val_2/trimmed}; done

# Align to TAIR10 in PE mode:
for f1 in *R1_trimmed.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && nohup STAR --genomeDir arab_STAR_genome --readFilesIn $f1 $f2 --runThreadN 4 --outFileNamePrefix ${f1/R1_trimmed.fq.gz/} --outSAMmultNmax 1 --alignEndsType EndToEnd --readFilesCommand zcat --outSAMtype BAM Unsorted --alignIntronMax 1 --alignMatesGapMax 1000 & done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Remove low MAPQ reads and reads not in proper pair:
for file in *bam; do echo $file && samtools view -hq 5 -f 2 $file | samtools sort -n - -o ${file/.bam/_mapq.bam}; done

# Run fixmate (required by markdup):
for file in *mapq.bam; do echo $file && samtools fixmate -m $file ${file/.bam/_fixmate.bam}; done

# Sort by coordinates:
for file in *fixmate.bam; do echo $file && samtools sort $file -o ${file/.bam/_sorted.bam}; done

# Deduplicate in PE mode:
for file in *sorted.bam; do echo $file && samtools markdup -r $file ${file/.bam/_dedup.bam}; done


# Make bw files for visualization:
ls  *dedup.bam  |xargs -i samtools index {} 
ls *dedup.bam |while read id;do
nohup bamCoverage  --normalizeUsing RPKM -b $id -o ${id%%.*}.dedup.bw & done

# Convert BAM files to BED:
for file in *dedup.bam; do bedtools bamtobed -i $file > ${file/bam/bed}; done

# Pool replicates:
cat Luc1*bed Luc2*bed| sort -k1,1 -k2,2n > Luc1_Luc2_merged.bed
cat Luc1*bed Luc3*bed| sort -k1,1 -k2,2n > Luc1_Luc3_merged.bed
cat Luc2*bed Luc3*bed| sort -k1,1 -k2,2n > Luc2_Luc3_merged.bed
cat Even*bed Odd*bed |sort -k1,1 -k2,2n > Even_Odd_merged.bed

# Run CCAT for peakcalling:
CCAT Even_Odd_merged.bed Luc1_Luc2_merged.bed TAIR10_chrom_sizes.txt ./example/config_histone.txt Even_Odd_vs_Luc1_Luc2_merged
CCAT Even_Odd_merged.bed Luc1_Luc3_merged.bed TAIR10_chrom_sizes.txt ./example/config_histone.txt Even_Odd_vs_Luc1_Luc3_merged
CCAT Even_Odd_merged.bed Luc2_Luc3_merged.bed TAIR10_chrom_sizes.txt ./example/config_histone.txt Even_Odd_vs_Luc2_Luc3_merged

# Convert output files to BED:
for file in *significant.peak; do echo $file && cut -f 1,3,4,5,6,7,8 $file | sort -k1,1 -k2,2n | sed 's/ChrM/Mt/;s/ChrC/Pt/;s/Chr//'> ${file}.bed; done
