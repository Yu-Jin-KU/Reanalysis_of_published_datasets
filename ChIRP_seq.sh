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
