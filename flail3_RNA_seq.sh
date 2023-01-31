#QC
ls ../01raw_data/*gz|xargs fastqc -t 10 -o  ./

# Quality and adapter trimming:
for f1 in *R1*; do f2=${f1//R1.fastq.gz/R2.fastq.gz} 
R1paired=${f1//.fastq/_paired.fastq.gz} 
R1unpaired=${f1//.fastq/_unpaired.fastq.gz} 
R2paired=${f2//.fastq/_paired.fastq.gz} 
R2unpaired=${f2//.fastq/_unpaired.fastq.gz} && 
trimmomatic PE -threads 16 -phred33 -trimlog logfile $f1 $f2 $R1paired $R1unpaired $R2paired $R2unpaired 
ILLUMINACLIP:/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:8:true  
SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:36;done

# Build index:
STAR --runThreadN 6 --runMode genomeGenerate 
--genomeDir arab_STAR_genome \ 
--genomeFastaFiles TAIR10_chr_all.fas \
--sjdbGTFfile Araport11_GFF3_genes_transposons.201606.gtf \
--sjdbOverhang 74

# Alignment to tair10
for f1 in /02clean_data/*_paired_R1.fastq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && nohup STAR \
--runThreadN 5 \
--genomeDir arab_STAR_genome \
--readFilesCommand zcat \
--readFilesIn $f1 $f2  \
--outFileNamePrefix ${f1/paired_R1.fastq.gz/} \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 5 \
--quantMode TranscriptomeSAM GeneCounts & done

# filter BAM files 
for file in *toTranscriptome.out.bam; do echo $file && samtools view -hu -q 10 $file -o ${file/.bam/_filted.bam}; done

#Rsem prepare reference
rsem-prepare-reference --gtf Araport11_GFF3_genes_transposons.201606.gtf TAIR10_chr_all.fas arab_RSEM/arab_rsem

# Make unstranded Bedgraph files (for visualization in genomic browsers):
for file in *Aligned.sortedByCoord.out.bam; do echo $file && bedtools genomecov -ibam $file -bg -split | sort -k1,1 -k2,2n | sed "1i track type=bedGraph" | gzip > ${file/bam/bedgraph.gz}; done

#RSEM
for f1 in *Aligned.toTranscriptome.out_filted.bam; 
do echo $f1  
nohup rsem-calculate-expression --paired-end --no-bam-output --alignments -p 15 -q $f1 arab_rsem  ${f1/Aligned.toTranscriptome.out_filted.bam/_rsem} & done 


#build data matrix 
rsem-generate-data-matrix *_rsem.genes.results > ../output.matrix


#delete unexpressed genes
awk 'BEGIN{printf "geneid\ta1\ta2\ta3\tb1\tb2\tb3\tc1\tc2\tc3\n"}{if($2+$3+$4+$5+$6+$7+$8+$9+$10>0)print $0}' output.matrix |cut -f 1,5,6,7,8,9,10 > kd_gene.txt
awk 'BEGIN{printf "geneid\ta1\ta2\ta3\tb1\tb2\tb3\tc1\tc2\tc3\n"}{if($2+$3+$4+$5+$6+$7+$8+$9+$10>0)print $0}' output.matrix |cut -f 1,2,3,4,5,6,7,8,9,10 > FLAIL_KD_OE_gene.txt
