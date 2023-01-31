#QC
ls ../01raw_data/*gz|xargs fastqc -t 10 -o  ./

# Quality and adapter trimming:
for f1 in *_1.fastq.gz; do f2=${f1/_1.fastq.gz/_2.fastq.gz} 
1paired=${f1//.fastq/_paired.fastq.gz} 
1unpaired=${f1//.fastq/_unpaired.fastq.gz} 
2paired=${f2//.fastq/_paired.fastq.gz} 
2unpaired=${f2//.fastq/_unpaired.fastq.gz} && 
trimmomatic PE -threads 16 -phred33 -trimlog logfile $f1 $f2 $1paired $1unpaired $2paired $2unpaired 
ILLUMINACLIP:/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:8:true  
SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:36;done

# Build index:
STAR --runThreadN 6 --runMode genomeGenerate \
--genomeDir arab_STAR_genome_149bp_12 \
--genomeFastaFiles /home/cfx756/RNA-Seq/00ref/TAIR10_chr_all.fas \
--sjdbGTFfile /home/cfx756/RNA-Seq/00ref/Araport11_GFF3_genes_transposons.201606.gtf \
--sjdbOverhang 149

# Alignment to tair10
for f1 in /02clean_data/*paired_R1.fastq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && STAR --runThreadN 5 \
--genomeDir arab_STAR_genome_149bp_12 \
--readFilesCommand zcat \
--readFilesIn $f1 $f2  \
--outFileNamePrefix ${f1/paired_R1.fastq.gz/} \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 5 \
--quantMode TranscriptomeSAM GeneCounts \
;done

# directly sort BAM files and filter
for file in *toTranscriptome.out.bam; do echo $file && samtools view -hu -q 10 $file | samtools sort - -o ${file/.bam/_sorted.bam} && rm $file; done
#Rsem prepare reference
rsem-prepare-reference --gtf Araport11_GFF3_genes_transposons.201606.gtf TAIR10_chr_all.fas arab_rsem

#RSEM
for f1 in *Aligned.toTranscriptome.out_filted.bam; 
do echo $f1  
nohup rsem-calculate-expression --paired-end --no-bam-output --alignments -p 15 -q $f1 arab_rsem  ${f1/Aligned.toTranscriptome.out_filted.bam/_rsem} & done 


#build data matrix 
rsem-generate-data-matrix *pcp_rep*_rsem.genes.results > ../08deseq_out/pcp_output.matrix

#delete unexpressed genes
awk 'BEGIN{printf "geneid\ta1\ta2\ta3\tb1\tb2\tb3\n"}{if($2+$3+$4+$5+$6+$7>0)print $0}' pcp_output.matrix |cut -f 1,2,3,4,5,6,7 > pcp_gene.txt

