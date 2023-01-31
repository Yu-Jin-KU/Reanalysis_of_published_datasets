#QC
ls ../01raw_data/*gz|xargs fastqc -t 10 -o  ./

# Quality and adapter trimming
ls raw/*fastq | while read f1; do nohup trim_galore --illumina -o /home/cfx756/flail/NSRa_RIPseq/clean $f1 &  done

# Alignment to tair10
for f1 in *fq; do nohup STAR \
--runThreadN 5 \
--genomeDir arab_STAR_genome_74bp_12 \
--readFilesCommand zcat \
--readFilesIn $f1   \
--outFileNamePrefix ${f1/_trimmed.fq/} \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 5 \
--quantMode TranscriptomeSAM GeneCounts & done


# sam to bam and filter bam files
for file in *toTranscriptome.out.bam; do echo $file && samtools view -hu -q 10 $file | samtools sort - -o ${file/.bam/_sorted.bam} && rm $file; done

#Rsem prepare reference
rsem-prepare-reference --gtf Araport11_GFF3_genes_transposons.201606.gtf TAIR10_chr_all.fas arab_rsem

#RSEM
for f1 in *Aligned.toTranscriptome.out_filted.bam; 
do echo $f1  
nohup rsem-calculate-expression --paired-end --no-bam-output --alignments -p 15 -q $f1 arab_rsem  ${f1/Aligned.toTranscriptome.out_filted.bam/_rsem} & done 


#build data matrix 
rsem-generate-data-matrix *_rsem.genes.results > ../08deseq_out/nsrab_RIP_output.matrix
