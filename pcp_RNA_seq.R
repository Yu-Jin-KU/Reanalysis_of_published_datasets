#pcp


input_data <- read.table("pcp_gene.txt", header=TRUE, row.names=1)

input_data <-round(input_data,digits = 0)


input_data <- as.matrix(input_data)
condition <- factor(c(rep("pcp", 3), rep("Control", 3)))
coldata <- data.frame(Sample = factor(c(rep("pcp",3), rep("Control",3)), 
                                      levels = c("Control", "pcp")),
                      row.names = colnames(input_data))
coldata
coldata$Sample


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=input_data,colData=coldata, design=~Sample)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05, independentFiltering=FALSE)
summary(res)
res <- res[order(res$padj),]
res
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized=TRUE)), 
                 by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"


# get diff_gene
diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized=F)), 
                 by="row.names", sort=FALSE)

head(resdata)
resdata$significant <- "unchanged_pcp"
resdata$significant[resdata$padj < 0.05 & resdata$log2FoldChange > 1 ] <- "upregulated_pcp"
resdata$significant[resdata$padj < 0.05 & resdata$log2FoldChange < -1 ] <- "downregulated_pcp"

write.table(resdata,file="pcp_diffexpr.txt",sep = "\t",quote = F, row.names = F)
