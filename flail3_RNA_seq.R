#kD(flail)
input_data <- read.table("kd_gene.txt", header=TRUE, row.names=1)

input_data <-round(input_data,digits = 0)


input_data <- as.matrix(input_data)
condition <- factor(c(rep("KD", 3), rep("Control", 3)))
coldata <- data.frame(Sample = factor(c(rep("KD",3), rep("Control",3)), 
                                      levels = c("Control", "KD")),
                      row.names = colnames(input_data))
coldata
coldata$Sample


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=input_data,colData=coldata, design=~Sample)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05)
summary(res)
res <- res[order(res$padj),]
res
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized=TRUE)), 
                 by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"


# get diff_gene
diff_gene <- subset(res, padj < 0.05 & abs(resdata$log2FoldChange) >= 1)

resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized=F)), 
                 by="row.names", sort=FALSE)

head(resdata)
resdata$significant <- "unchanged_KD"
resdata$significant[resdata$padj <= 0.05 & resdata$log2FoldChange >= 1 ] <- "upregulated_KD"
resdata$significant[resdata$padj <= 0.05 & resdata$log2FoldChange <= -1 ] <- "downregulated_KD"

write.table(resdata,file="flail_diffexpr.txt",sep = "\t",quote = F, row.names = F)


# Complementation_gene
input_data <- read.table("Complementation_gene.txt", header=TRUE, row.names=1)

input_data <-round(input_data,digits = 0)


input_data <- as.matrix(input_data)
condition <- factor(c(rep("OE", 2), rep("Control", 3)))
coldata <- data.frame(Sample = factor(c(rep("OE",2), rep("Control",3)), 
                                      levels = c("Control", "OE")),
                      row.names = colnames(input_data))
coldata
coldata$Sample


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=input_data,colData=coldata, design=~Sample)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05,independentFiltering=FALSE)
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

resdata$significant <- "unchanged_OE"
resdata$significant[resdata$padj < 0.05 & resdata$log2FoldChange > 1 ] <- "upregulated_OE"
resdata$significant[resdata$padj < 0.05 & resdata$log2FoldChange < -1 ] <- "downregulated_OE"

head(resdata)
write.table(resdata,file="Complementation_diffexpr.txt",sep = "\t",quote = F, row.names = F)

