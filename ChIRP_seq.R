library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library (dplyr)

Even_Odd_mergedvsluc1_luc2_merge <- readPeakFile( peakfile ="Even_Odd_mergedvsluc1_luc2_merged.significant.peak.bed", header= FALSE ) 
Even_Odd_mergedvsluc1_luc3_merge <- readPeakFile( peakfile ="Even_Odd_mergedvsluc1_luc3_merged.significant.peak.bed", header= FALSE ) 
Even_Odd_mergedvsluc2_luc3_merge <- readPeakFile( peakfile ="Even_Odd_mergedvsluc2_luc3_merged.significant.peak.bed", header= FALSE ) 


peakAnno17 <- annotatePeak(Even_Odd_mergedvsluc1_luc2_merge, tssRegion=c(-300, 300), 
                           TxDb=txdb)
peakAnno18 <- annotatePeak(Even_Odd_mergedvsluc1_luc3_merge, tssRegion=c(-300, 300), 
                           TxDb=txdb)
peakAnno19 <- annotatePeak(Even_Odd_mergedvsluc2_luc3_merge, tssRegion=c(-300, 300), 
                           TxDb=txdb)

annotation17 <- as.data.frame(peakAnno17)
annotation18 <- as.data.frame(peakAnno18)
annotation19 <- as.data.frame(peakAnno19)

merge15 <- full_join(annotation17, annotation18, by =c("geneId" = "geneId"))
merge16 <- full_join(merge15, annotation19, by =c("geneId" = "geneId"))
write.table(merge16,file="Even_Odd_mergedvsluc2_luc2_or_luc3_merge.txt",sep = "\t",quote = F, row.names = F)
