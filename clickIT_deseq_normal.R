# Import DESeq (without replicates it is necessary to use this).
#detach("package:DESeq2", unload=TRUE)
#detach("package:DESeq", unload=TRUE)
library( "DESeq" )
mutant <- c("Click")
mutant_id <- c("PulseClick_1, PulseClick_2, PulseClick_3")

countdata <- read.table("counts_all_totalclick.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)] 
colnames(countdata) <- gsub("\\_sorted.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("^...genome_alignment_filtered_trt.", "", colnames(countdata))
#coldata <- read.table("coldata.txt", header=T, row.names=1)

cds = newCountDataSet(countData = countdata, conditions = as.factor(c("Click", "Click", "Click", "Total", "Total", "Total")))
cds <- cds[rowSums(counts(cds))>1, ] 
cds = estimateSizeFactors(cds)
sizeFactors(cds)

#cds2 = cds[ ,c("SRR1104133_salt", mutant_id)]
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds, "Total", mutant)
head(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
plotMA(res[order(res$padj),])
res3 <-res[is.finite(res$foldChange) & is.finite(res$log2FoldChange), ]
max(res3$log2FoldChange)
min(res3$log2FoldChange)
plotMA(res3[order(res3$padj),], main=paste0("DESeq2 control (Mock) and ", mutant, " samples"), ylim = c(1.2*min(res3$log2FoldChange),1.2*max(res3$log2FoldChange)))

png(filename=paste0("maplot_Total_normal_", mutant, "_FDR_01.png"))
plotMA(res3[order(res3$padj),], main=paste0("DESeq2 control (Total) and ", mutant, " samples, padj < 0.1"), ylim = c(1.2*min(res3$log2FoldChange),1.2*max(res3$log2FoldChange)))
dev.off()

png(filename=paste0("hist_Total_normal_", mutant, "_", type, ".png"))
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

png(filename=paste0("disp_Total_normal_", mutant, "_", type, ".png"))
plotDispEsts(cds)
dev.off()


write.table(res[order(res$padj), ], file = paste0("deseq_Total_normal", mutant, "_all.txt"), quote =F)
#write.table((res2[which(res2$padj < 0.05), ]), file = paste0("deseq_Mock_", mutant, "_005.txt"), quote=F) 
#write.table((res2[which(res2$padj < 0.05 & res2$log2FoldChange> 1), ]), file = paste0("deseq_Mock_", mutant, "_fc1up_005.txt"), quote=F)
#write.table((res2[which(res2$padj < 0.1 & res2$log2FoldChange> 1), ]), file = paste0("deseq_Mock_", mutant, "_fc1up_01.txt"), quote=F)
