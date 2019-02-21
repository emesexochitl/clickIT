# Import and edit read count file.
countdata <- read.table("counts_all_totalclick.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)] 
colnames(countdata) <- gsub("\\_sorted.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("^...genome_alignment_filtered_trt.", "", colnames(countdata))

# Set up DESeq2 and run it.
mutant <- c("Click")
type <- c("exon")

coldata <- read.table("coldata.txt", header=T, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds))>1, ] 
dds$condition <- relevel(dds$condition, ref="Total") # for double-check.

#dispersions(dds) <- 0.1
dds <- DESeq(dds) 
res <- results(dds) 
res <- results(dds, contrast=c("condition", mutant, "Total")) # Note that the control is always the third value!
res 
summary(res) 
summary(res, alpha = 0.05)

resOrdered <- res[order(res$padj), ]
resSig01 <- res[which(res$padj < 0.1), ]
resSig005 <- res[which(res$padj < 0.05), ]
sum(res$padj < 0.05, na.rm = T)
head(res[which(res$padj < 0.1), ])

# Visualization.
#rld <- rlog(dds, blind = F) 
#data <- plotPCA(rld, intgroup="condition", returnData =T) 
#percentVar <- round(100*attr(data, "percentVar")) 
#ggplot(data, aes(PC1, PC2, color=group, shape=condition)) + geom_point(size=9) + geom_text(aes(label=condition), size =3, vjust=-1.75,colour = "black") + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + ggtitle("PCA plot of RNA-Seq samples (ClickIT)") +theme(legend.position="none", plot.title = element_text(hjust = 0.5))

plotMA(res[order(res$padj),], main=paste0("DESeq2 control (Total) and ", mutant, " samples, padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
plotMA(res[order(res$padj),], main=paste0("DESeq2 control (Total) and ", mutant, " samples, padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
hist(res$pvalue, breaks = 20, main = paste0("Histogram of pvals of DESeq2\ncontrol (Total) and ", mutant, " samples"), col="lightskyblue1", border = F, xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))

#png(filename=paste0("maplot_Total_", mutant, "_", type, "_01.png"))
#plotMA(res[order(res$padj),], main=paste0("DESeq2 control (Total) and ", mutant, " samples, padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
#dev.off()

#png(filename=paste0("maplot_Total_", mutant, "_", type, "_005.png"))
#plotMA(res[order(res$padj),], main=paste0("DESeq2 control (Total) and ", mutant, " samples, padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
#dev.off()

#png(filename=paste0("hist_Total_", mutant, "_", type, ".png"))
#hist(res$pvalue, breaks = 20, main = paste0("Histogram of pvals of DESeq2\ncontrol (Total) and ", mutant, " samples"), col="lightskyblue1", border = F, xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))
#dev.off()

# Write significant results to txt files.
#write.table((res[which(res$padj < 0.05), ]), file = paste0("deseq2_Total_", mutant, "_", type, "_005.txt"), quote=F) 
#write.table((res[which(res$padj < 0.1), ]), file = paste0("deseq2_Total_", mutant, "_", type, "_01.txt"), quote=F)
write.table(res[order(res$padj), ], file = paste0("deseq2_Total_", mutant, "_", type, "_all_nonoverlap.txt"), quote =F)
#write.table((res[which(res$padj < 0.05 & res$log2FoldChange> 0.59 | res$log2FoldChange < -0.59), ]), file = paste0("deseq2_Total_", mutant, "_", type, "_fc1_005.txt"), quote=F)
#write.table((res[which(res$padj < 0.01 & res$log2FoldChange> 0.59 | res$log2FoldChange < -0.59), ]), file = paste0("deseq2_Total_", mutant, "_", type, "_fc1_01.txt"), quote=F) 