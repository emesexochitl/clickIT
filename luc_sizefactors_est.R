counts(dds, normalized=TRUE)setwd("/media/emese/331109B65C43AFD8/projects_2017/philipp_clickit/time_course_de/")
countdata <- read.table("counts_time_course_only_luc.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)] 
coldata <- read.table("coldata.txt", header=T, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~time)
dds <- estimateSizeFactors(dds)
luc_sizefactors <- sizeFactors(dds)
