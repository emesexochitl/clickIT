setwd("/media/emese/331109B65C43AFD8/projects_2017/philipp_clickit/time_course_tpm/")
library("Cairo", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/Rlibs")
library("reshape2", lib.loc="~/Rlibs")
library("pheatmap", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
# 
# 
heatdata <- read.table("expression_norm1_vals_1.25_average_filtering.txt", sep ="\t", header=FALSE, row.names = 1)
## header?
colnames(heatdata) <- c(0,1,2,6,12,24)
heatdata <- data.matrix(heatdata, rownames.force = NA)
map = pheatmap(heatdata,cluster_rows = TRUE, cluster_cols = F, show_rownames = FALSE)
map = pheatmap(heatdata,cluster_rows = TRUE, cluster_cols = F, kmeans_k = 10)
set.seed(1)
map = pheatmap(heatdata,cluster_rows = TRUE, cluster_cols = F, kmeans_k = 5)
cluster_list <- data.frame(map$kmeans$cluster)

mynum = 1
#mylist <- cluster_list[cluster_list$map.kmeans.cluster== mynum, ,drop=FALSE]
mylist <- read.table(file = paste0("average_chron_cluster", mynum, "_clean.txt"), header = FALSE, row.names = 1)
set.seed(1)
#mylist <- mylist[sample(nrow(mylist), 6300), ,drop=FALSE]
#nrow(mylist)
# 
# #write.table(mylist, file = paste0("average_chron_cluster", mynum, ".txt"), quote = F)
myname <- rownames(mylist)
mymat <- which(rownames(heatdata) %in% myname)
myhits <- heatdata[mymat,]
colnames(myhits) <- c(0,1,2,6,12,24)
nrow(myhits)
hits_melt <- melt(myhits)
hits_melt[ ,c(2)] <- sapply(hits_melt[ ,c(2)], as.numeric)
#ggplot(hits_melt, aes(hits_melt$Var2, hits_melt$value)) + geom_point() + geom_smooth(method = "loess") + theme_classic() + xlab("Time points") + ylab("Values")  + ggtitle(paste0("Cluster ", mynum, ", Loess")) +scale_x_discrete(breaks=c(0:24), labels=c(0:24),limits=c(0:24)) + geom_hline(yintercept = 0.5, color = "red", linetype ="dashed")
#ggsave(paste0("average_chron_cluster", mynum, "_loess.eps"), device = cairo_ps)
p <- ggplot(hits_melt, aes(hits_melt$Var2, hits_melt$value)) + geom_point() + geom_smooth(method = "loess") + theme_classic() + xlab("Time points") + ylab("Values")  + ggtitle(paste0("Cluster ", mynum, ", Loess")) +scale_x_discrete(breaks=c(0:24), labels=c(0:24),limits=c(0:24)) + geom_hline(yintercept = 0.5, color = "red", linetype ="dashed")
v0 <- 0.5
df <- data.frame(x=ggplot_build(p)$data[[2]]["x"])
df <- within(df, y <- ggplot_build(p)$data[[2]]["y"])
df[ ,c(2)] <- sapply(df[ ,c(2)], as.numeric)
f2 <- approxfun(df$y, df$x)
f2(v0)
cross <- f2(v0)
#ggplot(hits_melt, aes(hits_melt$Var2, hits_melt$value)) + geom_point() + geom_smooth(method = "loess") + theme_classic() + xlab("Time points") + ylab("Values")  + ggtitle(paste0("Cluster ", mynum, ", Loess")) +scale_x_discrete(breaks=c(0:24), labels=c(0:24),limits=c(0:24)) + geom_hline(yintercept = 0.5, color = "red", linetype ="dashed") + annotate("text", x = 17, y=1.2, label = paste0("Half-life intercept at ", cross))
#ggsave(paste0("average_chron_cluster", mynum, "_loess.eps"), device = cairo_ps)
postscript(paste0("cluster_exp_lowess_", mynum, ".ps"))
#plot(hits_melt, main = paste0("lowess(cluster ", mynum, ")"))
plot(x=jitter(hits_melt$Var2), y=hits_melt$value, main = paste0("lowess and exponential decay fitting (cluster ", mynum, ")"), cex =0.1, xaxt="n", yaxt="n", bty="l")
lines(lowess(x=hits_melt$Var2, y=hits_melt$value, f=0.2), col = "red", lwd=1)
abline(h=0.5, col="red", lty=2, lwd=1)
fit <- nls(value ~ a*exp(-k * Var2), hits_melt, start=c(a=1, k=1)) 
pr.lm <- predict(fit)
lines(pr.lm~hits_melt$Var2, col="blue", lwd=1)
axis(1, at=c(0,1,2,6,12,24))
axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
mytext <- paste0("Lowess, predicted half-life: ", cross)
text(17, 0.9, mytext, col="red")
text(17, 0.8, "Exponential decay", col="blue") 
dev.off()
#ggplot(hits_melt, aes(hits_melt$Var2, hits_melt$value)) + geom_point() + geom_smooth(method = "auto") + theme_classic() + xlab("Time points") + ylab("Values")  + ggtitle(paste0("Cluster ", mynum, ", Auto"))
