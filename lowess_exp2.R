setwd("/media/emese/331109B65C43AFD8/projects_2017/philipp_clickit/time_course_tpm/")
library("Cairo", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/Rlibs")
library("reshape2", lib.loc="~/Rlibs")
library("pheatmap", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")

mynum = 1

mylist <- read.table(file = paste0("average_chron_cluster", mynum, "_clean.txt"), header = FALSE, row.names = 1)
#set.seed(1)
#mylist <- mylist[sample(nrow(mylist), 6300), ,drop=FALSE]
nrow(mylist)
myname <- rownames(mylist)
mymat <- which(rownames(heatdata) %in% myname)
myhits <- heatdata[mymat,]
colnames(myhits) <- c(0,1,2,6,12,24)
nrow(myhits)
hits_melt <- melt(myhits)
hits_melt[ ,c(2)] <- sapply(hits_melt[ ,c(2)], as.numeric)
v0 <- 0.5
lowess_pred <- lowess(x=hits_melt$Var2, y=hits_melt$value, f=0.1)
lowess_df <- data.frame(x=lowess_pred$x)
lowess_df <- within(lowess_df, y <- lowess_pred$y)
lowess_f2 <- approxfun(lowess_df$y, lowess_df$x)
lowess_f2(v0)
lowess_cross <- lowess_f2(v0)

fit <- nls(value ~ a*exp(-k * Var2), hits_melt, start=c(a=1, k=1)) 
pr.lm <- predict(fit)
exp_df <- data.frame(x=hits_melt$Var2)
exp_df <- within(exp_df, y <- pr.lm)
exp_f2 <- approxfun(exp_df$y, exp_df$x)
exp_f2(v0)
exp_cross <- exp_f2(v0)

#postscript(paste0("cluster_exp_lowess_", mynum, ".ps"))

plot(x=jitter(hits_melt$Var2), y=hits_melt$value, main = paste0("lowess and exponential decay fitting (cluster ", mynum, ")"), cex =0.1, xaxt="n", yaxt="n", bty="l")
lines(lowess(x=hits_melt$Var2, y=hits_melt$value, f=0.2), col = "red", lwd=1)
abline(h=0.5, col="red", lty=2, lwd=1)
#fit <- nls(value ~ a*exp(-k * Var2), hits_melt, start=c(a=1, k=1)) 
#pr.lm <- predict(fit)
lines(pr.lm~hits_melt$Var2, col="blue", lwd=1)
axis(1, at=c(0,1,2,6,12,24))
axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
lowess_mytext <- paste0("Lowess, predicted half-life: ", lowess_cross)
text(17, 0.9, lowess_mytext, col="red")
exp_mytext <- paste0("Exponential decay, predicted half-life: ", exp_cross)
text(17, 0.8, exp_mytext, col="blue") 
#dev.off()
