click_tpm_no0 <- read.table("gene_click_all_TPM_vals.txt", header=FALSE)
total_tpm_no0 <- read.table("gene_total_all_TPM_vals.txt", header=FALSE)
total_tpm_no0$category <- paste0("total_",total_tpm_no0$V6)
click_tpm_no0$category <- paste0("click_",click_tpm_no0$V6)
tpm_table_no0 <- rbind(click_tpm_no0, total_tpm_no0, stringsAsFactors=FALSE)
ggplot(tpm_table_no0, aes(x=V4)) + geom_density(aes(group=category, color=category)) + theme_classic() + scale_x_log10(expand=c(0,0))
+ ggtitle("All TPM density") + scale_y_continuous(expand = c(0, 0))

tpm_table_no0_minus <- tpm_table_no0[grep("minus", tpm_table_no0$category), ]
ggplot(tpm_table_no0_minus, aes(x=V4)) + geom_density(aes(group=category, color=category)) + theme_classic() + scale_x_log10(expand=c(0,0)) + ggtitle("All TPM density") + scale_y_continuous(expand = c(0, 0))
ggplot(intron_log2fc, aes(intron_log2fc$`Click log2TPM`, intron_log2fc$`Total log2 TPM`, color= intron_log2fc$Category)) + geom_point(shape=1) + geom_abline(intercept = 0, slope = 1,color = "red") +  ggtitle("log2TPMs of intronic features of 1st isoforms") +theme(plot.title = element_text(hjust = 0.5)) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank()) + xlab(label = "nascent RNA log2TPM of first isoform introns") + ylab(label = "total RNA log2TPM of first isoform introns") + theme(legend.title = element_text("blank")) + scale_color_manual(values=c('tomato2', 'lightgrey', 'darkgrey', 'lightblue2')) + theme(text=element_text(family="Arial", size=8)) + scale_color_manual(labels = c("total RNA reads only", "no reads", "normal", "nascent RNA reads only"), values = c('tomato2', 'lightgrey', 'darkgrey', 'lightblue2'))
