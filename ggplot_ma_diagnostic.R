sizeres2 = res
res2$significant = (res2$padj < .05)
res2$significant = as.factor(res2$significant)
res2$significant[is.na(res2$significant)] = F



ggplot(as.data.table(res2), aes(x=log2(baseMean), y=log2FoldChange, color=significant)) + geom_point() + geom_hline(color = "blue3", yintercept = 0) + stat_smooth(se = FALSE, method = "auto", color = "red3") + scale_color_manual(values=c("Black","Red"))

