# Fishers exact test

######################################################

obs=18393.29
total = 1000000
noncat= total-obs

bg = 28700.91
bg_total = 1000000
bg_noncat = total-bg
ct <- matrix(c(obs, noncat, bg, bg_noncat), nrow = 2, ncol = 2)

colnames(ct) <- c("y", "n")
rownames(ct) <- c("inv", "rest")
ct
fisher.test(ct, alternative = "two.sided")
fisher.test(ct, alternative = "greater")
fisher.test(ct, alternative = "less")


