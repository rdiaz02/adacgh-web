largesamp <- read.table("large.sample", header = TRUE, comment.char = "", sep = "\t")


m1 <- largesamp[sort(sample(1:42325, size = 5000, replace = FALSE)), ]
write.table(m1, file = "m1", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

m2 <- largesamp[sort(sample(1:42325, size = 10000, replace = FALSE)), ]
write.table(m1, file = "m2", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
