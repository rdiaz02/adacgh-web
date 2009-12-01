library(ADaCGH2)

data(cghE1)
tmpchr <- sub("chr", "", cghE1$Chromosome)
chrom.numeric <- as.numeric(as.character(tmpchr))
chrom.numeric[tmpchr == "X"] <- 23
chrom.numeric[tmpchr == "Y"] <- 24
rm(tmpchr)
### we need the data ordered
reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]


dir.to.create <- tempdir2()
setwd(dir.to.create)
dir.to.create

convertAndSave(cghE1$Name, chrom.numeric, cghE1$UG.Start, cghE1[, 5:7])

pSegmentHaarSeg("cghData.RData", "chromData.RData")
