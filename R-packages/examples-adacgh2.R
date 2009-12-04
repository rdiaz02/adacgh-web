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


haar.out <- pSegmentHaarSeg("cghData.RData", "chromData.RData")
open(haar.out[[1]]); open(haar.out[[2]])
summary(haar.out[[1]][,])
summary(haar.out[[2]][,])


glad.out <- pSegmentGLAD("cghData.RData", "chromData.RData")
open(glad.out[[1]]); open(glad.out[[2]])
summary(glad.out[[1]][,])
summary(glad.out[[2]][,])

cbs.out <- pSegmentDNAcopy("cghData.RData", "chromData.RData")
open(cbs.out[[1]]); open(cbs.out[[2]])
summary(cbs.out[[1]][,])
summary(cbs.out[[2]][,])

hmm.out <- pSegmentHMM("cghData.RData", "chromData.RData")
open(hmm.out[[1]]); open(hmm.out[[2]])
summary(hmm.out[[1]][,])
summary(hmm.out[[2]][,])


waves.out <- pSegmentWavelets("cghData.RData", "chromData.RData")
open(waves.out[[1]]); open(waves.out[[2]])
summary(waves.out[[1]][,])
summary(waves.out[[2]][,])

waves.out.ml <- pSegmentWavelets("cghData.RData", "chromData.RData", merging = "mergeLevels")
open(waves.out.ml[[1]]); open(waves.out.ml[[2]])
summary(waves.out.ml[[1]][,])
summary(waves.out.ml[[2]][,])

waves.out.nm <- pSegmentWavelets("cghData.RData", "chromData.RData", merging = "none")
open(waves.out.nm[[1]]); open(waves.out.nm[[2]])
summary(waves.out.nm[[1]][,])
summary(waves.out.nm[[2]][,])


cghseg.out <- pSegmentCGHseg("cghData.RData", "chromData.RData")
open(cghseg.out[[1]]); open(cghseg.out[[2]])
summary(cghseg.out[[1]][,])
summary(cghseg.out[[2]][,])

cghseg.out.ml <- pSegmentCGHseg("cghData.RData", "chromData.RData", merging = "mergeLevels")
open(cghseg.out.ml[[1]]); open(cghseg.out.ml[[2]])
summary(cghseg.out.ml[[1]][,])
summary(cghseg.out.ml[[2]][,])

cghseg.out.nm <- pSegmentCGHseg("cghData.RData", "chromData.RData", merging = "none")
open(cghseg.out.nm[[1]]); open(cghseg.out.nm[[2]])
summary(cghseg.out.nm[[1]][,])
summary(cghseg.out.nm[[2]][,])



biohmm.out <- pSegmentBioHMM("cghData.RData", "chromData.RData", "posData.RData")
open(biohmm.out[[1]]); open(biohmm.out[[2]])
summary(biohmm.out[[1]][,])
summary(biohmm.out[[2]][,])
