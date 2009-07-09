rm(list = ls())
library(ADaCGH)

dir.to.create <- tempdir2()
setwd(dir.to.create)
dir.to.create

library(Rmpi)
mpiInit(universeSize = mpi.universe.size())


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

cbs.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric)
 haar.out <- pSegmentHaarSeg(cghE1[, 5:7], chrom.numeric)

yminmax <- c(min(as.matrix(cghE1[, 5:7])),
                  max(as.matrix(cghE1[, 5:7])))

 segmentPlot(haar.out,
                 geneNames = cghE1[, 1],
                 yminmax = yminmax,
                 idtype = "ug",
                 organism = "Hs")

     segmentPlot(cbs.out, 
                 geneNames = cghE1[, 1],
                 yminmax = yminmax,
                 idtype = "ug",
                 organism = "Hs")
