par(mfrow = c(1, 3))
DNAcopyDiagnosticPlots(CNA.object,
                       smoothed.CNA.object)

par(mfrow = c(3, 3))
par(ask = TRUE)
DNAcopyDiagnosticPlots(CNA.object,
                       smoothed.CNA.object, array.chrom = TRUE,
                       chrom.numeric = chrom.numeric)

#####################

data(cghMCRe)

chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

par(mfrow = c(5, 5))
WaveletsDiagnosticPlots(cghMCRe[, 5:7], chrom.numeric)


################

mpiInit()
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- pSegmentDNAcopy(smoothed.CNA.object)


###******************

data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object)
mergedOut <- mergeDNAcopy(segment.smoothed.CNA.object)



###*********************
data(cghMCRe)
common <- cghMCRe[, -c(5:7)]
common$MidPoint <- common$Start + 0.5 * (common$End - common$Start)
colnames(common)[1] <- "ID"
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

## We use nIter = 100 for the example for speed reasons;
## you probably want 1000 or more.

psw.out <- pSegmentPSW(common, as.matrix(cghMCRe[, 5:7]), chrom.numeric,
                       sign = - 1, nIter = 100, prec = 100, p.crit = 0.10)

## so you don't want all the stuff in common? just use a smaller common.dat
common2 <- cghMCRe[, -c(5:7)]
psw2.out <- pSegmentPSW(common2, as.matrix(cghMCRe[, 5:7]), chrom.numeric,
                       sign = - 1, nIter = 100, prec = 100, p.crit = 0.10)


            
CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object)
mergedOut <- mergeDNAcopy(segment.smoothed.CNA.object)



###***************
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
wave.out <- pSegmentWavelets(cghMCRe[, 5:7], chrom.numeric)

#############
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
ace.out <- pSegmentACE(cghMCRe[, 5:7], chrom.numeric)

