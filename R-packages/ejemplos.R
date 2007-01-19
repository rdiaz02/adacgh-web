
#####################


################
setwd("/tmp/o3")
mpiInit()
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

## Recall: we must reorder the data by chromosome and
## by position within chromosome
reorder <- order(chrom.numeric,
                 cghMCRe$Start,
                 cghMCRe$End,
                 cghMCRe$Name)
cghMCRe <- cghMCRe[reorder, ]
chrom.numeric <- chrom.numeric[reorder]
CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- pSegmentDNAcopy(smoothed.CNA.object)

segmentPlot(segment.smoothed.CNA.object, arraynames = colnames(cghMCRe[, 5:7]),
            chrom.numeric = chrom.numeric, idtype = "ug", organism = "Hs",
            geneNames = rownames(cghMCRe),
            yminmax = c(min(as.matrix(cghMCRe[, 5:7])),
            max(as.matrix(cghMCRe[, 5:7]))),
            superimposed = FALSE)





mpiInit()
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[chrom.numeric == "X"] <- 23
chrom.numeric[chrom.numeric == "Y"] <- 24

## Recall: we must reorder the data by chromosome and
## by position within chromosome
reorder <- order(chrom.numeric,
                 cghMCRe$Start,
                 cghMCRe$End,
                 cghMCRe$Name)
cghMCRe <- cghMCRe[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

## chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
## chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
## chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

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
reorder <- order(chrom.numeric,
                 common$MidPoint,
                 cghMCRe$Start,
                 cghMCRe$End,
                 cghMCRe$Name)
cghMCRe <- cghMCRe[reorder, ]
chrom.numeric <- chrom.numeric[reorder]



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

o1 <- pSegmentDNAcopy(smoothed.CNA.object)
o2 <- pSegmentDNAcopy(smoothed.CNA.object,
                      merge = FALSE, nperm = 10000)

h1 <- hmmWrapper(cghMCRe[, 5],
                 Clone = 1:nrow(cghMCRe),
                 Chrom = chrom.numeric,
                 Pos = 1:nrow(cghMCRe))






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
summary(ace.out)


##################
data(cghE1)
tmpchr <- sub("chr", "", cghE1$Chromosome)
chrom.numeric <- as.numeric(as.character(tmpchr))
chrom.numeric[tmpchr == "X"] <- 23
chrom.numeric[tmpchr == "Y"] <- 24
rm(tmpchr)

## Recall: we must reorder the data by chromosome and by position within
## chromosome


reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

## ACE
ace.out <- pSegmentACE(cghE1[, 5:7], chrom.numeric)
ace.out.sum <- summary(ace.out)
segmentPlot(ace.out.sum, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(ace.out.sum, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

## DNA copy + merging
CNA.object <- CNA(as.matrix(cghE1[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghE1),
                  data.type = "logratio",
                  sampleid = colnames(cghE1[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
dnacopy.out <- segment(smoothed.CNA.object)
merged.out <- mergeDNAcopy(dnacopy.out)

segmentPlot(merged.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(merged.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

segmentPlot(dnacopy.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(dnacopy.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

## PSW
psw.loss.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]),
                        chrom.numeric,
                        sign = -1, nIter = 5000, prec = 100,
                        p.crit = 0.10)
psw.gain.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]),
                        chrom.numeric,
                        sign = 1, nIter = 5000, prec = 100,
                        p.crit = 0.10)

segmentPlot(psw.loss.out, geneNames = cghE1[, 1],
            idtype = "ug",
            organism = "Hs")
segmentPlot(psw.gain.out, geneNames = cghE1[, 1],
            idtype = "ug",
            organism = "Hs")


## wavelets
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)
segmentPlot(wave.out, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")
segmentPlot(wave.out, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

## acceptedIDTypes = ('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl', 'entrez', 'ug')
## acceptedOrganisms = ('None', 'Hs', 'Mm', 'Rn')









#################
data(cghE1)
tmpchr <- sub("chr", "", cghE1$Chromosome)
chrom.numeric <- as.numeric(as.character(tmpchr))
chrom.numeric[tmpchr == "X"] <- 23
chrom.numeric[tmpchr == "Y"] <- 24
rm(tmpchr)

## Recall: we must reorder the data by chromosome and by position within
## chromosome


reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

### DNA copy
CNA.object <- CNA(as.matrix(cghE1[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghE1),
                  data.type = "logratio",
                  sampleid = colnames(cghE1[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
dnacopy.out <- segment(smoothed.CNA.object)

## wavelets
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)



### plateau plots
par(ask = TRUE)
plateauPlot(dnacopy.out)
plateauPlot(wave.out, cghE1[, 5:7])




library(ADaCGH)
data(cghMCRe)
MidPoint <- cghMCRe$Start + 0.5 * (cghMCRe$End - cghMCRe$Start)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
reorder <- order(chrom.numeric,
                 MidPoint,
                 cghMCRe$Start,
                 cghMCRe$End,
                 cghMCRe$Name)
cghMCRe <- cghMCRe[reorder, ]
chrom.numeric <- chrom.numeric[reorder]
MidPoint <- MidPoint[reorder]

hmm.out <- pSegmentHMM(cghMCRe[, 5:7], chrom.numeric)
glad.out <- pSegmentGLAD(cghMCRe[, 5:7], chrom.numeric)
cghseg.out <- pSegmentCGHseg(cghMCRe[, 5:7], chrom.numeric)
ace.out <- pSegmentACE(cghMCRe[, 5:7], chrom.numeric)
wave.out <- pSegmentWavelets(cghMCRe[, 5:7], chrom.numeric)
cbs.out <- pSegmentDNAcopy(cghMCRe[, 5:7], chrom.numeric)
cbs.nm.out <- pSegmentDNAcopy(cghMCRe[, 5:7], chrom.numeric, merge = FALSE)

psw.pos.out <- pSegmentPSW(cghMCRe[, 5:7], chrom.numeric, sign = 1)
psw.neg.out <- pSegmentPSW(cghMCRe[, 5:7], chrom.numeric, sign = -1)

## BioHMM is the only one that uses distances
biohmm.out <- pSegmentBioHMM(cghMCRe[, 5:7], chrom.numeric, MidPoint)














############# Plots and priting.

setwd("/tmp/o3") ## all slaves need a common dir to read and write.
mpiInit()

data(cghE1)
tmpchr <- sub("chr", "", cghE1$Chromosome)
chrom.numeric <- as.numeric(as.character(tmpchr))
chrom.numeric[tmpchr == "X"] <- 23
chrom.numeric[tmpchr == "Y"] <- 24
rm(tmpchr)
reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]


## run all methods

hmm.out <- pSegmentHMM(cghE1[, 5:7], chrom.numeric)
glad.out <- pSegmentGLAD(cghE1[, 5:7], chrom.numeric)
cghseg.out <- pSegmentCGHseg(cghE1[, 5:7], chrom.numeric)
ace.out <- pSegmentACE(cghE1[, 5:7], chrom.numeric)
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)
wave.nm.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric, merge = FALSE)
cbs.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric)
cbs.nm.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric, merge = FALSE)
cbs.nm.ns.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric, merge = FALSE,
                                 smooth = FALSE)
psw.pos.out <- pSegmentPSW(cghE1[, 5:7], chrom.numeric, sign = 1)
psw.neg.out <- pSegmentPSW(cghE1[, 5:7], chrom.numeric, sign = -1)

## BioHMM is the only one that uses distances
## it is the slowest, so do only two
biohmm.out <- pSegmentBioHMM(cghE1[, 5:6], chrom.numeric, cghE1$UG.Start)

segmentPlot(hmm.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")


segmentPlot(glad.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")


segmentPlot(cghseg.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")


segmentPlot(wave.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")


segmentPlot(wave.nm.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(cbs.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

### Te following will fail, because it is now deprecated
## segmentPlot(cbs.nm.out,
##             geneNames = cghE1[, 1],
##             chrom.numeric = chrom.numeric,
##             cghdata = cghE1[, 5:7],
##             idtype = "ug",
##             organism = "Hs")

segmentPlot(psw.pos.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(psw.neg.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")



segmentPlot(biohmm.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:6],
            idtype = "ug",
            organism = "Hs")

## need to choose fdr

ace.out.sum <- summary(ace.out)
segmentPlot(ace.out.sum, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")





#### Writing
common <- cghE1[, -c(5:7)]

writeResults(hmm.out, cghE1[, 5:7], common)
writeResults(glad.out, cghE1[, 5:7], common)        
writeResults(cghseg.out, cghE1[, 5:7], common)   
writeResults(ace.out.sum, cghE1[, 5:7], common)
writeResults(wave.out, cghE1[, 5:7], common)
writeResults(wave.nm.out, cghE1[, 5:7], common)
writeResults(cbs.out, cghE1[, 5:7], common)
writeResults(psw.pos.out, cghE1[, 5:7], common)
writeResults(psw.neg.out, cghE1[, 5:7], common)
writeResults(biohmm.out, cghE1[, 5:6], common)
