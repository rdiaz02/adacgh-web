
#####################


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
summary(ace.out)


##################

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

segmentPlot(dnacopy.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")



segmentPlot(ace.out.sum, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)


## PSW
psw2.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]), chrom.numeric,
                       sign = - 1, nIter = 5000, prec = 100, p.crit = 0.10)

## wavelets
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)


acceptedIDTypes = ('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl', 'entrez', 'ug')
acceptedOrganisms = ('None', 'Hs', 'Mm', 'Rn')

