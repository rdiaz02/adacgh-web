
## PSW
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


psw.loss.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]),
                        chrom.numeric,
                        sign = -1, nIter = 5000, prec = 100,
                        p.crit = 0.10)
psw.gain.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]),
                        chrom.numeric,
                        sign = 1, nIter = 5000, prec = 100,
                        p.crit = 0.10)









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
library(ADaCGH)
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





#######################################################
#######################################################
#######################################################
###
###       Verify constructSegmOut works
###            
#######################################################
#######################################################
#######################################################


stretchCNAoutput <- function(object) {
## Verify we are doing it OK. Test against original version. We need
## to stretch the output of the original function
    if(!(inherits(object, "DNAcopy")))
        stop("This function can only be applied to DNAcopy objects")
    numarrays <- ncol(object$data) - 2
    stretched <- list()
    for(arraynum in 1:numarrays) {
        obs <- object$data[, 2 + arraynum]
        segmented <-
            object$output[object$output$ID ==
                          colnames(object$data)[2 + arraynum], ]
        smoothed <- object$data$maploc 
        for(i in 1:nrow(segmented)) {
            smoothed[(segmented[i,'loc.end'] >= smoothed) &
                     (segmented[i,'loc.start'] <= smoothed)] <-
                         segmented[i,'seg.mean']
        }
        stretched[[arraynum]] <- cbind(Observed = obs,
                                       Predicted = smoothed)
    }
    return(stretched)
}
rec.verif <- function(x, y) {
    t1 <- all.equal(x$output[, 3], y$output$loc.start)
    t2 <- all.equal(x$output[, 4], y$output$loc.end)
    t3 <- all.equal(x$output[, 6], round(y$output$seg.mean, 4))

    if(t1) print("Start equal")
    if(t2) print("End equal")
    if(t3) print("seg.mean equal") else t3
    if(t1 & t2 & t3) print("***Global result:  OK ***") else stop("Not equal")
}




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
common  <- common[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
out.original <- segment(smoothed.CNA.object, nperm = 10000,
                        undo.splits = "prune")

out1 <- pSegmentDNAcopy(cghMCRe[, 5:7], mergeSegs= FALSE, nperm = 10000,
                        chrom.numeric = chrom.numeric)

## recall the originall segment rounds output
verifCNA <- mapply(function(x, y) all.equal(x[, 2], round(y[, 2], 4)),
                    stretchCNAoutput(out.original), out1$segm) 

## verify reconstruction
reconstructed <- ADaCGH:::constructSegmObj(out1$segm,  chrom.numeric, cghMCRe[, 5:7],
                                  Pos = 1:nrow(cghMCRe))

rec.verif(out.original, reconstructed)



###### another, with "real" pos
CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = common$MidPoint,
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
out.original <- segment(smoothed.CNA.object, nperm = 10000,
                        undo.splits = "prune")
reconstructed <- ADaCGH:::constructSegmObj(out1$segm,  chrom.numeric, cghMCRe[, 5:7],
                                  Pos = common$MidPoint)
rec.verif(out.original, reconstructed)









#######################################################
#######################################################
#######################################################
###
###            Overall function
###            
#######################################################
#######################################################
#######################################################


library(ADaCGH)
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




SegmentPlotWrite(cghE1[, 5:7], chrom.numeric,
                 merge = FALSE, pos = cghE1$UG.Start,
                 idtype = "ug", organism = "Hs",
                 method = "Wavelets",
                 geneNames = cghE1[, 1],
                 commondata = cghE1[, 1:4])



for(mm in c("Wavelets", "DNAcopy", "GLAD", "HMM", "BioHMM", "CGHseg")) {

    cat("\n\n mm is ", mm, "\n\n")
    SegmentPlotWrite(cghE1[, 5:7], chrom.numeric,
                     merge = TRUE, Pos = cghE1$UG.Start,
                     idtype = "ug", organism = "Hs",
                     method = mm,
                     geneNames = cghE1[, 1],
                     commondata = cghE1[, 1:4])
}
