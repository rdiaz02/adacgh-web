##############################################################
###      
###      To verify our implementation of MCR is working
###
##############################################################


### Helper functions

stretchCNAoutput <- function(object) {
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






############ Do segmentation

library(ADaCGH)

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

CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
out.dnacopy <- segment(smoothed.CNA.object, nperm = 10000,
                       undo.splits = "prune")

out.adacgh <- pSegmentDNAcopy(cghMCRe[, 5:7], mergeSegs = FALSE,
                              nperm = 10000,
                              chrom.numeric = chrom.numeric)


### Verify output of DNAcopy and our segmentation are the same
verifCNA <- mapply(function(x, y) all.equal(x[, 2], round(y[, 2], 4)),
                    stretchCNAoutput(out.dnacopy), out.adacgh$segm) 

### Verify we can reconstruct the DNAcopy object properly

reconstructed <- ADaCGH:::constructSegmObj(out.adacgh$segm,
                                           chrom.numeric, cghMCRe[, 5:7],
                                           Pos = 1:nrow(cghMCRe))

rec.verif(out.dnacopy, reconstructed)



#########    MCR

class(out.dnacopy[[1]]) <- "data.frame"

ga <- 500
alo <- 0.3
ahi <- 0.7
recu <- 90

MCR(cghMCR(reconstructed,  gapAllowed = ga,  alteredLow = alo,
           alteredHigh = ahi, recurrence = recu))
MCR(cghMCR(out.dnacopy,  gapAllowed = ga,  alteredLow = alo,
           alteredHigh = ahi, recurrence = recu))

doMCR(out.adacgh$segm, chrom.numeric, cghMCRe[, 5:7],
      ga, alo, ahi, recu, Pos = 1:nrow(cghMCRe))


rec2 <- ADaCGH:::constructSegmObj(x, rep(1, 2000), y,
                                  Pos = Pos)

MCR(cghMCR(rec2,  gapAllowed = ga,  alteredLow = alo,
           alteredHigh = ahi, recurrence = recu))

doMCR(x, rep(1, 2000), y, ga, alo, ahi, recu)



o1 <- MCR(cghMCR(rec2,  gapAllowed = ga,  alteredLow = 0.1,
                 alteredHigh = 0.9, recurrence = 95))




doMCR(x, rep(1, 2000), y, 500, 0.03, 0.97, 75, fsink = "r1.txt", Pos = Pos)
MCR(cghMCR(rec2,  gapAllowed = 500,  alteredLow = 0.03,
                 alteredHigh = 0.97, recurrence = 75))
