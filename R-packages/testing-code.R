### Some code for testing and checking and catching bugs
rm(list = ls())
library(ADaCGH)
load("p1.RData")  ## run the code in SegmentPlotWirte, and save it here
dir.to.create <- tempdir2()
setwd(dir.to.create)
dir.to.create
mpiInit()


yminmax <- c(min(as.matrix(cghE1[, 5:7])),
             max(as.matrix(cghE1[, 5:7])))


segmentPlot(cbs.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(hmm.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(glad.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(cghseg.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(wave.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(haar.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")


try(segmentPlot(psw.pos.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs"))

try(segmentPlot(psw.neg.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs"))

try(segmentPlot(biohmm.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs"))

segmentPlot(ace.out.sum, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")


SegmentPlotWrite(cghE1[, 5:7], chrom.numeric,
                 merge = FALSE, Pos = cghE1$UG.Start,
                 idtype = "ug", organism = "Hs",
                 method = "Wavelets",
                 geneNames = cghE1[, 1],
                 commondata = cghE1[, 1:4],
                 html_js = TRUE,
                 genomewide_plot = TRUE,
                 superimp = TRUE)



## all other methods except PSW and ACE and BioHMM (this fails because
## of a problem in the library)
for(mm in c("DNAcopy", "GLAD", "HMM", "CGHseg", "HaarSeg", "Wavelets")) {

    cat("\n\n mm is ", mm, "\n\n")
    SegmentPlotWrite(cghE1[, 5:6], chrom.numeric,
                     merge = TRUE, Pos = cghE1$UG.Start,
                     idtype = "ug", organism = "Hs",
                     method = mm,
                     geneNames = cghE1[, 1],
                     commondata = cghE1[, 1:4])

    cat("\n\n with all stuff: mm is ", mm, "\n\n")
    SegmentPlotWrite(cghE1[, 5:6], chrom.numeric,
                     merge = TRUE, Pos = cghE1$UG.Start,
                     idtype = "ug", organism = "Hs",
                     method = mm,
                     geneNames = cghE1[, 1],
                     commondata = cghE1[, 1:4],
                     html_js = TRUE,
                     genomewide_plot = TRUE,
                     superimp = TRUE)


}

## Now try BioHMM

try(SegmentPlotWrite(cghE1[, 5:6], chrom.numeric,
                     merge = TRUE, Pos = cghE1$UG.Start,
                     idtype = "ug", organism = "Hs",
                     method = "BioHMM",
                     geneNames = cghE1[, 1],
                     commondata = cghE1[, 1:4]))
