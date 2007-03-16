\name{writeResults}
\alias{writeResults}
\alias{writeResults.CGH.PSW}
\alias{writeResults.CGH.ACE.summary}
\alias{writeResults.CGH.wave}
\alias{writeResults.DNAcopy}
\alias{writeResults.CGHseg}
\alias{writeResults.mergedHMM}
\alias{writeResults.adacghGLAD}
\alias{writeResults.mergedBioHMM}
\alias{writeResults.summaryACE}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write results to an external file}
\description{
  Write segmentation results to an external file, including original
  information (name, chromosome, position, log-ratio) and assigned state
  and/or smoothed mean (depending on what is available from the
  method). (In addition, a set of files with links to PaLS
  \url{http://pals.bioinfo.cnio.es} will be generated when there is
  state information.)
}
\usage{
writeResults(obj, ...)
writeResults.CGH.PSW(obj, commondata, file = "PSW.output.txt", ...)
writeResults.summaryACE(obj, acghdata, commondata, file = NULL, ...)
writeResults.CGH.wave(obj, acghdata, commondata,
                     file = "wavelet.output.txt", ...)
writeResults.DNAcopy(obj, acghdata, commondata, 
                     file = "CBS.output.txt", ...)
writeResults.CGHseg(obj, acghdata, commondata, 
                     file = "CGHseg.output.txt", ...)
writeResults.mergedHMM(obj, acghdata, commondata, 
                     file = "HMM.output.txt", ...)
writeResults.adacghGLAD(obj, acghdata, commondata,
                     file = "GLAD.output.txt", ...)
writeResults.mergedBioHMM(obj, acghdata, commondata, 
                     file = "BioHMM.output.txt", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ The output object from previous analysis; in all cases,
    the output from pSegement except for ACE, where it is the output
    from summary.}
  \item{file}{The name of the output file.}
  \item{acghdata}{The aCGH data set; see \code{\link{pSegmentWavelets}}.}
  \item{commondata}{The common data; see \code{common.data} in 
    \code{\link{pSegmentPSW}}.}
  \item{...}{Other arguments passed to the underlying functions.}
}

\value{
  A file is written in your system with columns corresponding to the
  original data (including ID, Chromosome, etc) and the results of the segmentations.
}
\references{
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
}
\author{ Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }
\seealso{
\code{\link{pSegmentACE}}, \code{\link{pSegmentWavelets}},
\code{\link{pSegmentPSW}},  \code{\link{pSegmentDNAcopy}},
\code{\link{summary.ACE}}, 
 }
 \examples{
setwd("/tmp/o3") ## all slaves need a common dir to read and write.
\dontrun{mpiInit()}

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

## need to choose fdr
ace.out.sum <- summary(ace.out)

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



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ IO }% __ONLY ONE__ keyword per line
\keyword{file}
