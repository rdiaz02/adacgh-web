\name{pSegmentDNAcopy}
\alias{pSegmentDNAcopy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized version of DNAcopy's segment }
\description{
  This function is a parallelized version (using papply and Rmpi) of \code{\link[DNAopy]{segment}}.
}
\usage{

pSegmentDNAcopy(x, alpha = 0.01, nperm = 10000,
                p.method = c("hybrid","perm"),
                kmax = 25, nmin = 200, eta = 0.05, window.size = NULL,
                overlap = 0.25, trim = 0.025,
                undo.splits = c("none","prune", "sdundo"),
                undo.prune = 0.05, undo.SD = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{See \code{\link[DNAcopy]{segment}}}
  \item{alpha}{See \code{\link[DNAcopy]{segment}} }
  \item{nperm}{See \code{\link[DNAcopy]{segment}} }
  \item{p.method}{See \code{\link[DNAcopy]{segment}} }
  \item{kmax}{See \code{\link[DNAcopy]{segment}} }
  \item{nmin}{See \code{\link[DNAcopy]{segment}} }
  \item{eta}{See \code{\link[DNAcopy]{segment}} }
  \item{window.size}{See \code{\link[DNAcopy]{segment}} }
  \item{overlap}{See \code{\link[DNAcopy]{segment}} }
  \item{trim}{See \code{\link[DNAcopy]{segment}} }
  \item{undo.splits}{See \code{\link[DNAcopy]{segment}} }
  \item{undo.prune}{See \code{\link[DNAcopy]{segment}} }
  \item{undo.SD}{See \code{\link[DNAcopy]{segment}} }
}
\value{
  See \code{\link[DNAcopy]{segment}}
}
\references{
  Hsu L, Self SG, Grove D, Randolph T, Wang K, Delrow JJ, Loo L, Porter
P. (2005) Denoising array-based comparative genomic hybridization data
using wavelets. \emph{Biostatistics}, \bold{6}:211-26.

  Olshen, A. B.  and Venkatraman, E. S.  and Lucito, R.  and Wigler,
  M. (2004) Circular binary segmentation for the analysis of array-based
  DNA copy number data. \emph{Biostatistics}, \bold{4}, 557--572.
  \url{http://www.mskcc.org/biostat/~olshena/research}.

  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.

}




\author{
A lot of code has been taken from the original code in
\code{\link[DNAcopy]{plot.DNAcopy}} by  E. S. Venkatraman and Adam Olshen \email{olshena@mskcc.org};
modifications and additions for wavelet-based method by by Ramon
Diaz-Uriarte \email{rdiaz02@gmail.com}
}


\seealso{\code{\link[DNAcopy]{segment}}
  \code{\link[DNAcopy]{plot.DNAcopy}}
  \code{\link{pSegmentWavelets}}}



\examples{
\dontrun{mpiInit()}

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

out1 <- pSegmentDNAcopy(smoothed.CNA.object, merge = FALSE, nperm = 10000)

## Yes, tilingArray has an identically named function, but DNAcopy does
## not have a namespace. Problems would be caught by different params.

out.original <- segment(smoothed.CNA.object, nperm = 10000,
                        undo.splits = "prune")

## Verify we are doing it OK. Test against original version. We need
## to stretch the output of the original function

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


## recall the originall segment rounds output
verifCNA <- mapply(function(x, y) all.equal(x[, 2], round(y[, 2], 4)),
                    stretchCNAoutput(out.original), out1) 

if(! all(verifCNA)) stop("Differences in original and parallelized DNAcopy.\n Rerun test")


}
\keyword{ nonparametric }

