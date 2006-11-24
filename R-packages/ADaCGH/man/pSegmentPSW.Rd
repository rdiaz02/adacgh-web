\name{pSegmentPSW}
\alias{pSegmentPSW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{   Parallelized version of \code{\link[cgh]{sw}} }
\description{
  Parallelized version of \code{\link[cgh]{sw}}, with some changes in
  the interfaz to make it more similar to other functions in this
  package and including the Smith-Waterman algorithm run, the
  thresholding, the permutation, and the
  robustness calculation.
}
\usage{
pSegmentPSW(common.dat, acghdata, chrom.numeric, sign, nIter, prec,
            p.crit = 0.10, name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{common.dat}{A data frame with columns ID, Chrosome (need not be
    numeric, can be a factor), Start (the initial position of a clone or
    gene), End (the final position of a clone or gene) and MidPoint
    (often, the mid point or whatever you regard as appropriate
    representative of the "center" of the clone).}
  \item{acghdata}{ The aCGH data }
  \item{chrom.numeric}{ A numeric vector with the chromosome indicator }
  \item{sign}{ See \code{\link[cgh]{sw.threshold}} }
  \item{nIter}{ See \code{\link[cgh]{sw.perm.test}} }
  \item{prec}{ See \code{\link[cgh]{sw.rob}} }
  \item{p.crit}{The largest p-value for which we want a region to be shown, in red,in the plot.}
  \item{name}{ Sample name or ID }
}
\details{
  You should only use this function if you are using two or more CPUs.

  You are expected to provide both a Chromosome in the common.dat and a
  numeric chromosome indicator. Given the first, we could obtain the
  second, but how we would do that might not be how you would do it, so
  we prefer you decide how you want the mapping from your "Chromosome"
  to a numeric representation.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (common.dat, acghdata, chrom.numeric, sign, nIter, prec, 
    name = NULL) 
{
    numarrays <- ncol(acghdata)
    ncrom <- length(unique(chrom.numeric))
    out$Data <- common.dat
    out$plotData <- list()
    for (i in 1:numarrays) {
        tmp <- my.sw3(acghdata[, i], chrom.numeric, sign = sign, 
            p.crit = p.crit, main = colnames(acghdata)[i], nIter = nIter, 
            prec = prec, name = paste(name, colnames(acghdata)[i], 
                sep = ""))
        out$Data <- cbind(out$Data, tmp$out)
        p.crit.bonferroni <- tmp$plotdat$p.crit/ncrom
        out$plotData[[i]] <- c(tmp$plotdat, p.crit.bonferroni = p.crit.bonferroni)
    }
    class(out) <- c(class(out), "CGH.PSW")
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
