\name{pSegmentDNAcopy}
\alias{pSegmentDNAcopy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized version of \code{\link[DNAopy]{segment}} }
\description{
  This function is a parallelized version (suing papply and Rmpi) of \code{\link[DNAopy]{segment}}.
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
\details{
  You should only use this functino if you are using two or more CPUs.
}
\value{
  See \code{\link[DNAcopy]{segment}}
}
\references{
 Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004).
     Circular binary segmentation for the analysis of array-based DNA
     copy number data.  _Biostatistics_ 5: 557-572. <URL:
     http://www.mskcc.org/biostat/~olshena/research.>
     }

\author{ E. S. Venkatraman and Adam Olshen \email{olshena@mskcc.org};
  parallelization by Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }

\seealso{\code{\link[DNAcopy]{segment}} }
\examples{

\dontrun{
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
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }

