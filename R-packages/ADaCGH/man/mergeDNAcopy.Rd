\name{mergeDNAcopy}
\alias{mergeDNAcopy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Merge smoothed levels from DNAcopy }
\description{
  Merges smoothed levels from DNAcopy's segementation, and assigns them
  to status gain/loss/no-change (as in Willenbrock and Fridlyand, 2005).
}
\usage{
mergeDNAcopy(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class"DNAcopy", as obtained from
    \code{\link[DNAcopy]{segment}} or \code{\link{pSegmentDNAcopy}}.}
}
\details{
  We use the functionality provided by \code{\link[aCGH]{mergeLevels}}
  and remap the merged segments to just three classes: no-change (the
  merged segment closest to 0), gained (those above the "no-change" one)
  and lost (those below "no-change").
}
\value{
A list with two components.
\item{chrom.numeric}{A vector with the chromosome identifier}
\item{segm}{A list with components with as many components as
  samples. Each component is a matrix of dimensions number of genes *
  3. The three columns are \code{merged.mean}, the predicted mean of the
  (merged) segments, \code{obs}, the original value of the log2 ratio,
  and \code{alteration}, the mapping to gained (1), lost (-1), or no
  change (0)}
}
\references{
  Olshen, A. B.  and Venkatraman, E. S.  and Lucito, R.  and Wigler,
  M. (2004) Circular binary segmentation for the analysis of array-based
  DNA copy number data. \emph{Biostatistics}, \bold{4}, 557--572.
  \url{http://www.mskcc.org/biostat/~olshena/research}.
  
  Willenbrock, H. and Fridlyand, J. (2005). A comparison study: applying
  segmentation to array CGH data for downstream
  analyses. \emph{Bioinformatics}, \bold{21}, 4084--4091.
  
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
  
}
\author{ Oscar M. Rueda \email{omrueda@cnio.es} }
\note{ We use the default parameters for
  \code{\link[aCGH]{mergeLevels}}; it would be simple to allow to pass
  them as parameters of this function.}

\seealso{\code{\link[DNAcopy]{segment}},  \code{\link{pSegmentDNAcopy}}, \code{\link[aCGH]{mergeLevels}}
}
\examples{

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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

