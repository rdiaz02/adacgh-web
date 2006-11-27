\name{pSegmentWavelets}
\alias{pSegmentWavelets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized version of the wavelet-based segmentation analysis
  of Hsu, Grove, et al.}
\description{
  This function wraps and parallelized the code of Hsu and Grove for 
  wavelet-based segmentation of aCGH data. }
\usage{
pSegmentWavelets(acghdata, chrom.numeric, minDiff = 0.25, thrLvl = 3,
                 init.cluster.levels = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acghdata}{The aCGH data  }
  \item{chrom.numeric}{ The (numeric) chromosome indicator }
  \item{minDiff}{ MINimum (absolute) DIFFerence between the medians of
    two adjacent clusters for them to be considered truly different.
    Clusters "closer" together than this are collapsed together to form
    a single cluster. }
  \item{thrLvl}{ The level used for the wavelet thresholding.}
  \item{initClusterLevels}{The initial number of clusters to form. }
}
\details{
  See the original paper for details. See also the
  \code{\link[waveslim]} package.
}
\value{A list with two components:
  \item{Predicted}{A matrix with genes in rows and samples in columns,
    with the "predicted" (smoothed and clustered + merged) values.}
  \item{State}{A matrix with genes in rows and samples in columns,
    with the "State", a categorical state indicator starting at
1. Please note that a state of 2 in one sample might denote a completely
different thing from a state of 2 in another sample. These are the
states for each subject/array, but are not necessarily comparable among
subjects/arrays.}
}
\references{

  Hsu L, Self SG, Grove D, Randolph T, Wang K, Delrow JJ, Loo L, Porter
P. (2005) Denoising array-based comparative genomic hybridization data
using wavelets. \emph{Biostatistics}, \bold{6}:211-26.

  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}. }
\author{ Original code by L. Hsu and D. Grove (see reference for full
  details). Some wrapping and parallelization by Ramon Diaz-Uriarte.}

\note{
  We have only wrapped the code that was kindly provided by Hsu and
Grove, and parallelized a few calls. Their original code is included in
the sources of the package.
}
\seealso{ \code{\link{segmentPlot}} }
\examples{
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
wave.out <- pSegmentWavelets(cghMCRe[, 5:7], chrom.numeric)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
