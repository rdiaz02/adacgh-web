\name{pSegmentACE}
\alias{pSegmentACE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized version of Analysis of Copy Errors method }
\description{
  An R interface (with parallelization) to a C implementation of the ACE method for aCGH
  arrays. 
}
\usage{
pSegmentACE(acghdata, chrom.numeric, echo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acghdata}{ The aCGH data }
  \item{chrom.numeric}{The (numeric) chromosome indicator }
  \item{echo}{ More verbose intermediate output }
}
\details{
  The C code is based on the original Java code by Lingjaerde, as
  provided in his Java implementation of ACE.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...  }

\references{

  Lingjærde OC, Baumbusch LO, Liestøl K, Glad I,
Børresen-Dale AL. (2005). CGH-Explorer: a program for analysis of
CGH-data. \emph{Bioinformatics}, \bold{21}: 821--822.

Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.

}

\author{ Oscar M. Rueda \email{omrueda@cnio.es} }

\seealso{ \code{\link{segmentPlot}}, \code{\link{summary.ACE}} }
\examples{
data(cghMCRe)
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
ace.out <- pSegmentACE(cghMCRe[, 5:7], chrom.numeric)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
