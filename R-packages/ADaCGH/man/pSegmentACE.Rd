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
  A list with as many components as samples. For each sample, a list
  with as many sublists as chromosomes. Each chromosomes list has as
  components several vectors: "sample1", "FDR", "calledGenes","Sdev",
  "called", "ACEPgene", "first","last" . This output, however is rarely
  of immediate use. This output is later used by the function \code{\link{summary.ACE}}.
}

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

## Recall: we must reorder the data by chromosome and
## by position within chromosome
reorder <- order(chrom.numeric,
                 cghMCRe$Start,
                 cghMCRe$End,
                 cghMCRe$Name)
cghMCRe <- cghMCRe[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

ace.out <- pSegmentACE(cghMCRe[, 5:7], chrom.numeric)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
