\name{summary.ACE}
\alias{summary.ACE}
\alias{summary.ACE.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gains/loss and FDR from an ACE object}
\description{
  Shows the number of genes detected as gains/loss at different levels
  of FDR, allows setting the desired FDR, and returns the assigned
  status of each genes for the selected FDR.
}
\usage{
\method{summary}{ACE}(object, fdr = NULL, html = TRUE,
            outhtml = "ace.fdrtable.html", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class ACE or ACE.array, as returned from  \code{\link{pSegmentACE}} }
  \item{fdr}{ The desired FDR. The actual value used is the closest one
    the desired. If NULL, the FDR selected is the closest one to 0.15. }
  \item{html}{ If TRUE, produce also an HTML table }
  \item{outhtml}{ The name of the file with the HTML output }
  \item{...}{Not used.}
}
\value{
  Calling this function always returns (prints) a table with the number
  of Genes called gains/losses at each of the available FDRs and the
  index (row) corresponding to the selected FDR (or the closest FDR to
  0.15). Additionally, a list with as many components as samples. Each
  list component is a data frame with genes/clones as rows and three
  columns: the chromosome identifier, the original log-ratio, and the
  status (0 for no change, -1 for loss, +1 for gain).
}
\references{
  Lingjaerde OC, Baumbusch LO, Liestol K, Glad I,
Borresen-Dale AL. (2005). CGH-Explorer: a program for analysis of
CGH-data. \emph{Bioinformatics}, \bold{21}: 821--822.

Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
 }
\author{ Oscar M. Rueda \email{omrueda@cnio.es} }

\seealso{ \code{\link{pSegmentACE}}, \code{\link{segmentPlot}} }
\examples{
\dontrun{
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
summary(ace.out, fdr = 0.01)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
