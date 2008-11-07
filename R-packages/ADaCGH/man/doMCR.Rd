\name{doMCR}
\alias{doMCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find minimal common regions}
\description{
  A wrapper for the functionality implemented in the \code{cghMCR}
  package suitable for objets such as those returned by our
  functions. Output is returned as HTML and text files.
}
\usage{
doMCR(x, chrom.numeric, data, MCR.gapAllowed = 500,
      MCR.alteredLow = 0.03, MCR.alteredHigh = 0.97,
      MCR.recurrence = 75, fsink = "results.txt",
      hsink = "mcr.results.html", Pos, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The "segm" component of an object returned by our pSegment* functions, such as
    \code{\link{pSegmentDNAcopy}}. For instance
    \code{whatevername$semg}, where \code{whatevername} is the result of
    a   previous invocation of \cdoe{pSegment}.}
  \item{chrom.numeric}{See \code{\link{pSegmentDNAcopy}}} 
  \item{data}{The original aCGH data;  genes in rows, subjects or arrays in
    columns. Data should be ordered by chromosome and position within chromosome.}
  \item{MCR.gapAllowed}{ See \code{\link[cghMCR]{MCR}}}
  \item{MCR.alteredLow}{ See \code{\link[cghMCR]{MCR}}}
  \item{MCR.alteredHigh}{ See \code{\link[cghMCR]{MCR}}}
  \item{MCR.recurrence}{ See \code{\link[cghMCR]{MCR}}}
  \item{fsink}{ Name of file where text results are written to.}
  \item{hsink}{ Name of file where HTML results are written to. }
  \item{Pos}{The location (e.g., position in kbases) of each probe in
    the chromosome. A vector of the same length as the number of rows in
  \code{data}}
  \item{\dots}{ Additional arguments (not used now). }
}
\value{
  Used only for its side effects of writing out files.
}





\references{
  Aguirre, Brennan, Bailey, Sinha, Feng, Leo, Zhang, Zhang, Gans,
Bardeesy, Cauwels, Cordon-Cardo, Redston, Depinho,
Chin. (2004). High-resolution characterization of the pancreatic
adenocarcinoma genome. PNAS, 101: 9067--9072.

   Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh2.bioinfo.cnio.es}.

}





\author{Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }
\seealso{ \code{\link[cghMCR]{MCR}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

