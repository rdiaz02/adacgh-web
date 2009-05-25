\name{tempdir2}
\alias{tempdir2}
\title{Create a directory for output}
\description{
  Creates a (unique) temporary directory for storing all output.
  A simple utility function.
}
\usage{
tempdir2()
}
%- maybe also 'usage' for other objects documented here.
% \arguments{
% }
\details{
Just a convenience wrapper.
}
\value{
  This function is used for its side effects. Produce a temporary
  directory.
  It returns the name of the temporary directory, so you can then move
  all further work to that directory.
}

\references{
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
}

\author{ Ramon Diaz-Uriarte \email{rdiaz02@gmail.com}}


\examples{
## See examples under \link{pSegmentDNAcopy}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }% __ONLY ONE__ keyword per line
