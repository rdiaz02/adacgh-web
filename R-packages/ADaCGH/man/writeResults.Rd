\name{writeResults}
\alias{writeResults, writeResults.CGH.PSW, writeResults.CGH.ACE.summary,
  writeResults.CGH.wave,  writeResults.DNAcopy,}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write results to an external file}
\description{
  Write segmentation results to an external file, including original
  information (name, chromosome, position, log-ratio) and assigned state
  or smoothed mean.
}
\usage{
writeResults(obj, ...)
writeResults.CGH.PSW(obj, file = "PSW.output.txt")
writeResults.CGH.ACE.summary(obj, file = "ACE.output.txt")
writeResults.CGH.wave(obj, acghdata, commondata,
                                  file = "wavelet.output.txt")
writeResults.DNAcopy(obj, acghdata, commondata, merged = NULL,
                                 output = "CBS.output.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ The output object from previous analysis; in all cases,
    the output from pSegement except for ACE, where it is the output
    from summary.}
  \item{\dots}{ Optionally, the file name. For the wavelet-based method,
    you need to provide the commondata (see commondata in, e.g.,
    \code{\link{pSegmentPSW}}) and the acghdata; for DNAcopy you also
    need, in addition to the commondata and acghdata, the output from
    \code{\link{mergeDNAcopy}} if you want the output to reflect the merging. }
}

\value{
  A file is written in your system with columns corresponding to the
  original data (including ID, Chromosome, etc) and the results of the segmentations.
}
\references{
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
}
\author{ Ramon diaz-Uriarte \email{rdiaz02@gmail.com} }
\seealso{
\code{\link{pSegmentACE}}, \code{\link{pSegmentWavelets}},
\code{\link{pSegmentPSW}},  \code{\link{pSegmentDNAcopy}},
\code{\link{mergeDNAcopy}}, \code{\link{summary.ACE}}, 
 }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ IO }% __ONLY ONE__ keyword per line
\keyword{file}
