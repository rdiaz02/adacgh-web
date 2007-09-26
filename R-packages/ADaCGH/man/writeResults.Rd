\name{writeResults}
\alias{writeResults}
\alias{writeResults.CGH.PSW}
\alias{writeResults.CGH.ACE.summary}
\alias{writeResults.CGH.wave}
\alias{writeResults.DNAcopy}
\alias{writeResults.CGHseg}
\alias{writeResults.mergedHMM}
\alias{writeResults.adacghGLAD}
\alias{writeResults.mergedBioHMM}
\alias{writeResults.summaryACE}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write results to an external file}
\description{
  Write segmentation results to an external file, including original
  information (name, chromosome, position, log-ratio) and assigned state
  and/or smoothed mean (depending on what is available from the
  method). (In addition, a set of files with links to PaLS
  \url{http://pals.bioinfo.cnio.es} will be generated when there is
  state information.)
}
\usage{
writeResults(obj, ...)
\method{writeResults}{CGH.PSW}(obj, commondata, file = "PSW.output.txt", ...)
\method{writeResults}{summaryACE}(obj, acghdata, commondata, file = NULL, ...)
\method{writeResults}{CGH.wave}(obj, acghdata, commondata,
                     file = "wavelet.output.txt", ...)
\method{writeResults}{DNAcopy}(obj, acghdata, commondata, 
                     file = "CBS.output.txt", ...)
\method{writeResults}{CGHseg}(obj, acghdata, commondata, 
                     file = "CGHseg.output.txt", ...)
\method{writeResults}{mergedHMM}(obj, acghdata, commondata, 
                     file = "HMM.output.txt", ...)
\method{writeResults}{adacghGLAD}(obj, acghdata, commondata,
                     file = "GLAD.output.txt", ...)
\method{writeResults}{mergedBioHMM}(obj, acghdata, commondata, 
                     file = "BioHMM.output.txt", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ The output object from previous analysis; in all cases,
    the output from pSegement except for ACE, where it is the output
    from summary.}
  \item{file}{The name of the output file.}
  \item{acghdata}{The aCGH data set; see \code{\link{pSegmentWavelets}}.}
  \item{commondata}{The common data; see \code{common.data} in 
    \code{\link{pSegmentPSW}}.}
  \item{...}{Other arguments passed to the underlying functions.}
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
\author{ Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }
\seealso{
\code{\link{pSegmentACE}}, \code{\link{pSegmentWavelets}},
\code{\link{pSegmentPSW}},  \code{\link{pSegmentDNAcopy}},
\code{\link{summary.ACE}}, 
 }
 \examples{
## See examples under \code{\link{pSegmentDNAcopy}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ IO }% __ONLY ONE__ keyword per line
\keyword{file}
