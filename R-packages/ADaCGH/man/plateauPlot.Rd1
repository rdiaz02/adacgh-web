\name{plateauPlot}
\alias{plateauPlot}
\alias{plateauPlot.DNAcopy}
\alias{plateauPlot.CGH.wave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Plateau plots for aCGH segmentation }
\description{
  Following Olshen et al. and the DNAcopy package, these are plateau
  plots for the DNAcopy and wavelet-based methods of aCGH analsysis.
}
\usage{
plateauPlot(obj, ...)
plateauPlot.CGH.wave(obj, cghdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ An object of class "DNAcopy" (as from
    \code{\link{pSegmentDNAcopy}}), or "CGH.wave" (as from
    \code{\link{pSegmentWavelets}}).} 
  \item{cghdata} {Only for objects of type CGH.wave; the aCGH data.}
  \item{\dots}{ Additional arguments to plotting functions (unused now) }
}
\details{
}
\value{
  A series of plateau plots. First one by sample, next superimposing all
  (but just a superimposition of all the by-sample) and finally one
  where all data are put together (as if they came from a single array)
  and the data plotted again.
}
\references{
  Hsu L, Self SG, Grove D, Randolph T, Wang K, Delrow JJ, Loo L, Porter
P. (2005) Denoising array-based comparative genomic hybridization data
using wavelets. \emph{Biostatistics}, \bold{6}:211-26.

  Olshen, A. B.  and Venkatraman, E. S.  and Lucito, R.  and Wigler,
  M. (2004) Circular binary segmentation for the analysis of array-based
  DNA copy number data. \emph{Biostatistics}, \bold{4}, 557--572.
  \url{http://www.mskcc.org/biostat/~olshena/research}.

  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}. }

 \author{
A lot of code has been taken from the original code in
\code{\link[DNAcopy]{plot.DNAcopy}} by  E. S. Venkatraman and Adam Olshen \email{olshena@mskcc.org};
 modifications and additions for wavelet-based method by by Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }

\seealso{\code{\link[DNAcopy]{segment}},
  \code{\link[DNAcopy]{plot.DNAcopy}},
\code{\link{pSegmentWavelets}}}

\note{  Plateau plots are not provided for merged data, as platau plots are
  used to help in kind of a manual merging.
}

\examples{
data(cghE1)
tmpchr <- sub("chr", "", cghE1$Chromosome)
chrom.numeric <- as.numeric(as.character(tmpchr))
chrom.numeric[tmpchr == "X"] <- 23
chrom.numeric[tmpchr == "Y"] <- 24
rm(tmpchr)
## Recall: we must reorder the data by chromosome and by position within
## chromosome

reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

### DNA copy
CNA.object <- CNA(as.matrix(cghE1[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghE1),
                  data.type = "logratio",
                  sampleid = colnames(cghE1[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
dnacopy.out <- segment(smoothed.CNA.object)

## wavelets
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)


### plateau plots
par(ask = TRUE)
plateauPlot(dnacopy.out)
plateauPlot(wave.out, cghE1[, 5:7])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ misc }% __ONLY ONE__ keyword per line
