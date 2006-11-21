\name{DNAcopyDiagnosticPlots}
\alias{DNAcopyDiagnosticPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic plots for DNAcopy (circular binary segmentation) }
\description{
  Following a suggestion in the original paper (Olshen et al., 2004),
plots data before and after CNA smoothing.
}
\usage{
DNAcopyDiagnosticPlots(CNA.object, CNA.smoothed.object, array.chrom = FALSE, chrom.numeric = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNA.object}{A CNA object obtained with a call to \code{\link[DNAcopy]{CNA}}.}
  \item{CNA.smoothed.object}{A CNA object after applying \code{\link[DNAcopy]{smooth.CNA}}.}
  \item{array.chrom}{If \code{TRUE}, produce plots by chromosome (so you
  will get number of arrays * number of chromosomes plots).}
  \item{chrom.numeric}{ If \code{array.chrom = TRUE}, you need to
    provide the vector with the numeric chromosome indicator.}
}

\value{ Use for it side effects of producing plots.}
\references{
  Olshen, A. B.  and Venkatraman, E. S.  and Lucito, R.  and Wigler,
  M. (2004) Circular binary segmentation for the analysis of array-based
DNA copy number data. \emph{Biostatistics}, \bold{4}, 557--572.

Díaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
web-based application for the analysis of aCGH data. Tech. report.
\url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
\url{http://adacgh.bioinfo.cnio.es}.
}

\author{Ramon Diaz-Uriarte \email{rdiaz02@gmail.com}}
\note{Do not trust these plots blindly; on the other hand, if there are
  too many dots outside the diagonal, and/or you start seeing patterns,
  double check your data and the settings of the smoothing.}
\seealso{\code{\link[DNAcopy]{smooth.CNA}}, \code{\link{WaveletsDiagnosticPlots}}}
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


par(mfrow = c(1, 3))
DNAcopyDiagnosticPlots(CNA.object,
                       smoothed.CNA.object)

par(mfrow = c(3, 3))
par(ask = TRUE)
DNAcopyDiagnosticPlots(CNA.object,
                       smoothed.CNA.object, array.chrom = TRUE,
                       chrom.numeric = chrom.numeric)
} 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
