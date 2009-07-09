\name{segmentPlot}
\alias{segmentPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Segment plots for aCGH }
\description{
  A (more or less) uniform interface and uniform output to segment plots
for the available  aCGH methods in this package. By the default, it
produces html files (and associated png) with an image map with links to
additional information. You need Python for all this to work.
}
\usage{
segmentPlot(x, geneNames, yminmax, 
            idtype = "ug", organism = "Hs",
            arrays = NULL, chroms = NULL,
            colors = c("orange", "red", "green", "blue", "black"),
            html_js = FALSE,
            superimp = FALSE,
            imgheight = 500,
            genomewide_plot = FALSE,
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The fitted object of the appropriate class. }
  \item{geneNames}{ A vector of gene names. }
    \item{yminmax}{ A vector of the form (min, max) for the minimum and
    maximum used in the figures. }
%  \item{chrom.numeric}{ A (numeric) vector of chromosome indicators. }
%   \item{cghdata}{ The aCGH data. }
%  \item{arraynames}{ A vector of array names. }
  \item{idtype}{ The type of id of the gene name; one of
    ('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl', 'entrez', 'ug')
    corresponding, respectively to None, CNIO ID, Affymetrix, Clone ID,
    Accession Number, Ensembl ID, Entrez ID, UniGene. }
  \item{organism}{ The organism; one of ('None', 'Hs', 'Mm', 'Rn')
    corresponding, respectively, to None, Homo sapiens, Mus musculus
    (mouse) and Rattus norvegicus (rat). }
  \item{arrays}{ The arrays to plot (a vector with the numbers, the
    indices, of the arrays), or NULL for plotting all arrays.}
  \item{chroms}{ The chromosomes to plot (a vector with the numbers, the
    indices, of the chromsomes), or NULL for plotting all chromosomes.}
  
  \item{colors}{Colors used for plotting: a five element vector denoting
    colors for non-altered points, gained, lost, the smoothed line and the
    horizontal line at 0.}
  \item{html_js}{Produce the html figures including lots of Javascript?
    Default is FALSE.}
  \item{superimp}{Produce the figures that superimpose all arrays into a
    single plot? Default is FALSE. Superimposed plots will soon NOT be available!}
  \item{imgheight}{The image height when generating png. Default is
    500.}
  \item{genomewide_plot}{Do a genomewide plot. Default is
    FALSE. Genomewide plots will soon NOT be available.}
  \item{\dots}{ Other arguments, passed to the underlying plotting functions. }
}
\value{
  Used only for its side effects of producing plots. Beware that you
  will get from a few to a lot of HTML and png files written to your
  local file system.
}
\references{
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}. 
 }
\author{ Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }
\note{ 
  Most of the arguments are not really needed in most cases. See
  examples below.

  You NEED a working Python for this function to work and produce the
  imagemaps.
  
}
\seealso{ \code{\link{pSegmentWavelets}},
  \code{\link{pSegmentACE}},
  \code{\link{pSegmentDNAcopy}},
  \code{\link{pSegmentPSW}},
}
\examples{
## See examples under \link{pSegmentDNAcopy}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ IO }% __ONLY ONE__ keyword per line
