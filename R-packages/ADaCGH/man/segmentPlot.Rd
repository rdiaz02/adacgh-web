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
segmentPlot(x, geneNames, chrom.numeric = NULL, cghdata = NULL,
            arraynames = NULL, idtype = "ug", organism = "Hs",
            yminmax = NULL, numarrays = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The fitted object of the appropriate class. }
  \item{geneNames}{ A vector of gene names. }
  \item{chrom.numeric}{ A (numeric) vector of chromosome indicators. }
  \item{cghdata}{ The aCGH data. }
  \item{arraynames}{ A vector of array names. }
  \item{idtype}{ The type of id of the gene name; one of
    ('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl', 'entrez', 'ug')
    corresponding, respectively to None, CNIO ID, Affymetrix, Clone ID,
    Accession Number, Ensembl ID, Entrez ID, UniGene. }
  \item{organism}{ The organism; one of ('None', 'Hs', 'Mm', 'Rn')
    corresponding, respectively, to None, Homo sapiens, Mus musculus
    (mouse) and Rattus norvegicus (rat). }
  \item{yminmax}{ A vector of the form (min, max) for the minimum and
    maximum used in the figures. }
  \item{numarrays}{ The number of arrays to plot. }
  \item{\dots}{ Other arguments, passed to the underlying plotting functions. }
}
\details{
}
\value{
  Used only for its side effects of producing plots. Beware that you
  will get from a few to a lot of HTML and png files writtend to your
  local file system.
}
\references{
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}. }
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

\dontrun{

## This will not work unless you have a working Python installation

setwd("/tmp/o3") ## all slaves need a common dir to read and write.
mpiInit()

data(cghE1)
tmpchr <- sub("chr", "", cghE1$Chromosome)
chrom.numeric <- as.numeric(as.character(tmpchr))
chrom.numeric[tmpchr == "X"] <- 23
chrom.numeric[tmpchr == "Y"] <- 24
rm(tmpchr)
### we need the data ordered
reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]


## run all methods
hmm.out <- pSegmentHMM(cghE1[, 5:7], chrom.numeric)
glad.out <- pSegmentGLAD(cghE1[, 5:7], chrom.numeric)
cghseg.out <- pSegmentCGHseg(cghE1[, 5:7], chrom.numeric)
ace.out <- pSegmentACE(cghE1[, 5:7], chrom.numeric)
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)
wave.nm.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric, merge = FALSE)
cbs.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric)
cbs.nm.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric, merge = FALSE)
cbs.nm.ns.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric, merge = FALSE,
                                 smooth = FALSE)
psw.pos.out <- pSegmentPSW(cghE1[, 5:7], chrom.numeric, sign = 1)
psw.neg.out <- pSegmentPSW(cghE1[, 5:7], chrom.numeric, sign = -1)

## BioHMM is the only one that uses distances
## it is the slowest, so do only two
biohmm.out <- pSegmentBioHMM(cghE1[, 5:6], chrom.numeric, cghE1$UG.Start)


### plot all
segmentPlot(hmm.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(glad.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(cghseg.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(wave.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(wave.nm.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(cbs.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(psw.pos.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(psw.neg.out,
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(biohmm.out, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:6],
            idtype = "ug",
            organism = "Hs")

## need to choose fdr
ace.out.sum <- summary(ace.out)
segmentPlot(ace.out.sum, 
            geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ IO }% __ONLY ONE__ keyword per line
