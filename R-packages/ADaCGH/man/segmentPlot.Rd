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
segmentPlot(x, geneNames, chrom.numeric = NULL, cghdata = NULL, arraynames = NULL, idtype = "ug", organism = "Hs", superimposed = FALSE, html = TRUE, yminmax = NULL, numarrays = NULL, ...)
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
  \item{superimposed}{ If FALSE (the dafault) you get a figure (and thus
  HTML file) for each sample by chromosome. If TRUE, all the samples are
  collapse, so you get as many figures as chromosomes. Might be useful
  to see overall patterns.}
  \item{html}{ If TRUE, produce HTML output.}
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

## ACE
ace.out <- pSegmentACE(cghE1[, 5:7], chrom.numeric)
ace.out.sum <- summary(ace.out)
segmentPlot(ace.out.sum, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(ace.out.sum, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

## DNA copy + merging
CNA.object <- CNA(as.matrix(cghE1[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghE1),
                  data.type = "logratio",
                  sampleid = colnames(cghE1[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)
dnacopy.out <- segment(smoothed.CNA.object)
merged.out <- mergeDNAcopy(dnacopy.out)

segmentPlot(merged.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(merged.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

segmentPlot(dnacopy.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")

segmentPlot(dnacopy.out, geneNames = cghE1[, 1],
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

## PSW
psw.loss.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]),
                        chrom.numeric,
                        sign = -1, nIter = 5000, prec = 100,
                        p.crit = 0.10)
psw.gain.out <- pSegmentPSW(cghE1[, -c(5:7)], as.matrix(cghE1[, 5:7]),
                        chrom.numeric,
                        sign = 1, nIter = 5000, prec = 100,
                        p.crit = 0.10)

segmentPlot(psw.loss.out, geneNames = cghE1[, 1],
            idtype = "ug",
            organism = "Hs")
segmentPlot(psw.gain.out, geneNames = cghE1[, 1],
            idtype = "ug",
            organism = "Hs")


## wavelets
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)
segmentPlot(wave.out, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs")
segmentPlot(wave.out, geneNames = cghE1[, 1],
            chrom.numeric = chrom.numeric,
            cghdata = cghE1[, 5:7],
            idtype = "ug",
            organism = "Hs",
            superimposed = TRUE)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ IO }% __ONLY ONE__ keyword per line
