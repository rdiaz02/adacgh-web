\name{SegmentPlotWrite}
\alias{SegmentPlotWrite}
\title{Segment aCGH data, plot, and write out results}
\description{
  Wrapper to do in a single call the segmentation, plotting, and
  writing out of results from segmentation.
}
\usage{
SegmentPlotWrite(data, chrom, mergeSegs, Pos, idtype,
                 organism, method, geneNames, commondata,
                 colors = c("orange", "red", "green", "blue", "black"),
                 html_js = TRUE,
                 superimp = TRUE,
                 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The aCGH data; genes in rows, subjects or arrays in
    columns. Data should be ordered by chromosome and position within chromosome.}
  \item{chrom}{The (numeric) vector with the chromosome
    indicator}
  \item{mergeSegs}{Merge (genome-wide) the resulting segments using \code{\link[aCGH]{mergeLevels}}?}
  \item{Pos}{The location (e.g., position in kbases) of each probe in
    the chromosome. A vector of the same length as the number of rows in
  \code{data}}


  \item{idtype}{ The type of id of the gene name; one of
    ('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl', 'entrez', 'ug')
    corresponding, respectively to None, CNIO ID, Affymetrix, Clone ID,
    Accession Number, Ensembl ID, Entrez ID, UniGene. }
  \item{organism}{ The organism; one of ('None', 'Hs', 'Mm', 'Rn')
    corresponding, respectively, to None, Homo sapiens, Mus musculus
    (mouse) and Rattus norvegicus (rat). }
  \item{method}{The segmentation method. One of "ACE", "HMM", "BioHMM",
    "GLAD", "DNAcopy", "CGHseg", "PSW", "Wavelets".}
  \item{geneNames}{A vector with the names of the genes or probes, as
    you want to use them to label figures and output.}
  \item{commondata}{A data frame with columns that contain the common
    information about the clones/genes. We do nothing with this, except
    it is part of the output. It is common practice that this data frame
    contains columns ID, Chrosome (need not be
    numeric, can be a factor), Start (the initial position of a clone or
    gene), End (the final position of a clone or gene) and MidPoint
    (often, the mid point or whatever you regard as appropriate
    representative of the "center" of the clone). But you can use
    something else.}
  \item{colors}{Colors used for plotting: a five element vector denoting
    colors for non-altered points, gained, lost, the smoothed line and the
    horizontal line at 0.}
  \item{html_js}{Produce the html figures including lots of Javascript?
    Default is TRUE.}
  \item{superimp}{Produce the figures that superimpose all arrays into a
    single plot? Default is TRUE.}

  \item{\dots}{Additional arguments passed to the underlying pSegment*
    and segmentPlot functions.}
}
\details{
Just a convenience wrapper.
}
\value{
  This function is used for its side effects. Produce plots and write
  out results to a file.  See \code{\link{writeResults}} and \code{\link{segmentPlot}}.
}



\references{
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
}

\author{ Ramon Diaz-Uriarte \email{rdiaz02@gmail.com}}

\seealso{
\code{\link{writeResults}}, \code{\link{segmentPlot}}, \code{\link{pSegment}}
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
reorder <- order(chrom.numeric,
                 cghE1$UG.Start,
                 cghE1$UG.End,
                 cghE1$Name)
cghE1 <- cghE1[reorder, ]
chrom.numeric <- chrom.numeric[reorder]



## one example
SegmentPlotWrite(cghE1[, 5:7], chrom.numeric,
                 merge = FALSE, Pos = cghE1$UG.Start,
                 idtype = "ug", organism = "Hs",
                 method = "Wavelets",
                 geneNames = cghE1[, 1],
                 commondata = cghE1[, 1:4])



## all other methods except PSW and ACE
for(mm in c("DNAcopy", "GLAD", "HMM", "BioHMM", "CGHseg")) {

    cat("\n\n mm is ", mm, "\n\n")
    SegmentPlotWrite(cghE1[, 5:7], chrom.numeric,
                     merge = TRUE, Pos = cghE1$UG.Start,
                     idtype = "ug", organism = "Hs",
                     method = mm,
                     geneNames = cghE1[, 1],
                     commondata = cghE1[, 1:4])
}


}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{ IO }% __ONLY ONE__ keyword per line
