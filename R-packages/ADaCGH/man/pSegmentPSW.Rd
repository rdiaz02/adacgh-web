\name{pSegmentPSW}
\alias{pSegmentPSW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{   Parallelized version of \code{\link[cgh]{sw}} }
\description{
  Parallelized version of \code{\link[cgh]{sw}}, with some changes in
  the interfaz to make it more similar to other functions in this
  package and including the Smith-Waterman algorithm run, the
  thresholding, the permutation, and the
  robustness calculation.
}
\usage{
pSegmentPSW(common.dat, acghdata, chrom.numeric, sign = -1,
            nIter = 1000, prec = 100,
            p.crit = 0.10, name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{common.dat}{A data frame with columns that contain the common
    information about the clones/genes. We do nothing with this, except
    it is part of the output. It is common practice that this data frame
    contains columns ID, Chrosome (need not be
    numeric, can be a factor), Start (the initial position of a clone or
    gene), End (the final position of a clone or gene) and MidPoint
    (often, the mid point or whatever you regard as appropriate
    representative of the "center" of the clone). But you can use
    something else.}
  \item{acghdata}{ The aCGH data }
  \item{chrom.numeric}{ A numeric vector with the chromosome indicator }
  \item{sign}{ See \code{\link[cgh]{sw.threshold}} }
  \item{nIter}{ See \code{\link[cgh]{sw.perm.test}} }
  \item{prec}{ See \code{\link[cgh]{sw.rob}} }
  \item{p.crit}{The largest p-value for which we want a region to be shown, in red,in the plot.}
  \item{name}{ Sample name or ID }
}
\details{
  You should only use this function if you are using two or more CPUs.

  You are expected to provide both a Chromosome in the common.dat and a
  numeric chromosome indicator. Given the first, we could obtain the
  second, but how we would do that might not be how you would do it, so
  we prefer you decide how you want the mapping from your "Chromosome"
  to a numeric representation.
}
\value{A list with two lists as components.
  \item{Data}{A data frame. The first columns are those from the data
    frame you used as common.dat. The next columns are trios of columns
    where each three successive columns correspond to the "Sign", "Robust", and
    "p.value" components of each sample. }
  \item{plotData}{A list with as many sublists as samples. This is the
    information used for the plots. Each sublist has components
\code{logratio},\code{sign},\code{rob},\code{swt.run},\code{swt.perm},\code{p.crit},\code{chrom},\code{p.crit.bonferroni}.
The meaning of all except \code{swt.run},\code{swt.perm} is evident;
\code{swt.perm} is obtained from a call to
\code{\link[cgh]{sw.perm.test}}. \code{swt.run} is itself a list, with
two sublists, \code{length} and \code{start}, that denote the starting
position and lenght of each of the segments considered by the algorithm.
  }

} \references{

  Price TS, Regan R, Mott R, Hedman A, Honey B, Daniels RJ,
Smith L, Greenfield A, Tiganescu A, Buckle V, Ventress N, Ayyub H,
Salhan A, Pedraza-Diaz S, Broxholme J, Ragoussis J, Higgs DR, Flint J,
Knight SJ. (2005) SW-ARRAY: a dynamic programming solution for the
identification of copy-number changes in genomic DNA using array
comparative genome hybridization data. \emph{Nucleic Acids Res.},
\bold{33}:3455-64. 


  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.
 }
 \author{
   Original code by T.S.Price;  parallelization and wrapping by Ramon Diaz-Uriarte \email{rdiaz02@gmail.com}
 }
 
 \seealso{
   \code{\link[cgh]{sw.threshold}},
   \code{\link[cgh]{sw.perm.test}},
   \code{\link[cgh]{sw.rob}},
   \code{\link{segmentPlot}}
 }
 
 \examples{
data(cghMCRe)
common <- cghMCRe[, -c(5:7)]
common$MidPoint <- common$Start + 0.5 * (common$End - common$Start)
colnames(common)[1] <- "ID"
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

## We use nIter = 100 for the example for speed reasons;
## you probably want 1000 or more.

psw.out <- pSegmentPSW(common, as.matrix(cghMCRe[, 5:7]), chrom.numeric,
                       sign = - 1, nIter = 100, prec = 100, p.crit = 0.10)

## so you don't want all the stuff in common? just use a smaller common.dat
common2 <- cghMCRe[, -c(5:7)]
psw2.out <- pSegmentPSW(common2, as.matrix(cghMCRe[, 5:7]), chrom.numeric,
                       sign = - 1, nIter = 100, prec = 100, p.crit = 0.10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

