\name{pSegmentDNA}
\alias{pSegmentDNAcopy}
\alias{pSegmentACE}
\alias{pSegmentHMM}
\alias{pSegmentGLAD}
\alias{pSegmentBioHMM}
\alias{pSegmentCGHseg}
\alias{pSegmentPSW}
\alias{pSegmentWavelets}
\alias{pSegment}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized versions of several aCGH segementation algorithms/methods}

\description{
  These functions parallelize several segmentation algorithms. is a parallelized version (using papply and Rmpi) of
\code{\link[DNAopy]{segment}}.  }


\usage{

pSegmentDNAcopy(x, chrom.numeric, mergeSegs = TRUE, smooth = TRUE,
                alpha = 0.01, nperm = 10000,
                kmax = 25, nmin = 200, eta = 0.05, 
                overlap = 0.25, trim = 0.025,
                undo.prune = 0.05, undo.SD = 3,
                merge.pv.thresh = 1e-04, merge.ansari.sign = 0.05,
                merge.thresMin = 0.05, merge.thresMax = 0.5, ...)

pSegmentACE(x, chrom.numeric, ...)

pSegmentHMM(x, chrom.numeric, ...)

pSegmentGLAD(x, chrom.numeric, ...)

pSegmentBioHMM(x, chrom.numeric, Pos, ...)

pSegmentCGHseg(x, chrom.numeric, ...)

pSegmentPSW(x, chrom.numeric, common.data,
            sign = -1,
            nIter = 1000, prec = 100,  p.crit = 0.10,
            name = NULL, ...)

pSegmentWavelets(x, chrom.numeric, mergeSegs = TRUE,
                minDiff = 0.25,
                minMergeDiff = 0.05,
                thrLvl = 3, initClusterLevels = 10, ...)



}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The aCGH data; genes in rows, subjects or arrays in
    columns. Data should be ordered by chromosome and position within chromosome.}
  \item{chrom.numeric}{The (numeric) vector with the chromosome
    indicator}
  \item{Pos}{The location (e.g., position in kbases) of each probe in
    the chromosome. A vector of the same length as the number of rows in
  \code{x}}
  \item{mergeSegs}{Merge (genome-wide) the resulting segments using \code{\link[aCGH]{mergeLevels}}?}
  \item{smooth}{Smooth data before segmentation (using \code{\link[DNAcopy]{smoothCNA}})?}
  \item{alpha}{See \code{\link[DNAcopy]{segment}} }
  \item{nperm}{See \code{\link[DNAcopy]{segment}} }
  \item{kmax}{See \code{\link[DNAcopy]{segment}} }
  \item{nmin}{See \code{\link[DNAcopy]{segment}} }
  \item{eta}{See \code{\link[DNAcopy]{segment}} }
  \item{overlap}{See \code{\link[DNAcopy]{segment}} }
  \item{trim}{See \code{\link[DNAcopy]{segment}} }
  \item{undo.prune}{See \code{\link[DNAcopy]{segment}} }
  \item{undo.SD}{See \code{\link[DNAcopy]{segment}} }
  \item{merge.pv.thresh}{See \code{\link[aCGH]{mergeLevels}}}
  \item{merge.ansari.sign}{See \code{\link[aCGH]{mergeLevels}}}
  \item{merge.thresMin}{See \code{\link[aCGH]{mergeLevels}}}
  \item{merge.thresMax}{See \code{\link[aCGH]{mergeLevels}}}
  \item{common.data}{A data frame with columns that contain the common
    information about the clones/genes. We do nothing with this, except
    it is part of the output. It is common practice that this data frame
    contains columns ID, Chrosome (need not be
    numeric, can be a factor), Start (the initial position of a clone or
    gene), End (the final position of a clone or gene) and MidPoint
    (often, the mid point or whatever you regard as appropriate
    representative of the "center" of the clone). But you can use
    something else.}
  \item{sign}{See \code{\link[cgh]{sw.threshold}}}
  \item{nIter}{See \code{\link[cgh]{sw.perm.test}}}
  \item{prec}{See \code{\link[cgh]{sw.rob}}}
  \item{p.crit}{The largest p-value for which we want a region to be
    shown, in red, in the plot.}
  \item{name}{Sample name or ID}
  \item{minDiff}{ Minimum (absolute) difference between the medians of
    two adjacent clusters for them to be considered truly different.
    Clusters "closer" together than this are collapsed together to form
    a single cluster. }
\item{minMergeDiff}{Used only when doing merging in the wavelet method.
  The finall call as to which segments go together is done by a
  \code{mergeLevels}  approach, but an initial collapsing of very close
  values is performed (otherwise, we could end up passing to mergeLevels
  as many initial levels as there are points). }
  \item{thrLvl}{ The level used for the wavelet thresholding.}
  \item{initClusterLevels}{The initial number of clusters to form. }
  \item{...}{Additional arguments; not used.}
}



\details{

In most cases, these are wrappers to the original code, with
modifications for parallelization. For BioHMM and HMM the results are
merged, by the mergeLevels algorithm, as recommended in
\cite{Willenbrock and Fridlyand, 2005}. Merging is also done in GLAD
(with GLAD's own merging algorithm). Merging is the default for DNAcopy
(but this can be set by the user). Merging can also be used with the
wavelet-based method; please note that the later is an experimental
feature implemented by us, and there is no evidence of its performance. 



For PSW, the method implemented here is a small departure from the one
described originally in Price et al.: here we use all the data for a
subject/array for the thresholding, but subsequent analyses (finding
islands and robustness analysis and permutation tests) are carried out
chromosome by chromosome. This should allow to detect cases where a
complete chromosome is gained or lost, as well as cases where smaller
regions are gained or lost. Note that this approach is probably not
advisable for the X chromosome, at least for males. You probably want to
analyze the X separately, since you only expect one copy of the genomic
DNA in males. Finally, the function for calculating the threshold (the
default used by T. Price originally, which we also use ---median + 0.2 *
MAD---) was validated by Price and collaborators on a 75-probe dataset
that spanned the telomeric 2MB of chromosome 16p, using subjects with
well-characterized deletions in the region and healthy controls. But it
has not been validated with whole-genome chips or other DNA sources. (If
you want to play around with the code, the sources are available from
http://www.well.ox.ac.uk/~tprice/cgh/ and as the R package
\code{cgh}). (I thank Tom Price for his discussion on these issues; some
sentences here are copied almost verbatim from his emails.)

}


\value{
  All the functions return an object of class "adacgh.generic.out". In
  addition, depending on the function, the returned object can also be
  of class "CGH.wave.merged","CGH.wave", "mergedHMM", "adacghGLAD", "mergedBioHMM",
  "CGHseg", "CGH.PSW".

  
  For methods DNAcopy, HMM, BioHMM, GLAD, CGHseg, and wavelet smoothing
  the output object is a list with two components:
  \item{segm}{The output from the segmentation}
  \item{chrom.numeric}{The same value as chrom.numeric in the function
    call.}


  The \code{segm} component is a list, with as many components as
  subjects or arrays. Each of this lists is a matrix with three
  columns. The first are the observed data, the second the
  predicted/smoothed/merged data, and the third the state/alteration. When
  the method performs a merging
  step, we call the third column "Alteration", and it can only take
  three values, -1 (loss of genomic DNA), 0, +1 (gain of genomic DNA).
  

  Method PSW's output contains two components. The second,
  \code{plotData} with data that is used for plotting (and the details
  are not docummented, since you are not supposed to use them, as this
  is subject to change). The first component, \code{Data}, is a list
  with as many components as arrays or subjects. Each of this is an
  array with columns with original data, the sign, the value from the
  robustness test, and the p-value from the permutation test. See
  details in  \code{\link[cgh]{sw.threshold}},
  \code{\link[cgh]{sw.perm.test}}, \code{\link[cgh]{sw.rob}}.


  Method ACE's output is a list. None of the elements of this list are
  to be used directly by a user, and are later used by the function
  \code{\link{summary.ACE}}.
  
}






\references{


Fridlyand, Jane   and Snijders, Antoine  M.  and Pinkel, Dan   and
Albertson, Donna  G. (2004). Hidden Markov models approach to the
analysis of array CGH data. \emph{Journal of Multivariate Analysis},
\bold{90}: 132--153. 
  
  Hsu L, Self SG, Grove D, Randolph T, Wang K, Delrow JJ, Loo L, Porter
P. (2005) Denoising array-based comparative genomic hybridization data
using wavelets. \emph{Biostatistics}, \bold{6}:211-26.


Hupe, P.  and Stransky, N.  and Thiery, J. P.  and Radvanyi, F.  and
Barillot, E. (2004). Analysis of array CGH data: from signal ratio to
gain and loss of DNA regions. \emph{Bioinformatics}, \bold{20}:
3413--3422. 


  Lingjærde OC, Baumbusch LO, Liestøl K, Glad I,
Børresen-Dale AL. (2005). CGH-Explorer: a program for analysis of
CGH-data. \emph{Bioinformatics}, \bold{21}: 821--822.



Marioni, J. C.  and Thorne, N. P.  and Tavare, S. (2006). BioHMM: a
heterogeneous hidden Markov model for segmenting array CGH
data. \emph{Bioinformatics}, \bold{22}: 1144--1146.





  Olshen, A. B.  and Venkatraman, E. S.  and Lucito, R.  and Wigler,
  M. (2004) Circular binary segmentation for the analysis of array-based
  DNA copy number data. \emph{Biostatistics}, \bold{4}, 557--572.
  \url{http://www.mskcc.org/biostat/~olshena/research}.


  Picard, F.  and Robin, S.  and Lavielle, M.  and Vaisse, C.  and
  Daudin, J. J. (2005). A statistical approach for array CGH data
  analysis. \emph{BMC Bioinformatics}, \bold{6},
  27. \url{http://dx.doi.org/10.1186/1471-2105-6-27}. 


  
  Price TS, Regan R, Mott R, Hedman A, Honey B, Daniels RJ,
Smith L, Greenfield A, Tiganescu A, Buckle V, Ventress N, Ayyub H,
Salhan A, Pedraza-Diaz S, Broxholme J, Ragoussis J, Higgs DR, Flint J,
Knight SJ. (2005) SW-ARRAY: a dynamic programming solution for the
identification of copy-number changes in genomic DNA using array
comparative genome hybridization data. \emph{Nucleic Acids Res.},
\bold{33}:3455-64. 


  Willenbrock, H. and Fridlyand, J. (2005). A comparison study: applying
  segmentation to array CGH data for downstream
  analyses. \emph{Bioinformatics}, \bold{21}, 4084--4091.
  

  
  Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
  web-based application for the analysis of aCGH data. Tech. report.
  \url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
  \url{http://adacgh.bioinfo.cnio.es}.

}





\author{

The code for DNAcopy, HMM, BioHMM, PSW, and GLAD are basically
wrappers around the original functions by their corresponding authors,
with some modiffications for parallelization (packages
\code{DNAcopy}, \code{aCGH}, \code{snapCGH}, \code{cgh}, \code{GLAD},
respectively). The CGHseg method uses package \code{tilingArray}.


For the wavelet-based method we have only wrapped the code that was
kindly provided by L. Hsu and D. Grove, and parallelized a few calls. Their
original code is included in the sources of the package.

The code for ACE is based on the original Java code for the CGHExplorer
package. It was turned into C and R code by Oscar M. Rueda \email{omrueda@cnio.es}.

Parallelization and other modifications and additions are by Ramon Diaz-Uriarte \email{rdiaz02@gmail.com}

}

 \seealso{
%   \code{\link[DNAcopy]{segment}}
%   \code{\link[DNAcopy]{plot.DNAcopy}}
   \code{\link{segmentPlot}}}
 }


 \examples{
\dontrun{
setwd("/tmp/o3") ## all slaves need a common dir to read and write.
mpiInit()}

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

## need to choose fdr
ace.out.sum <- summary(ace.out)

#### Writing
common <- cghE1[, -c(5:7)]
writeResults(hmm.out, cghE1[, 5:7], common)
writeResults(glad.out, cghE1[, 5:7], common)        
writeResults(cghseg.out, cghE1[, 5:7], common)   
writeResults(ace.out.sum, cghE1[, 5:7], common)
writeResults(wave.out, cghE1[, 5:7], common)
writeResults(wave.nm.out, cghE1[, 5:7], common)
writeResults(cbs.out, cghE1[, 5:7], common)
writeResults(psw.pos.out, cghE1[, 5:7], common)
writeResults(psw.neg.out, cghE1[, 5:7], common)
writeResults(biohmm.out, cghE1[, 5:6], common)



########### Comparing our parallel output with original sequential in DNAcopy

data(cghMCRe)
common <- cghMCRe[, -c(5:7)]
common$MidPoint <- common$Start + 0.5 * (common$End - common$Start)
colnames(common)[1] <- "ID"
chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24
reorder <- order(chrom.numeric,
                 common$MidPoint,
                 cghMCRe$Start,
                 cghMCRe$End,
                 cghMCRe$Name)
cghMCRe <- cghMCRe[reorder, ]
chrom.numeric <- chrom.numeric[reorder]

CNA.object <- CNA(as.matrix(cghMCRe[, 5:7]),
                  chrom = chrom.numeric,
                  maploc = 1:nrow(cghMCRe),
                  data.type = "logratio",
                  sampleid = colnames(cghMCRe[, 5:7]))
smoothed.CNA.object <- smooth.CNA(CNA.object)

out.merge <- pSegmentDNAcopy(smoothed.CNA.object)

out1 <- pSegmentDNAcopy(smoothed.CNA.object, merge = FALSE,
                        nperm = 10000)

## Yes, tilingArray has an identically named function, but DNAcopy does
## not have a namespace. Problems would be caught by different params.

out.original <- segment(smoothed.CNA.object, nperm = 10000,
                        undo.splits = "prune")


## We need to stretch the output of the original function

stretchCNAoutput <- function(object) {
    if(!(inherits(object, "DNAcopy")))
        stop("This function can only be applied to DNAcopy objects")
    numarrays <- ncol(object$data) - 2
    stretched <- list()
    for(arraynum in 1:numarrays) {
        obs <- object$data[, 2 + arraynum]
        segmented <-
            object$output[object$output$ID ==
                          colnames(object$data)[2 + arraynum], ]
        smoothed <- object$data$maploc 
        for(i in 1:nrow(segmented)) {
            smoothed[(segmented[i,'loc.end'] >= smoothed) &
                     (segmented[i,'loc.start'] <= smoothed)] <-
                         segmented[i,'seg.mean']
        }
        
        stretched[[arraynum]] <- cbind(Observed = obs,
                                       Predicted = smoothed)
        
    }
    return(stretched)
}


## recall the originall segment rounds output
verifCNA <- mapply(function(x, y) all.equal(x[, 2], round(y[, 2], 4)),
                    stretchCNAoutput(out.original), out1$segm) 

if(! all(verifCNA)) stop("Differences in original and parallelized DNAcopy.\n Rerun test")




## verify reconstruction
rec.verif <- function(x, y) {
    t1 <- all.equal(x$output[, 3], y$output$loc.start)
    t2 <- all.equal(x$output[, 4], y$output$loc.end)
    t3 <- all.equal(x$output[, 6], round(y$output$seg.mean, 4))

    if(t1) print("Start equal")
    if(t2) print("End equal")
    if(t3) print("seg.mean equal") else t3
    if(t1 & t2 & t3) print("***Global result:  OK ***") else stop("Not equal")
}
reconstructed <- ADaCGH:::constructSegmObj(out1$segm,  chrom.numeric, cghMCRe[, 5:7],
                                  Pos = 1:nrow(cghMCRe))

rec.verif(out.original, reconstructed)




}


\keyword{ nonparametric }

