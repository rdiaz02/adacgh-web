\name{pSegment}
\alias{pSegmentDNAcopy}
\alias{pSegmentHMM}
\alias{pSegmentGLAD}
\alias{pSegmentBioHMM}
\alias{pSegmentCGHseg}
\alias{pSegmentWavelets}
\alias{pSegmentHaarSeg}
\alias{pSegment}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized/"unified" versions of several aCGH segementation algorithms/methods}

\description{
  These functions parallelize several segmentation algorithms or (for
  HaarSeg) make their calling use the same conventions as for other methods.
 }


\usage{

pSegmentDNAcopy(cghRDataName, chromRDataName, merging = "mergeLevels",
                mad.threshold = 3, smooth = TRUE,
                alpha=0.01, nperm=10000,
                p.method = "hybrid",
                min.width = 2,
                kmax=25, nmin=200,
                eta = 0.05,trim = 0.025,
                undo.splits = "none",
                undo.prune=0.05, undo.SD=3,
                ...)

pSegmentHaarSeg(cghRDataName, chromRDataName,
                merging = "MAD", mad.threshold = 3, 
                W = vector(),
                rawI = vector(), 
                breaksFdrQ = 0.001,			  
                haarStartLevel = 1,
                haarEndLevel = 5, ...)

pSegmentHMM(cghRDataName, chromRDataName,
            merging = "mergeLevels", mad.threshold = 3,
            aic.or.bic == "AIC", ...)



pSegmentBioHMM(cghRDataName, chromRDataName, PosRDataName,
               merging = "mergeLevels", mad.threshold = 3,
               aic.or.bic = "AIC",
               ...)

pSegmentCGHseg(cghRDataName, chromRDataName, CGHseg.thres = -0.05,
               merging = "MAD", mad.threshold = 3, ...)


pSegmentGLAD(cghRDataName, chromRDataName,
             deltaN = 0.10,
             forceGL = c(-0.15, 0.15),
             deletion = -5,
             amplicon = 1,
             ...)


pSegmentWavelets(cghRDataName, chromRDataName, merging = "MAD",
                 mad.threshold = 3,
                 minDiff = 0.25,
                 minMergeDiff = 0.05,
                 thrLvl = 3, initClusterLevels = 10, ...)




}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{cghRDataName}{The Rdata file name that contains the
    \code{\link[ff]{ffdf}} with the aCGH data. This file can be created
    using \code{\link[ff]{as.ffdf}} with a data frame with genes
    (probes) in rows and subjects or arrays in columns. Function
    \code{\link{inputDataToADaCGHData}} produces these type of files.}
  
  
  \item{chromRDataName}{The RData file name with the ff (short integer)
    vector with the chromosome indicator. Function
    \code{\link{inputDataToADaCGHData}} produces these type of files. }
  
  \item{PosRDataName}{The RData file name with the ff double vector with
    the location (e.g., position in kbases) of each probe in
    the chromosome. Function
    \code{\link{inputDataToADaCGHData}} produces these type of files.}
  
  \item{merging}{Merging method; for most methods one of "MAD" or
    "mergeLevels". For CBS (pSegmentDNAcopy) also "none". This options
    does not apply to GLAD (which has its own merging-like approach). See
    details.}
  
  \item{mad.threshold}{If using \code{merging = "MAD"} the value such
    that all segments where abs(smoothed value) > m*MAD will be declared
    aberrant ---see p. i141 of Ben-Yaacov and Eldar. No effect if
    merging = "mergeLevels" (or "none").}

  

  \item{smooth}{For DNAcopy only. If TRUE (default) carry out smoothing as
    explained in \code{\link{smooth.CNA}}.}

  \item{alpha}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{nperm}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{p.method}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{min.width}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{kmax}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{eta}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{trim}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{undo.splits}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{undo.prune}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 
  \item{undo.SD}{For DNAcopy only. See \code{\link[DNAcopy]{segment}}.} 


  \item{W}{For HaarSeg:  Weight matrix, corresponding to quality of measurment. 
    Insert 1/(sigma**2) as weights if your platform output sigma as
    the quality of measurment. W must have the same size as I.}
  \item{rawI}{For HaarSeg. Mininum of the raw red and raw green measurment, before the log.
    rawI is used for the non-stationary variance compensation. 
    rawI must have the same size as I. }
  \item{breaksFdrQ}{For HaarSeg.  The FDR q parameter. Common used
    values are 0.05, 0.01, 0.001. Default value is 0.001.}
  \item{haarStartLevel}{For HaarSeg.  The detail subband from which we
    start to detect peaks. The higher  this value is, the less sensitive
    we are to short segments. The default is value is 1, corresponding to
    segments of 2 probes.}
  \item{haarEndLevel}{For HaarSeg. The detail subband until which we
    use to detect peaks. The higher  this value is, the more sensitive we
    re to large trends in the data. This value DOES NOT indicate the
    largest possible segment that can be detected. The default is value is
    5, corresponding to step of 32 probes in each direction.}

  \item{aic.or.bic}{For HMM and BioHMM. One of "AIC" or "BIC". See
  \option{criteria} in \code{\link[snapCGH]{runBioHMM}}.} 


  \item{CGHseg.thres}{The threshold for the adaptive penalization in
    Picard et al.'s CGHseg. See p. 13 of the original paper. Must be a
    negative number. The default value used in the original reference
    is -0.5. However, our experience with the simulated data in
    Willenbrock and Fridlyand (2005) indicates that for those data
    values around -0.005 are more appropriate. We use here -0.05 as
    default.}


  \item{deltaN}{Only for GLAD. See \option{deltaN} in \code{\link[GLAD]{daglad}}.}
  \item{forceGL}{Only for GLAD. See \option{forceGL} in \code{\link[GLAD]{daglad}}.}
  \item{deletion}{Only for GLAD. See \option{deletion} in \code{\link[GLAD]{daglad}}.}
  \item{amplicon}{Only for GLAD. See \option{amplicon} in \code{\link[GLAD]{daglad}}.}
  


  \item{minMergeDiff}{Used only when doing merging in the wavelet
    method of Hsu et al..  The finall call as to which segments go together is done by
    a \code{mergeLevels} approach, but an initial collapsing of very
    close values is performed (otherwise, we could end up passing to
    mergeLevels as many initial levels as there are points). }

  \item{minDiff}{For Wavelets (Hsu et al.). Minimum (absolute)
    difference between the medians of two adjacent clusters for them to
    be considered truly different.  Clusters "closer" together than this
    are collapsed together to form a single cluster. }

  \item{thrLvl}{ The level used for the wavelet thresholding in Hsu et
    al.}

  \item{initClusterLevels}{For Wavelets (Hsu et al.). The initial number of clusters to form. }
 
  
  \item{...}{Additional arguments; not used.}
}



\details{

In most cases, these are wrappers to the original code, with
modifications for parallelization and for using \code{\link[ff]{ff}}
objects. The functions will not work if you try to use them with the
regular R data frames, matrices, and vectors.


HaarSeg has been turned into an R package, available from
\url{https://r-forge.r-project.org/projects/haarseg/}. That package uses
the same C code as we do. We have not used the available R package for
historical reasons (we used the C and R code in the former ADaCGH
package, before a proper R package was available).


. For CGHseg, BioHMM and HMM the
results are merged, by the mergeLevels algorithm, as recommended in
\cite{Willenbrock and Fridlyand, 2005}. Merging is also done in GLAD
(with GLAD's own merging algorithm). Merging is the default for DNAcopy
(but this can be set by the user). Merging can also be used with the
wavelet-based method of Hsu et al.; please note that the later is an
experimental feature implemented by us, and there is no evidence of its
performance. 

Using the output of mergeLevels, we then map the results to states of
  "Alteration", so that we categorize each probe as taking one, and only
  one, of three possible values, -1 (loss of genomic DNA), 0 (no change
  in DNA content), +1 (gain of genomic DNA). We have made the
  assumption, in this mapping, that the "no change" class is the one
  that has the absolute value closest to zero, and any other classes are
  either gains or losses. When the data are normalized, the "no change"
  class should be the most common one.
  
For HaarSeg, we use the approach in page i141 of Ben-Yaacov
and Eldar, section 2.3, "Determining aberrant intervals": a MAD (per
their definition) is computed and any segment with absolute value larger
than m * MAD is considered aberrant.

  
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

Beware that, with DNAcopy, parallelization over arrays by chromosome can
lead to slightly different results than parallelization over arrays
(and, thus, to the original function). This is because of the way the
smoothing is carried out, which affects the timmed.SD used.

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
  subjects or arrays. Each of these lists is a matrix with three
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


  Lingjaerde OC, Baumbusch LO, Liestol K, Glad I,
Borresen-Dale AL. (2005). CGH-Explorer: a program for analysis of
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

  Ben-Yaacov, E. and Eldar, Y.C. (2008). A Fast and Flexible Method for
the Segmentation of aCGH Data, \emph{Bioinformatics}, \bold{24}: i139-i145.

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
   \code{\link{segmentPlot}}
 }



 \examples{
\dontrun{
## all slaves need a common dir to read and write.
dir.to.create <- tempdir2()
setwd(dir.to.create)
dir.to.create
mpiInit()  ## this is used in SegmentPlotWrite, which is run first.
}

## library(Rmpi)
## try(mpi.spawn.Rslaves(nslaves = 1), silent = TRUE)

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
cbs.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric)
cbs.ns.out <- pSegmentDNAcopy(cghE1[, 5:7], chrom.numeric, smooth = FALSE)
glad.out <- pSegmentGLAD(cghE1[, 5:7], chrom.numeric)
cghseg.out <- pSegmentCGHseg(cghE1[, 5:7], chrom.numeric, CGHseg.thres = -0.05)
ace.out <- pSegmentACE(cghE1[, 5:7], chrom.numeric)

hmm.out <- pSegmentHMM(cghE1[, 5:7], chrom.numeric)
wave.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric)
wave.nm.out <- pSegmentWavelets(cghE1[, 5:7], chrom.numeric, mergeSegs = FALSE)

haar.out <- pSegmentHaarSeg(cghE1[, 5:7], chrom.numeric)

## Next onew we try, as it could fail when it should not.
try(psw.pos.out <- pSegmentPSW(cghE1[, 5:7], chrom.numeric, sign = 1))
try(psw.neg.out <- pSegmentPSW(cghE1[, 5:7], chrom.numeric, sign = -1))

## BioHMM is the only one that uses distances
## it is the slowest, so do only two
## We use try, as snapCGH sometimes fails and does not load its compiled C code
try(biohmm.out <- pSegmentBioHMM(cghE1[, 5:6], chrom.numeric, cghE1$UG.Start))

## With ACE we need to choose fdr
ace.out.sum <- summary(ace.out)

#### Writing
common <- cghE1[, 1:4]
writeResults(hmm.out, cghE1[, 5:7], common)
writeResults(glad.out, cghE1[, 5:7], common)        
writeResults(cghseg.out, cghE1[, 5:7], common)   
writeResults(ace.out.sum, cghE1[, 5:7], common)
writeResults(wave.out, cghE1[, 5:7], common)
writeResults(wave.nm.out, cghE1[, 5:7], common)
writeResults(cbs.out, cghE1[, 5:7], common)
writeResults(haar.out, cghE1[, 5:7], common)

try(writeResults(psw.pos.out, common))
try(writeResults(psw.neg.out, common))
try(writeResults(biohmm.out, cghE1[, 5:6], common))





### plot all

\dontrun{

## This will not work unless you have a working Python installation

yminmax <- c(min(as.matrix(cghE1[, 5:7])),
             max(as.matrix(cghE1[, 5:7])))

segmentPlot(hmm.out, 
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(glad.out,
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(cghseg.out,
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(wave.out, 
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(wave.out, 
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(haar.out,
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

segmentPlot(cbs.out, 
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

try(segmentPlot(psw.pos.out,
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs"))

try(segmentPlot(psw.neg.out,
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs"))

try(segmentPlot(biohmm.out, 
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs"))

segmentPlot(ace.out.sum, 
            geneNames = cghE1[, 1],
            yminmax = yminmax,
            idtype = "ug",
            organism = "Hs")

}

## To verify we are doing the right thing with Picard et al.'s method,
## compare with their figure 1. Need to get the txt files with the
## original data.
## We got the data from the GLAD BioC package, and then exported the txt.
## For convenience, it is in the "/inst/txt-data" directory of our package.

## This verifies the output (look at ADaCGH.R for further details and other
## comparisons):

## coriel.data <- read.table("gm03563.txt", header = TRUE)
## cd3 <- coriel.data[coriel.data$Chromosome == 3, 3]
    
## out.our <- CGHsegWrapper(cd3, rep(1, length(cd3)), s = -0.5)  ## k = 2
## out.our

## so our implementation seems correct. Look at where the
## breakpoint is located, etc, and it is like figure 1 of
## Picard's paper.




}


\keyword{ nonparametric }

