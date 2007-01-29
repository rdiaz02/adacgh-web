#### Benchmarking: makes no sense for ACE or Wavelets, since no other
##   implementation.

##   For PSW, we are doing different things (some anal. are over genome),
##   so hard to compare.

##   With Piccard, we are also doing different things wrt to Huber
##   (getting best k).


## Wait: with all PSW, wavelets, and Piccard can show timings, with MPI,
## even if no comparison possible

##  Therefore, compare: DNAcopy, Hmm, BioHMM, GLAD.


#### Prepare data; only needed once

library(DNAcopy)
data(coriell)


## Prepare data
cordat <- na.omit(unlist(coriell[, c(4, 5)]))
names(cordat) <- NULL
cordat <- matrix(c(cordat, sample(cordat, 353)), ncol = 2)

dataList <- list()
dataList[[1]] <- list()
dataList[[1]]$acgh <- cordat
dataList[[1]]$nr <- 2271
dataList[[1]]$nc <- 2
dataList[[1]]$Chrom <- coriell$Chromosome
dataList[[1]]$Pos <- coriell$Positions

load("long.data.RData")


## give ordered data
start2 <- unlist(tapply(posdat2$start, posdat2$chrom, sort))

## verify ordered
cucu <- tapply(start2, posdat2$chrom, order)
lapply(cucu, function(z) max(abs(z - (1:length(z)))))

dataList[[2]] <- list()
dataList[[2]]$acgh <- long.long
dataList[[2]]$nr <- 42325
dataList[[2]]$nc <- 11
dataList[[2]]$Chrom <- posdat2$chrom
dataList[[2]]$Pos <- start2

save(file = "dataList.RData", dataList)





##################################################33


rm(list = ls())
load("dataList.RData")
library(Rmpi)
library(papply)
library(ADaCGH)



mpiSetup <- function(nslaves) {
    ## setup a clean Rmpi universe
    ## (note: we make no changes to the underlying LAM/MPI)
    try(mpi.close.Rslaves())
    mpi.spawn.Rslaves(nslaves= nslaves)
    mpi.setup.rngstream() 
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(ADaCGH))
}


make.cghdata <- function(ind, samps) {
    ## data live in a list with as many elements
    ## as data sets.
    ## each list component has a acgh matrix, and vectors
    ## nr (number of rows, genes), nc (columns or arrays)
    ## chromosome and position

    ## this function returns a new acgh data set with
    ## the specified columns (samps). The data are either
    ## samples from original or that + random generation.
    
    nr <- dataList[[ind]]$nr
    nc <- dataList[[ind]]$nc
    tmp <- dataList[[ind]]$acgh
    if(samps <= nc) {
        inds <- sample(1:nc, samps, replace = TRUE)
        return(tmp[ , inds])
    } else {
        inds <- sample(1:nc, (samps - nc), replace = TRUE)
        rans <- matrix(rnorm(nr * (samps - nc)), ncol = (samps - nc))
        return(cbind(tmp, tmp[, inds] + rans))
    }
}

## the next two, needed to feed the functions the appropriate
## type of object

make.cbs.obj <- function(ind, samps) {
    tmp <- make.cghdata(ind, samps)
    cna.obj <- CNA(tmp, dataList[[ind]]$Chrom,
                         dataList[[ind]]$Pos,
                         data.type = "logratio")
    return(cna.obj)
}

make.hmm.obj <- function(ind, samps) {
    Clone <- 1:dataList[[ind]]$nr
    tmp <- make.cghdata(ind, samps)
    obj <- create.aCGH(data.frame(tmp),
                       data.frame(Clone = Clone,
                                  Chrom = dataList[[ind]]$Chrom,
                                  kb     = dataList[[ind]]$Pos)
                       )
    return(obj)
}


## dowhatever: run the nonMPI method whatever


doCBS <- function(data) {
    smoothed <- smooth.CNA(data)
    segmented <- segment(smoothed, undo.splits = "prune", nperm = 10000)
}

doHMM <- function(data, nsamp) {
    res <- find.hmm.states(data, aic = TRUE, bic = FALSE)
    hmm(data) <- res
    m <- matrix(NA, nrow = nrow(data$hmm$states.hmm[[1]]),
                ncol = nsamp)
    for(j in 1:ncol(data)) {
        m[, j] <- mergeLevels(data$hmm$states.hmm[[1]][, 2 + (6 * n)],
                              data$hmm$states.hmm[[1]][, 2 + (6 * n) - 2])$vecMerged
    }
}
        
doBioHMM <- function(data, Chrom, Pos) {
    uchrom <- unique(Chrom)
    m <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    for (j in 1:ncol(data)) {
        dataj <- data[, j]
        smoothed <- NULL
        
        for (ic in uchrom) {
            ydat <- dataj[Chrom == ic]
            n <- length(ydat)
            res <- fit.model(sample = 1, chrom = ic, dat = matrix(ydat, ncol = 1),
                             datainfo = data.frame(Name = 1:n, Chrom = rep(ic, n),
                             Position = Pos[Chrom == ic]))
            
            smoothed <- c(smoothed, res$out.list$mean)
        }
        m[, j] <- mergeLevels(data[, j], smoothed)$vecMerged
    }
}


doGLAD <- function(data, Chrom, Pos) {
    m <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    Pos <- if (is.null(Pos)) (1:nrow(data)) else Pos
    for(j in 1:ncol(data)) {
        x <- data[, j]
        tmpf <- data.frame(LogRatio = x,
                           PosOrder = Pos,
                           Chromosome = Chrom)
        tmpf <- list(profileValues = tmpf)
        class(tmpf) <- "profileCGH"
        outglad <- glad.profileCGH(tmpf)
        m[, j] <- outglad$profileValues$Smoothing
    }    
}


## fwhatever: return the time from running method whatever
##  on data set "ind" with "samples" number of arrays.
##  Note: data generation not included in timing
##        (which gives slight advantage to some non-mpi
##         implementations, that have cumbrsome data
##         preparation setups).

fHMM <- function(ind, samples) {
    data <- make.hmm.obj(ind, samples)
    return(unix.time(trash <- doHMM(data, samples))[3])
}

fCBS <- function(ind, samples) {
    data <- make.cbs.obj(ind, samples)
    return(unix.time(trash <- doCBS(data))[3])
}

fBioHMM <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    pos <- dataList[[ind]]$Pos
    return(unix.time(trash <- doBioHMM(data, chr, pos))[3])
}

fGLAD <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    pos <- dataList[[ind]]$Pos
    return(unix.time(trash <- doGLAD(data, chr, pos))[3])
}    


fHMM.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentHMM(data, chr))[3])
}

fCBS.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentDNAcopy(data, chr,
                                              mergeSegs = FALSE))[3])
}

fBioHMM.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    pos <- dataList[[ind]]$Pos
    return(unix.time(trash <- pSegmentBioHMM(data, chr, pos))[3])
}    

fGLAD.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentGLAD(data, chr))[3])
}

fACE.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentACE(data, chr))[3])
}

fCGHseg.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentCGHseg(data, chr))[3])
}

fPSW.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentPSW(data, chr))[3])
}

fWave.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentWavelets(data, chr, mergeSegs = TRUE))[3])
}



getTimes <- function(dataind,
                     nsamps,
                     reps,
                     mpisizes) {
    ## we call gc() a number of times
    gcnum <- 5

    Method <- c("GLAD", "HMM", "BioHMM", "CBS")
    MPI    <- c("None", mpisizes)
    method2 <- c("ACE", "CGHseg", "PSW", "Wave")
    mpi2 <- as.character(mpisizes)

    designm1 <- expand.grid(Rep = 1:reps,
                           NumberArrays = nsamps,
                           Method = Method,
                           MPI    = MPI)
    designm2 <- expand.grid(Rep = 1:reps,
                            NumberArrays = nsamps,
                            Method = method2,
                            MPI    = mpi2)
    designm <- rbind(designm1, designm2)
    ## Make sure no order effects
    nd <- nrow(designm)
    designm <- designm[sample(1:nd), ]
    designm$number <- 1:nd

    cat("\n\n **** getTimes will need to run for ",
        nd, "  times. **** \n\n")

    f1 <- function(nsamp, method, mpi, number) {
        cat("\n\n Doing case number ", number,
            ". Method = ", method, ". MPI = ", mpi,
            ". Number of samples = ", nsamp, "\n\n")

        for(gci in 1:gcnum) print(gc())
        
        if(mpi == "None") {
            return(do.call(paste("f", method, sep = ""),
                           list(ind = dataind, samples = nsamp)))
        } else {
            mpiSetup(as.numeric(as.character(mpi)))
            return(do.call(paste("f", method, ".mpi", sep = ""),
                           list(ind = dataind, samples = nsamp)))
            
        }
    }      

    out <- mapply(f1, designm$NumberArrays, designm$Method,
                  designm$mpi, designm$number)
    out <- cbind(designm, out)
    
}


#numberArrays <- c(5, 10, 20, 40, 80, 120)
#mpiSizes     <- c(1, 4, 20, 60, 120)
#reps <- 5

numberArrays <- c(5, 13)
mpiSizes     <- c(1, 120)
reps <- 2

largeTiming <- getTimes(2, nsamps = numberArrays,
                        reps = reps,
                        mpisizes = mpiSizes)

smallTiming <- getTimes(1, nsamps = numberArrays,
                        reps = reps,
                        mpisizes = mpiSizes)
