## selected <- sort(sample(1:42325, 20000))
## dataList[[4]] <- list()
## dataList[[4]]$acgh <- dataList[[2]]$acgh[selected,]
## dataList[[4]]$nr <- 20000
## dataList[[4]]$nc <- 11
## dataList[[4]]$Chrom <- dataList[[2]]$Chrom[selected]
## dataList[[4]]$Pos <- dataList[[2]]$Pos[selected]
## save(file = "dataList.RData", dataList)



rm(list = ls())
load("dataList.RData")
library(Rmpi)
library(papply)
library(ADaCGH)
library(DNAcopy) ## try not to use segment from tilingArray


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
        dd <- tmp[ , inds]
    } else {
        inds <- sample(1:nc, (samps - nc), replace = TRUE)
        rans <- matrix(rnorm(nr * (samps - nc)), ncol = (samps - nc))
        dd <- cbind(tmp, tmp[, inds] + rans)
    }
    save.image(file = "after.make.cghdata.RData")
    save(file = "dd.RData", dd)
    return(dd)
}

## the next two, needed to feed the functions the appropriate
## type of object


make.cbs.obj <- function(ind, samps) {
    tmp <- as.matrix(make.cghdata(ind, samps))
    cna.obj <- CNA(tmp, dataList[[ind]]$Chrom,
                   1:nrow(tmp),
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
    segmented <- segment(smoothed, undo.splits = "none", nperm = 10000)
}

doHMM <- function(data, nsamp) {
    res <- find.hmm.states(data, aic = TRUE, bic = FALSE)
    hmm(data) <- res
    m <- matrix(NA, nrow = nrow(data$hmm$states.hmm[[1]]),
                ncol = nsamp)
    for(j in 1:ncol(data)) {
        m[, j] <- mergeLevels(data$hmm$states.hmm[[1]][, 2 + (6 * j)],
                              data$hmm$states.hmm[[1]][, 2 + (6 * j) - 2])$vecMerged
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
    tt <- try(unix.time(trash <- doBioHMM(data, chr, pos))[3])
    if(inherits(tt, "try-error")) {
        return(NA)
    } else {
        return(tt)
    }
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
    return(unix.time(trash <- pSegmentDNAcopy(data, chr))[3])
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

    
    Method <- c("GLAD", "HMM", "CBS") ## BioHMM excluded because of its crashes
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
        browser()
        method <- as.character(method)
        mpi <- as.character(mpi)
        cat("\n\n Doing case number ", number,
            ". Method = ", method, ". MPI = ", mpi,
            ". Number of samples = ", nsamp, "\n\n")

        for(gci in 1:gcnum) print(gc())
        
        if(mpi == "None") {
            return(do.call(paste("f", method, sep = ""),
                           list(ind = dataind, samples = nsamp)))
        } else {
            mpiSetup(as.numeric(mpi))
            return(do.call(paste("f", method, ".mpi", sep = ""),
                           list(ind = dataind, samples = nsamp)))
            
        }
    }      

    out <- mapply(f1, designm$NumberArrays, designm$Method,
                  designm$MPI, designm$number)
    out <- cbind(designm, out)
    
}







getTimesB <- function(dataind,
                      nsamps,
                      reps,
                      Methods,
                      MPI) {
    ## starting and stopping MPI repeatedly leads to problems
    ## we call gc() a number of times
    ## like getTimes, but methods are passed as parameters, which allows
    ## finer contol if need to repeat.
    save.image(file = "upToHere.RData") ## up to here
    gcnum <- 1
    
     designm <- expand.grid(Rep = 1:reps,
                           NumberArrays = nsamps,
                           Method = Methods,
                           MPI    = MPI)
    if(MPI != "None") {
        mpiSetup(MPI)
    }
    
    ## Make sure no order effects
    nd <- nrow(designm)
    designm <- designm[sample(1:nd), ]
    designm$number <- 1:nd

    cat("\n\n **** getTimes will need to run for ",
        nd, "  times. **** \n\n")

    f1 <- function(nsamp, method, mpi, number) {
        method <- as.character(method)
        mpi <- as.character(mpi)
        cat("\n\n Doing case number ", number,
            ". Method = ", method, ". MPI = ", mpi,
            ". Number of samples = ", nsamp, "\n\n")

        for(gci in 1:gcnum) print(gc())
        
        if(mpi == "None") {
            return(do.call(paste("f", method, sep = ""),
                           list(ind = dataind, samples = nsamp)))
        } else {
            return(do.call(paste("f", method, ".mpi", sep = ""),
                           list(ind = dataind, samples = nsamp)))
            
        }
    }      

    out <- mapply(f1, designm$NumberArrays, designm$Method,
                  designm$MPI, designm$number)
    out <- cbind(designm, out)
    out
}


fHMM_A.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentHMM_A(data, chr))[3])
}

fCBS_A.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentDNAcopy_A(data, chr))[3])
}

fBioHMM_A.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    pos <- dataList[[ind]]$Pos
    return(unix.time(trash <- pSegmentBioHMM_A(data, chr, pos))[3])
}    

fACE_A.mpi <- function(ind, samples) {
    data <- make.cghdata(ind, samples)
    chr <- dataList[[ind]]$Chrom
    return(unix.time(trash <- pSegmentACE_A(data, chr))[3])
}


numberArrays <- c(10, 20, 50, 100, 150)
mpiSizes     <- c(60)
reps <- 4

MethodsA <- c("GLAD", "HMM", "CBS", "BioHMM")
MethodsB <- c("ACE", "CGHseg", "Wave")

### Recall: in dataList,
##          4: 20000 genes: m
##          3: 10000 genes: s
##          2: 42325 genes: b



#### We first run the parallelized. The sequential ones can be run
#### several at a time from different machines as there is no interferecne.

#### Sure, we could write a function. But this allows for
##   simpler error recovery if something  breaks.



#### First those with a sequential counterpart

smallTiming60 <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = 60, Methods = MethodsA)
save(file = "smallTiming60.RData", smallTiming60)

smallTiming30 <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = 30, Methods = MethodsA)
save(file = "smallTiming30.RData", smallTiming30)

smallTiming10 <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = 10, Methods = MethodsA)
save(file = "smallTiming10.RData", smallTiming10)



mediumTiming60 <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = 60, Methods = MethodsA)
save(file = "mediumTiming60.RData", mediumTiming60)

mediumTiming30 <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = 30, Methods = MethodsA)
save(file = "mediumTiming30.RData", mediumTiming30)

mediumTiming10 <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = 10, Methods = MethodsA)
save(file = "mediumTiming10.RData", mediumTiming10)




largeTiming60 <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = 60, Methods = MethodsA)
save(file = "largeTiming60.RData", largeTiming60)

largeTiming30 <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = 30, Methods = MethodsA)
save(file = "largeTiming30.RData", largeTiming30)

largeTiming10 <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = 10, Methods = MethodsA)
save(file = "largeTiming10.RData", largeTiming10)




#### Now those without sequential counterpart


smallTimingNoSeq60 <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = 60, Methods = MethodsB)
save(file = "smallTimingNoSeq60.RData", smallTimingNoSeq60)

smallTimingNoSeq30 <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = 30, Methods = MethodsB)
save(file = "smallTimingNoSeq30.RData", smallTimingNoSeq30)

smallTimingNoSeq10 <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = 10, Methods = MethodsB)
save(file = "smallTimingNoSeq10.RData", smallTimingNoSeq10)



mediumTimingNoSeq60 <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = 60, Methods = MethodsB)
save(file = "mediumTimingNoSeq60.RData", mediumTimingNoSeq60)

mediumTimingNoSeq30 <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = 30, Methods = MethodsB)
save(file = "mediumTimingNoSeq30.RData", mediumTimingNoSeq30)

mediumTimingNoSeq10 <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = 10, Methods = MethodsB)
save(file = "mediumTimingNoSeq10.RData", mediumTimingNoSeq10)




largeTimingNoSeq60 <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = 60, Methods = MethodsB)
save(file = "largeTimingNoSeq60.RData", largeTimingNoSeq60)

largeTimingNoSeq30 <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = 30, Methods = MethodsB)
save(file = "largeTimingNoSeq30.RData", largeTimingNoSeq30)

largeTimingNoSeq10 <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = 10, Methods = MethodsB)
save(file = "largeTimingNoSeq10.RData", largeTimingNoSeq10)



### Now the sequential part


smallTimingNone <- getTimesB(3, nsamps = numberArrays,
                           reps = reps,
                           MPI = "None", Methods = MethodsA)
save(file = "smallTimingNone.RData", smallTimingNone)

mediumTimingNone <- getTimesB(4, nsamps = numberArrays,
                           reps = reps,
                           MPI = "None", Methods = MethodsA)
save(file = "mediumTimingNone.RData", mediumTimingNone)

largeTimingNone <- getTimesB(2, nsamps = numberArrays,
                           reps = reps,
                           MPI = "None", Methods = MethodsA)
save(file = "largeTimingNone.RData", largeTimingNone)
