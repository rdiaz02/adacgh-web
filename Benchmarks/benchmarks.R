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

## library(DNAcopy)
## data(coriell)


## ## Prepare data
## cordat <- na.omit(unlist(coriell[, c(4, 5)]))
## names(cordat) <- NULL
## cordat <- matrix(c(cordat, sample(cordat, 353)), ncol = 2)

## dataList <- list()
## dataList[[1]] <- list()
## dataList[[1]]$acgh <- cordat
## dataList[[1]]$nr <- 2271
## dataList[[1]]$nc <- 2
## dataList[[1]]$Chrom <- coriell$Chromosome
## dataList[[1]]$Pos <- coriell$Position

## load("long.data.RData")


## ## give ordered data
## start2 <- unlist(tapply(posdat2$start, posdat2$chrom, sort))

## ## verify ordered
## cucu <- tapply(start2, posdat2$chrom, order)
## lapply(cucu, function(z) max(abs(z - (1:length(z)))))

## dataList[[2]] <- list()
## dataList[[2]]$acgh <- long.long
## dataList[[2]]$nr <- 42325
## dataList[[2]]$nc <- 11
## dataList[[2]]$Chrom <- posdat2$chrom
## dataList[[2]]$Pos <- start2


## selected <- sort(sample(1:42325, 10000))

## dataList[[3]] <- list()
## dataList[[3]]$acgh <- dataList[[2]]$acgh[selected,]
## dataList[[3]]$nr <- 10000
## dataList[[3]]$nc <- 11
## dataList[[3]]$Chrom <- dataList[[2]]$Chrom[selected]
## dataList[[3]]$Pos <- dataList[[2]]$Pos[selected]

## save(file = "dataList.RData", dataList)


## nohup R --no-save --no-restore --slave < benchmarks.R > benchmarksB.Rout &


##################################################33


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
    segmented <- segment(smoothed, undo.splits = "prune", nperm = 10000)
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
    
    gcnum <- 5
    
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
    
}




#numberArrays <- c(5, 10, 20, 40, 80, 120)
#mpiSizes     <- c(1, 4, 20, 60, 120)
#reps <- 5

numberArrays <- c(5, 10, 20, 50, 100, 150)
mpiSizes     <- c(10, 30, 60, 120)
reps <- 10
Methods <- c("GLAD", "HMM", "CBS")

smallTiming120 <- getTimesB(1, nsamps = numberArrays,
                            reps = reps,
                            MPI = 120, Methods = Methods)
save(file = "smallTiming120_10reps.RData", smallTiming120)


smallTiming60 <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = 60, Methods = Methods)
save(file = "smallTiming60_10reps.RData", smallTiming60)


smallTiming30 <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = 30, Methods = Methods)
save(file = "smallTiming30_10reps.RData", smallTiming30)



smallTimingNone <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = "None", Methods = Methods)
save(file = "smallTimingNone_10reps.RData", smallTimingNone)



smallTiming10 <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = 10, Methods = Methods)
save(file = "smallTiming10_10reps.RData", smallTiming10)




mediumTiming120 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 120, Methods = Methods)
save(file = "mediumTiming120_10reps.RData", mediumTiming120)


mediumTiming60 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 60, Methods = Methods)
save(file = "mediumTiming60_10reps.RData", mediumTiming60)


mediumTiming30 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 30, Methods = Methods)
save(file = "mediumTiming30_10reps.RData", mediumTiming30)


mediumTimingNone <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = "None", Methods = Methods)
save(file = "mediumTimingNone_10reps.RData", mediumTimingNone)

mediumTiming10 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 10, Methods = Methods)
save(file = "mediumTiming10_10reps.RData", mediumTiming10)



















### Now only methods that have no sequential counterpart




numberArrays <- c(5, 10, 20, 50)
mpiSizes     <- c(10, 30, 60, 120)
reps <- 3
MethodsP <- c("ACE", "CGHseg", "PSW", "Wave")

## Already done
smallTimingNoSeq120 <- getTimesB(1, nsamps = numberArrays,
                            reps = reps,
                            MPI = 120, Methods = MethodsP)
save(file = "smallTimingNoSeq120_10reps.RData", smallTimingNoSeq120)


smallTimingNoSeq60 <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = 60, Methods = MethodsP)
save(file = "smallTimingNoSeq60_10reps.RData", smallTimingNoSeq60)


smallTimingNoSeq30 <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = 30, Methods = MethodsP)
save(file = "smallTimingNoSeq30_10reps.RData", smallTimingNoSeq30)



smallTimingNoSeq10 <- getTimesB(1, nsamps = numberArrays,
                        reps = reps,
                        MPI = 10, Methods = MethodsP)
save(file = "smallTimingNoSeq10_10reps.RData", smallTimingNoSeq10)



mediumTimingNoSeq120 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 120, Methods = MethodsP)
save(file = "mediumTimingNoSeq120_10reps.RData", mediumTimingNoSeq120)


mediumTimingNoSeq60 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 60, Methods = MethodsP)
save(file = "mediumTimingNoSeq60_10reps.RData", mediumTimingNoSeq60)


mediumTimingNoSeq30 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 30, Methods = MethodsP)
save(file = "mediumTimingNoSeq30_10reps.RData", mediumTimingNoSeq30)



mediumTimingNoSeq10 <- getTimesB(3, nsamps = numberArrays,
                        reps = reps,
                        MPI = 10, Methods = MethodsP)
save(file = "mediumTimingNoSeq10_10reps.RData", mediumTimingNoSeq10)
















largeTiming120 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 120, Methods = Methods)
save(file = "largeTiming120_10reps.RData", largeTiming120)


largeTiming60 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 60, Methods = Methods)
save(file = "largeTiming60_10reps.RData", largeTiming60)


largeTiming30 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 30, Methods = Methods)
save(file = "largeTiming30_10reps.RData", largeTiming30)


largeTimingNone <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = "None", Methods = Methods)
save(file = "largeTimingNone_10reps.RData", largeTimingNone)

largeTiming10 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 10, Methods = Methods)
save(file = "largeTiming10_10reps.RData", largeTiming10)









largeTimingNoSeq120 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 120, Methods = MethodsP)
save(file = "largeTimingNoSeq120_10reps.RData", largeTimingNoSeq120)


largeTimingNoSeq60 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 60, Methods = MethodsP)
save(file = "largeTimingNoSeq60_10reps.RData", largeTimingNoSeq60)


largeTimingNoSeq30 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 30, Methods = MethodsP)
save(file = "largeTimingNoSeq30_10reps.RData", largeTimingNoSeq30)



largeTimingNoSeq10 <- getTimesB(2, nsamps = numberArrays,
                        reps = reps,
                        MPI = 10, Methods = MethodsP)
save(file = "largeTimingNoSeq10_10reps.RData", largeTimingNoSeq10)





