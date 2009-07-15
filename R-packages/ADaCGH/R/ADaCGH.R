## .__ADaCGH_WEB_APPL <- TRUE in web appl!
if(exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv))
{
  warningsForUsers <- vector()
#  running.as.web.adacgh <- TRUE
} else if (exists(".__ADaCGH_SERVER_APPL", env = .GlobalEnv)) {
  warningsForUsers <- vector()
#  running.as.web.adacgh <- FALSE
} else {
#  running.as.web.adacgh <- FALSE
  warningsForUsers <- warning
}





## BEWARE: at least with DNAcopy we use some non-user callable
## functions, such as changepoints. If things do not work, check
## arguments of functions match.


mydcat <- function(x) {
  cat("\n", x, "\n")
}


mydcat2 <- function(x) {
  cat("\n ")
  cat(deparse(substitute(x)))
  cat("\n", x, "\n")
}

mydcat3 <- function() {
  ## to understand where are arguments from, etc
  ## with papply in slaves
  cat("\n parent frame \n")
  print(parent.frame())
  cat("\n environment \n")
  print(environment())
  cat("\n parent.env \n")
  print(parent.env(environment()))

  cat("\n PID \n")
  print(Sys.getpid())
}

  



## mydcat <- function(x){}


names.formals.changepoints.1.17 <- c("genomdat",
                                     "data.type",
                                     "alpha",
                                     "sbdry",
                                     "sbn",
                                     "nperm",
                                     "p.method",
                                     "min.width",
                                     "kmax",
                                     "nmin",
                                     "trimmed.SD",
                                     "undo.splits",
                                     "undo.prune",
                                     "undo.SD",
                                     "verbose",
                                     "ngrid",
                                     "tol")
names.formals.changepoints.1.18 <- c("genomdat",
                                     "data.type",
                                     "alpha",
                                     "weights",
                                     "sbdry",
                                     "sbn",
                                     "nperm",
                                     "p.method",
                                     "min.width",
                                     "kmax",
                                     "nmin",
                                     "trimmed.SD",
                                     "undo.splits",
                                     "undo.prune",
                                     "undo.SD",
                                     "verbose",
                                     "ngrid",
                                     "tol")

## in v.1.18.0 we take advantage weights has default of NULL in changepoints

vDNAcopy <- package_version(packageDescription("DNAcopy")$Version)
if (vDNAcopy >= "1.17.1")
  {
    cat("Setting adacgh_changepoints to DNAcopy:::changepoints\n")
    adacgh_changepoints <- DNAcopy:::changepoints
    cat("Setting adacgh_trimmed.variance to DNAcopy:::trimmed.variance\n")
    adacgh_trimmed.variance <- DNAcopy:::trimmed.variance
  } else {
    cat("Setting adacgh_changepoints to changepoints\n")
    adacgh_changepoints <- changepoints
    cat("Setting adacgh_trimmed.variance to trimmed.variance\n")
    adacgh_trimmed.variance <- trimmed.variance
  }

if(vDNAcopy >= "1.18.0") {
  names.formals.changepoints <- names.formals.changepoints.1.18
} else {
  names.formals.changepoints <- names.formals.changepoints.1.17
}

if(!identical(names.formals.changepoints, names(formals(adacgh_changepoints)))) {
  m1 <- "Arguments to DNAcopy function changepoints have changed.\n"
  m2 <- "Either your version of DNAcopy is newer than ours, or older.\n"
  m3 <- "If your version is different from 1.16.0 or 1.17.1 or 1.18.0 or 1.19.0\n please let us know of this problem.\n"
  m4 <- "We are assuming you are using DNAcopy version 1.16.0 or 1.17.1 or 1.18.0 or 1.19.0 ,\n"
  m5 <- "the ones for the former  BioConductor release (v. 2.3), current stable release (v. 2.4)\n and the devel releas (v. 2.5).\n"
  m6 <- paste("Your version of DNAcopy is ", packageDescription("DNAcopy")$Version, ".\n")
  mm <- paste(m1, m2, m3, m4, m5, m6)
  stop(mm)
}


## As of v. 1.12 at least snapCGH finally has a namespace. So now we have
## to do
if(package_version(packageDescription("snapCGH")$Version) > "1.11") {
  myfit.model <- snapCGH:::fit.model
} else {
  myfit.model <- fit.model
}
## becasue even if fit.model is documented, it is NOT exported.
## Well, they get away with it because there are no executable examples
## in the help for fit.model.


## Other vars
getbdry <- DNAcopy:::getbdry



## where do we live? to call the python script
.python.toMap.py <- system.file("Python", "toMap.py", package = "ADaCGH")
    
##############################################


###  Visible stuff

mpiInit <- function(wdir = getwd(), minUniverseSize = 15,
                    universeSize = NULL, exit_on_fail = FALSE) {
    trythis <- try({
        if(! is.null(universeSize))
            minUniverseSize <- universeSize
        require(Rmpi)
### FIXME: this should be a warning, and only an error on "exit_on_fail"
        if(mpi.universe.size() < minUniverseSize) {
            if(exit_on_fail)
                stop("MPI problem: universe size < minUniverseSize")
            else
                warning("MPI problem: universe size < minUniverseSize")
        }
        ##    mpi.spawn.Rslaves(nslaves= mpi.universe.size())
        if(! is.null(universeSize)) {
            mpi.spawn.Rslaves(nslaves = universeSize)
        } else {
            mpi.spawn.Rslaves(nslaves = mpi.universe.size())
        }
        ## mpi.setup.rngstream() ## or 
        mpi.setup.sprng()
        mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
        require(papply)
##        mpi.remote.exec(require(ADaCGH, quietly = TRUE))
## Note that if ADaCGH is NOT already installed, such as in R CMD check
        ### a new, with multiple nodes, this will fail!!!!
        mpi.remote.exec(library(ADaCGH, verbose = TRUE))
        mpi.bcast.Robj2slave(wdir)
        mpi.remote.exec(setwd(wdir))    
    })
    if(inherits(trythis, "try-error")) {
        cat("\nRmpi error\n", file = "Status.msg")
        if(exit_on_fail) quit(save = "yes", status = 12, runLast = FALSE)
    }
}


pSegmentHaarSeg <- function(x, chrom.numeric, HaarSeg.m = 3,
                            W = vector(),
                            rawI = vector(), 
                            breaksFdrQ = 0.001,			  
                            haarStartLevel = 1,
                            haarEndLevel = 5, ...) {
  stop.na.inf(x)
  stop.na.inf(chrom.numeric)
  warn.too.few.in.chrom(chrom.numeric)
  rle.chr <- rle(chrom.numeric)
  chr.end <- cumsum(rle.chr$lengths)
  chr.start <- c(1, chr.end[-length(chr.end)] + 1)
  chromPos <- cbind(chr.start, chr.end)
  out <- list()
  out$segm <- list()
##  dots <- as.list(substitute(list(...)))[-1]
  
  for(subj in 1:ncol(x)) {
    cat("\n      running subject or column ", subj)
    haarout <- ad_HaarSeg(x[, subj], chromPos = chromPos,
                          W = W, rawI = rawI,
                          breaksFdrQ = breaksFdrQ,
                          haarStartLevel = haarStartLevel,
                          haarEndLevel = haarEndLevel)[[2]]
    mad.subj <- median(abs(x[, subj] - haarout))/0.6745
    thresh <- HaarSeg.m * mad.subj
    ## alteration <- rep(0, nvalues)
    ## alteration[haarout > thresh] <- 1
    ## alteration[haarout < -thresh] <- -1
    
    out$segm[[subj]] <- cbind(Observed = x[, subj],
                              Smoothed = haarout,
                              Alteration =
                              ifelse( abs(haarout > thresh), 1, 0) * sign(haarout))
  }
  cat("\n")
  out$chrom.numeric <- chrom.numeric
  out <- add.names.as.attr(out, colnames(x))
  class(out) <- c("adacgh.generic.out", "adacghHaarSeg")
  return(out)
}

pSegmentACE <- function(x, chrom.numeric, parall = "auto", ...) {
  stop.na.inf(x)
  stop.na.inf(chrom.numeric)
  warn.too.few.in.chrom(chrom.numeric)
  if (parall == "auto")
        parall <- ifelse(ncol(x) > 75, "axc", "chr")
  if (parall == "chr") {
      cat("\n    running chr version \n")
      return(ACE_C(x, chrom.numeric, echo = FALSE,
                   coefs = file.aux, Sdev = 0.2))
  }
  if (parall == "axc") {
      cat("\n    running axc version \n")
      return(ACE(x, chrom.numeric, echo = FALSE,
                 coefs = file.aux, Sdev = 0.2))
  }
}

pSegmentHMM <- function(x, chrom.numeric, parall = "auto", ...) {
  stop.na.inf(x)
  stop.na.inf(chrom.numeric)
  warn.too.few.in.chrom(chrom.numeric)

  if (parall == "auto") parall <- "arr"
  if (parall == "arr") {
      cat("\n    running arr version \n")
      return(pSegmentHMM_A(x, chrom.numeric, ...))
  }
  if (parall == "axc") {
      cat("\n    running axc version \n")
      return(pSegmentHMM_axc(x, chrom.numeric, ...))
  }
}

pSegmentGLAD <- function(x, chrom.numeric, ...) {
  stop.na.inf(x)
  stop.na.inf(chrom.numeric)
  warn.too.few.in.chrom(chrom.numeric)
  
  require("GLAD") || stop("Package not loaded: GLAD")
  out <- papply(data.frame(x),
                function(z) gladWrapper(z,
                                        Chrom = slave_chrom),
                papply_commondata = list(
                  slave_chrom = chrom.numeric))
  outl <- list()
  outl$segm <- out
  outl$chrom.numeric <- chrom.numeric
  outl <- add.names.as.attr(outl, colnames(x))
  class(outl) <- c("adacgh.generic.out", "adacghGLAD")
  return(outl)
}    


pSegmentBioHMM <- function(x, chrom.numeric, Pos, parall = "auto", ...) {
  stop.na.inf(x)
  stop.na.inf(chrom.numeric)
  stop.na.inf(Pos)
  warn.too.few.in.chrom(chrom.numeric)

    if (parall == "auto")
        parall <- ifelse(ncol(x) > 110, "arr", "axc")
    if (parall == "arr") {
        cat("\n    running arr version \n")
        return(pSegmentBioHMM_A(x, chrom.numeric, Pos, ...))
    }
    if (parall == "axc") {
        cat("\n    running axc version \n")
        return(pSegmentBioHMM_axc(x, chrom.numeric, Pos, ...))
    }
}



pSegmentHMM_axc<- function(x, chrom.numeric, ...) {
  nsample <- ncol(x)
  nchrom <- unique(chrom.numeric)
    datalist <- list()
    klist <- 1
    for(i in 1:nsample) {
        for (j in nchrom) {
            datalist[[klist]] <- x[chrom.numeric == j, i]
            klist <- klist + 1
        }
    }
    out0 <- papply(datalist, hmmWrapper)
    matout0 <- matrix(unlist(out0), ncol = nsample)
    rm(datalist)
    rm(out0)
    
    datalist <- list()
    for(i in 1:nsample) {
        datalist[[i]] <- list()
        datalist[[i]]$logr <- x[, i]
        datalist[[i]]$pred <- matout0[, i]
    }
    out <- papply(datalist, function(z) ourMerge(z$logr, z$pred))
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    outl <- add.names.as.attr(outl, colnames(x))
    class(outl) <- c("adacgh.generic.out","mergedHMM")
    return(outl)
}

    
pSegmentBioHMM_axc <- function(x, chrom.numeric, Pos, ...) {
    nsample <- ncol(x)
    nchrom <- unique(chrom.numeric)
    datalist <- list()
    klist <- 1
    for(i in 1:nsample) {
        for (j in nchrom) {
            datalist[[klist]] <- list()
            datalist[[klist]]$logr <- x[chrom.numeric == j, i]
            datalist[[klist]]$pos <- Pos[chrom.numeric == j]
            klist <- klist + 1
        }
    }
    out0 <- papply(datalist,
                  function(z) BioHMMWrapper(z$logr, Pos = z$pos))

    te <- unlist(unlist(lapply(out0, function(x) inherits(x, "try-error"))))

    if(any(te)) {
        m1 <- "The BioHMM code occassionally crashes (don't blame us!)."
        m2 <- "You can try rerunning it a few times."
        m3 <- "You can also tell the original authors that you get the error(s): \n\n "
        mm <- paste(m1, m2, m3, paste(out0[which(te)], collapse = "    \n   "))
        caughtError(mm)
    }
    matout0 <- matrix(unlist(out0), ncol = nsample)
    rm(datalist)
    rm(out0)

    datalist <- list()
    klist <- 1
    for(i in 1:nsample) {
        datalist[[i]] <- list()
        datalist[[i]]$logr <- x[, i]
        datalist[[i]]$pred <- matout0[, i]
    }
    
    out <- papply(datalist, function(z) ourMerge(z$logr, z$pred))
   
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    outl$pos <- Pos
    outl <- add.names.as.attr(outl, colnames(x))
    class(outl) <- c("adacgh.generic.out","mergedBioHMM")
    return(outl)
}



pSegmentCGHseg <- function(x, chrom.numeric, CGHseg.thres = -0.05, ...) {
  stop.na.inf(x)
  stop.na.inf(chrom.numeric)
  warn.too.few.in.chrom(chrom.numeric)

  datalist <- data.frame(x)
  out0 <- papply(datalist, function(z) CGHsegWrapper(z,
                                                     chrom.numeric = slave_chrom,
                                                     s = slave_CGHseg.thres),
                  list(slave_chrom = chrom.numeric, slave_CGHseg.thres = CGHseg.thres))
  outl <- list()
  outl$segm <- out0
  outl$chrom.numeric <- chrom.numeric
  outl <- add.names.as.attr(outl, colnames(x))
  
  class(outl) <- c("adacgh.generic.out", "CGHseg")
  return(outl)
}


pSegmentPSW <- function(x, chrom.numeric, common.data,
                        sign = -1,
                        nIter = 1000, prec = 100,  p.crit = 0.10,
                        name = NULL, ...) {
  stop.na.inf(chrom.numeric)
  stop.na.inf(x)
  warn.too.few.in.chrom(chrom.numeric)

  slname <- name
  
  numarrays <- ncol(x)
    ncrom <- length(unique(chrom.numeric))
    out <- list()
    out$plotData <- list()
  
##  if (running.as.web.adacgh) { ## send to PaLS
  if (exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) { ## send to PaLS
    print("former testing value of .AD....")
    palsVect <- vector()
    palsL <- list()
  }

  
    datalist <- list()
    for (i in 1:ncol(x)) {
        datalist[[i]] <- list()
        datalist[[i]]$data <- x[, i]
        datalist[[i]]$num  <- i
    }
  cnamesdata <- colnames(x)
    papout <- papply(datalist,
                  function(z) my.sw3b(z$data,
                                      chrom = slave_chrom.numeric,
                                      sign = slave_sign,
                                      p.crit = slave_p.crit,
                                      nIter = slave_nIter,
                                      prec = slave_prec,
                                      name = paste(slave_name, cnamesdata[z$num], sep = ""),
                                      highest = FALSE),
                      list(slave_chrom.numeric = chrom.numeric,
                           slave_sign = sign, slave_p.crit = p.crit,
                           slave_nIter = nIter,
                           slave_prec = prec,
                           slave_name = slname)
                           )
    if(any(unlist(lapply(papout, function(z) z == "swt.perm.try-error")))) {
        m1 <- "There was a problem in the PSW routine; this is \n"
        m2 <- "probably related to the global thresholding + within \n"
        m3 <- "chromosome perm test with your data.\n"
        m4 <- "You might want to try another method, or the original \n"
        m5 <- " (thresholding within chromosome) PSW. \n "
        mm <- c(m1, m2, m3, m4, m5)
        caughtOtherError(mm)
        return(NULL) ## but should not get here
    }    
    for(i in 1:ncol(x)) {
        out$Data <- cbind(out$Data, papout[[i]]$out)
        p.crit.bonferroni <- papout[[i]]$plotdat$p.crit / ncrom
        out$plotData[[i]] <- c(papout[[i]]$plotdat,
                              p.crit.bonferroni = p.crit.bonferroni)
       
        if (exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) { ## send to PaLS
            selectedGenes <-
              as.character(common.data$ID[which(papout[[i]]$out[, 4] <= p.crit.bonferroni)])
            palsVect <- c(palsVect, paste("#", colnames(x)[i], sep = ""),
                          selectedGenes)
            palsL[[i]] <- selectedGenes 
        }
    }
    class(out) <- c(class(out), "CGH.PSW")
    if (exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) { ## send to PaLS
        namef <- ifelse(sign == -1,
                        "Lost_for_PaLS.txt",
                        "Gained_for_PaLS.txt")
        write(palsVect, file = namef)
        names(palsL) <- colnames(x)
        assign(paste(".__PSW_PALS.", namef, sep = ""),
               palsL, env = .GlobalEnv)
    }
  out$array.names <- colnames(x)
  return(out)
}


pSegmentWavelets <- function(x, chrom.numeric, mergeSegs = TRUE,
                             minDiff = 0.25,
                             minMergeDiff = 0.05,
                             thrLvl = 3, initClusterLevels = 10, ...) {
  stop.na.inf(chrom.numeric)
  stop.na.inf(x)
  warn.too.few.in.chrom(chrom.numeric)
  
##     ncloneschrom <- tapply(x[, 1], chrom.numeric, function(x) length(x))
##     if((thrLvl == 3) & ((max(ncloneschrom) > 1096) | (min(ncloneschrom) < 21)))
##         warningsForUsers <-
##             c(warningsForUsers,
##               paste("The number of clones/genes is either",
##                     "larger than 1096 or smaller than 21",
##                     "in at least one chromosome. The wavelet",
##                     "thresholding of 3 might not be appropriate."))
    Nsamps  <- ncol(x)
    uniq.chrom <- unique(chrom.numeric)

    datalist <- list()
    klist <- 1
    for(i in 1:Nsamps) {
        for (j in uniq.chrom) {
            datalist[[klist]] <- x[chrom.numeric == j, i]
            klist <- klist + 1
        }
    }
    thismdiff <- if(mergeSegs) minMergeDiff else minDiff
##    force(thismdiff)
    funwv <- function(ratio){ # , thrLvl, minDiff, initClusterLevels) {
        wc   <- modwt(ratio, "haar", n.levels=slave_thrLvl)
        thH  <- our.hybrid(wc, max.level=slave_thrLvl, hard=FALSE)
        recH <- imodwt(thH)
        ## cluster levels
        pred.ij <- segmentW(ratio, recH, minDiff=slave_minDiff,
                            n.levels = slave_initClusterLevels)
        labs <- as.character(1:length(unique(pred.ij)))
        state <- as.integer(factor(pred.ij, labels=labs))
        return(cbind(Observed = ratio,
                     Smoothed = pred.ij,
                     State = state))
    }
    out0 <- papply(datalist, funwv,
                    list(slave_thrLvl = thrLvl,
                         slave_minDiff = thismdiff,
                         slave_initClusterLevels = initClusterLevels))

    ## list with one entry per array
    out <- list()
    klist <- 1
    for(i in 1:Nsamps) {
        out[[i]] <- out0[[klist]]
        for(j in uniq.chrom[-1]) {
            klist <- klist + 1
            out[[i]] <- rbind(out[[i]], out0[[klist]])
        }
    }

    ## if merging, call ourMerge
    if(!mergeSegs) {
        outl <- list()
        outl$segm <- out
        outl$chrom.numeric <- chrom.numeric
        outl <- add.names.as.attr(outl, colnames(x))
        class(outl) <- c(class(out), "CGH.wave", "adacgh.generic.out")
        return(outl)
    } else {
        datalist <- list()
        for (i in 1:ncol(x)) {
            datalist[[i]] <- list()
            datalist[[i]]$logr <- x[, i]
            datalist[[i]]$pred  <- out[[i]][, 2]
        }
        papout <- papply(datalist,
                         function(z)  ourMerge(z[[1]], z[[2]]))
        outl <- list()
        outl$segm <- papout
        outl$chrom.numeric <- chrom.numeric
        outl <- add.names.as.attr(outl, colnames(x))
        class(outl) <- c(class(out), "CGH.wave", "CGH.wave.merged",
                         "adacgh.generic.out")
        return(outl)
    }
        
}


pSegmentDNAcopy <- function(x, chrom.numeric, parall = "arr", ...) {
  stop.na.inf(chrom.numeric)
  stop.na.inf(x)
  warn.too.few.in.chrom(chrom.numeric)

##  mydcat3()
  
 
    if (parall == "auto")
        parall <- ifelse(ncol(x) > 75, "arr", "axc")
    if (parall == "arr") {
        cat("\n    running arr version \n")
        return(pSegmentDNAcopy_A(x, chrom.numeric, ...))
    }
    if (parall == "axc") {
        cat("\n    running axc version. Beware smoothing is not the same as original! \n")
        return(pSegmentDNAcopy_axc(x, chrom.numeric, ...))
    }
}



pSegmentDNAcopy_axc <- function(x, chrom.numeric, smooth = TRUE,
                            alpha=0.01, nperm=10000, kmax=25, nmin=200,
                            eta = 0.05, overlap=0.25, trim = 0.025,
                            undo.prune=0.05, undo.SD=3, ...) {
    nsample <- ncol(x)
    nchrom <- unique(chrom.numeric)
    datalist <- list()
    klist <- 1
    for(i in 1:nsample) {
        for (j in nchrom) {
            datalist[[klist]] <- x[chrom.numeric == j, i]
            klist <- klist + 1
        }
    }
    
    if (nperm == 10000 & alpha == 0.01 & eta == 0.05) {
        sbdry <- default.DNAcopy.bdry
    } else {
        max.ones <- floor(nperm * alpha) + 1
        sbdry <- getbdry(eta, nperm, max.ones)
    }
    sbn <- length(sbdry)

    papply_common <- list(slave_alpha        = alpha,
                          slave_nperm        = nperm,
                          slave_kmax         = kmax,
                          slave_nmin         = nmin,
                          slave_overlap      = overlap,
                          slave_trim         = trim,
                          slave_undo.prune   = undo.prune,
                          slave_undo.SD      = undo.SD,
                          slave_sbdry        = sbdry,
                          slave_sbn          = sbn)

    
    if(smooth) {
        out0 <- papply(datalist, function(z)
                       wrapperDNAcopySmooth(z,
                                            alpha =         slave_alpha,     
                                            nperm =         slave_nperm,    
                                            kmax =          slave_kmax,      
                                            nmin =          slave_nmin,      
                                            overlap =       slave_overlap,   
                                            trim =          slave_trim,      
                                            undo.prune =    slave_undo.prune,
                                            undo.SD =       slave_undo.SD,   
                                            sbdry =         slave_sbdry,     
                                            sbn =           slave_sbn),
                       papply_common)
    } else {
        out0 <- papply(datalist, function(z)
                       wrapperDNAcopyNoSmooth(z,
                                            alpha =         slave_alpha,     
                                            nperm =         slave_nperm,    
                                            kmax =          slave_kmax,      
                                            nmin =          slave_nmin,      
                                            overlap =       slave_overlap,   
                                            trim =          slave_trim,      
                                            undo.prune =    slave_undo.prune,
                                            undo.SD =       slave_undo.SD,   
                                            sbdry =         slave_sbdry,     
                                            sbn =           slave_sbn),
                       papply_common)
            }
    matout0 <- matrix(unlist(out0), ncol = nsample)
    rm(datalist)
    rm(out0)
    datalist <- list()
    for(i in 1:nsample) {
        datalist[[i]] <- list()
        datalist[[i]]$logr <- x[, i]
        datalist[[i]]$pred <- matout0[, i]
    }
    ### FIXME: eh???? we always merge here!!!!!
    out <- papply(datalist, function(z) ourMerge(z$logr, z$pred))
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    outl <- add.names.as.attr(outl, colnames(x))
    class(outl) <- c("DNAcopy", "adacgh.generic.out")
    return(outl)
}


segmentPlot <- function (x, geneNames, yminmax,
                         idtype = "ug", organism = "Hs",
                         arrays = NULL,
                         chroms = NULL,
                         colors = c("orange", "red", "green", "blue", 
                           "black"), html_js = FALSE, superimp = FALSE,
                         imgheight = 500,
                         genomewide_plot = FALSE,
                         ...) {
  if(length(yminmax) != 2) {
    stop("yminmax must exist, and must be a two-element vector")
  }
  if (is.null(arrays)) {
      arrays <- 1:length(x$segm)
  }
  if (inherits(x, c("adacgh.generic.out"))) {
    if(length(geneNames) != length(x$chrom.numeric)) {
      stop("lenght of geneNames must equal length of x$chrom.numeric")
    }

    geneLoc <- if (inherits(x, "mergedBioHMM")) x$pos else NULL
    if (inherits(x, "CGH.ACE.summary")) {
      original.pos <- 2
      segment.pos <- NULL
    }
    else {
      original.pos <- 1
      segment.pos <- 2
    }
    if (inherits(x, "CGH.wave") & (!inherits(x, "CGH.wave.merged"))) {
      colors <- c(rep(colors[1], 3), colors[4], colors[5])
    }
    
    if (superimp) {
      ## Does not use arrays or chroms parameters
      ## The superimp option will soon be deprecated
      ## I move it here to feel free to delete stuff later
      plot.cw.superimpA(x$segm, x$chrom.numeric, geneNames = geneNames, 
                        main = "All_arrays", colors = colors, ylim = yminmax, 
                        idtype = idtype, organism = organism, geneLoc = geneLoc, 
                        html_js = html_js, imgheight = imgheight)
      cat("\n gc after plot.cw.superimpose \n")
      print(gc())
      plot.gw.superimp(res = x$segm,
                       chrom = x$chrom.numeric,
                       ylim = yminmax,
                       geneNames = geneNames,
                       imgheight = imgheight,
                       main = "All_arrays", colors = colors,  
                       geneLoc = geneLoc)
      cat("\n gc after plot.gw.superimp \n")
      print(gc())
    }

    pappl_common <- list(slave_cnum = x$chrom.numeric, 
                         slave_geneNames = geneNames, 
                         slave_idtype = idtype,  
                         slave_organism = organism, 
                         slave_colors = colors,
                         slave_geneLoc = geneLoc, 
                         slave_yminmax = yminmax, 
                         slave_html_js = html_js, 
                         slave_imgheight = imgheight, 
                         slave_genomewide_plot = genomewide_plot,
                         slave_chroms = chroms)

    tmp_papout <- papply(x$segm[arrays], function(z) {
                         plot.adacgh.nonsuperimpose(res = z,
                                                    main = attributes(z)$ArrayName,
                                                    chrom = slave_cnum, 
                                                    colors = slave_colors,
                                                    ylim = slave_yminmax, 
                                                    geneNames = slave_geneNames,
                                                    idtype = slave_idtype,
                                                    organism = slave_organism, 
                                                    geneLoc = slave_geneLoc,
                                                    html_js = slave_html_js,
                                                    imgheight = slave_imgheight,
                                                    genomewide_plot = slave_genomewide_plot,
                                                    chromsplot = slave_chroms)},
                          papply_commondata = pappl_common)
    
    cat("\n gc after plot.adacgh.nonsuperimpose \n")
    print(gc())
  } else if (inherits(x, "CGH.PSW")) {
    if (x$plotData[[1]]$sign < 0) {
      mainsl <- "Losses."
    }  else {
      mainsl <- "Gains."
    }
    l1 <- list()
    for (i in 1:length(x$plotData)) {
      l1[[i]] <- list()
      l1[[i]]$lratio <- x$plotData[[i]]$logratio
      l1[[i]]$sign <- x$plotData[[i]]$sign
      l1[[i]]$swt.perm <- x$plotData[[i]]$swt.perm
      l1[[i]]$rob <- x$plotData[[i]]$rob
      l1[[i]]$swt.run <- x$plotData[[i]]$swt.run
      l1[[i]]$p.crit <- x$plotData[[i]]$p.crit.bonferroni
      l1[[i]]$chrom <- x$plotData[[i]]$chrom
      l1[[i]]$arrayname <- x$array.names[i]
    }
    tmp_papout <- papply(l1,
                         function(z) {
                           cat("\n Doing sample ", z$arrayname, "\n")
                           sw.plot3(logratio = z$lratio,
                                    sign = z$sign,
                                    swt.perm = z$swt.perm, 
                                    rob = z$rob,
                                    swt.run = z$swt.run,
                                    p.crit = z$p.crit, 
                                    chrom = z$chrom,
                                    main = paste(slave_main, z$arrayname, 
                                      sep = ""),
                                    geneNames = slave_geneNames,
                                    idtype = slave_idtype, 
                                    organism = slave_organism,
                                    html_js = slave_html_js,
                                    imgheight = slave_imgheight)},
                          papply_commondata = list(slave_main= mainsl,  
                            slave_geneNames = geneNames, slave_idtype = idtype,
                            slave_organism = organism, 
                            slave_html_js = html_js,
                            slave_imgheight = imgheight))
    
  }
  else {
    stop("No plotting for this class of objects")
  }
}

SegmentPlotWrite <- function(data, chrom,
                             mergeSegs, Pos,
                             idtype, organism,
                             method,
                             geneNames,
                             commondata,
                             colors = c("orange", "red", "green", "blue", "black"),
                             html_js = FALSE,
                             superimp = FALSE,
                             genomewide_plot = FALSE,
                             imgheight = 500,
                             ...) {
    ymax <- max(data)
    ymin <- min(data)
    
    fseg <- get(paste("pSegment", method, sep = ""))
    trySegment <- try(
                   segmres <- fseg(data, chrom,
                                   Pos = Pos,
                                   mergeSegs = mergeSegs, ...)
                   )
    if(inherits(trySegment, "try-error"))
        caughtOurError(trySegment)
    cat("\n\n Segmentation done \n\n")

    save.image()
    cat("\n\n Saving segmentation results as segmres.RData \n\n")
    save(segmres, file = "segmres.RData")

    tryMCR <- try(doMCR(segmres$segm, chrom.numeric = chrom,
                         data = data,
                         Pos = Pos, ...))
    if(inherits(tryMCR, "try-error"))
        caughtOurError(tryMCR)
    if(inherits(segmres, "DNAcopy") & (mergeSegs == FALSE)) {
      class(segmres) <- c(class(segmres), "adacgh.generic.out")
      warning("Forcing plotting of DNAcopy object with merge = FALSE.",
              " But this might not be what you want.")
    }
    yminmax <- c(min(as.matrix(data)),
                 max(as.matrix(data)))
    tryPlot <- try(segmentPlot(segmres,
                               geneNames = geneNames,
                               yminmax = yminmax,
                               idtype = idtype,
                               organism = organism,
                               colors = colors,
                               html_js = html_js,
                               superimp = superimp,                               
                               imgheight = imgheight,
                               genomewide_plot = genomewide_plot))
    if(inherits(tryPlot, "try-error"))
        caughtOurError(tryPlot)
    cat("\n\n Plotting done \n\n")

    tryWrite <- try(writeResults(segmres,
                                data, commondata))
    if(inherits(tryWrite, "try-error"))
        caughtOurError(tryWrite)
}                                









#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
###
###              Analysis
###                 
###
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################





#######################################################
#######################################################
#######################################################
###
###              HMM and BioHMM
###                 
###
#######################################################
#######################################################
#######################################################

hmmWrapper0 <- function(logratio, Chrom, Pos = NULL) {
    ## Fit HMM, and do mergeLevels
    Clone <- 1:length(logratio)
    if(is.null(Pos)) Pos <- Clone
    obj.aCGH <- create.aCGH(data.frame(logratio),
                            data.frame(Clone = Clone,
                                       Chrom = Chrom,
                                       kb = Pos))
    res <- find.hmm.states(obj.aCGH, aic = TRUE, bic = FALSE)
    hmm(obj.aCGH) <- res
    out <- ourMerge(obj.aCGH$hmm$states.hmm[[1]][, 8],
                      obj.aCGH$hmm$states.hmm[[1]][, 6]) 
    return(out)
}



hmmWrapper <- function(logratio) {
    ## Fit HMM, and return the predicted
    ## we do not pass Chrom since we only fit by Chrom.
    Pos <- Clone <- 1:length(logratio)
    Chrom <- rep(1, length(logratio))
    obj.aCGH <- create.aCGH(data.frame(logratio),
                            data.frame(Clone = Clone,
                                       Chrom = Chrom,
                                       kb = Pos))
    res <- find.hmm.states(obj.aCGH, aic = TRUE, bic = FALSE)
    hmm(obj.aCGH) <- res
    return(obj.aCGH$hmm$states.hmm[[1]][, 6])
}


BioHMMWrapper <- function(logratio, Pos) {
    cat("\n       .... running BioHMMWrapper \n")
    ydat <- matrix(logratio, ncol=1)
    n <- length(ydat)
    res <- try(myfit.model(sample = 1, chrom = 1, dat = matrix(ydat, ncol = 1),
                         datainfo = data.frame(Name = 1:n, Chrom = rep(1, n),
                         Position = Pos)))
    if(inherits(res, "try-error")) {
        return(res)
    } else {
        return(res$out.list$mean)
    }
}

#######################################################
#######################################################
#######################################################
###
###                 Merging
###                 
###
#######################################################
#######################################################
#######################################################



ourMerge <- function(observed, predicted,
                       merge.pv.thresh = 1e-04,
                       merge.ansari.sign = 0.05,
                       merge.thresMin = 0.05,
                       merge.thresMax = 0.5) {

    cat("\n        Starting merge \n")
    segmentus2 <-
        mergeLevelsB(vecObs  = observed,
                    vecPred = predicted,
                    pv.thres = merge.pv.thresh,
                    ansari.sign = merge.ansari.sign,
                    thresMin = merge.thresMin,
                    thresMax = merge.thresMax)$vecMerged

    classes.ref <- which.min(abs(unique(segmentus2)))
    classes.ref <- unique(segmentus2)[classes.ref]
    ref <- rep(0, length(segmentus2))
    ref[segmentus2 > classes.ref] <- 1
    ref[segmentus2 < classes.ref] <- -1
    cat("\n        Done  merge \n")
    return(cbind(Observed = observed,
                 MergedMean = segmentus2,
                 Alteration = ref))
}



#######################################################
#######################################################
#######################################################
###
###                 GLAD
###                 Hupe et al.
###
#######################################################
#######################################################
#######################################################



gladWrapper <- function(x, Chrom, Pos = NULL) {
  Pos <- if (is.null(Pos)) (1:length(x)) else Pos
  tmpf <- data.frame(LogRatio = x,
                     PosOrder = Pos,
                     Chromosome = Chrom)
  tmpf <- list(profileValues = tmpf)
  class(tmpf) <- "profileCGH"
  outglad <- glad.profileCGH(tmpf)
  return(cbind(Observed = x, Smoothed = outglad$profileValues$Smoothing,
               State = outglad$profileValues$ZoneGNL))
}




#######################################################
#######################################################
#######################################################
###
###                 CGHseg
###                 Piccard et al. using tilingArray
###
#######################################################
#######################################################
#######################################################


piccardsKO <- function(loglik, n, s) {
    ## return the optimal number of segments, as in
    ## piccard et al., p. 13. k is number of segments, not breakponts.
    dks <- c(NA, diff(loglik, lag = 1, differences = 2))
    dkthresh <- s * n
    okdk <- which(dks < dkthresh)
    if(length(okdk) > 0) {
        return(max(okdk))
    } else {
        return(1)
    }
}


piccardsStretch01 <- function(obj, k, n, logratio) {
    if(k > 1) {
        poss <- obj@breakpoints[[k]]
        start <- c(1, poss)
        end <- c(poss - 1, n)
        smoothedC <- mapply(function(start, end) mean(logratio[start: end]), start, end)
        reps <- diff(c(start, n + 1))
        smoothed <- rep(smoothedC, reps)
        state <- rep(1:k, reps)
    } else { ## only one segment
        smoothed <- rep(mean(logratio), n)
        state <- rep(1, n)
    }
    return(cbind(smoothed = smoothed, state = state))
}

CGHsegWrapper <- function(logratio, chrom.numeric,
                          s, maxseg = NULL,
                          maxk = NULL,
                          verbose = TRUE,
                          domergeLevels = TRUE) {
    ## Using merge levels now if domergeLevels = TRUE
    uchrom <- unique(chrom.numeric)
    segmeans <- NULL
    segstates <- NULL
    for (ic in uchrom) {
        y <- logratio[chrom.numeric == ic]
        n <- length(y)
        obj1 <- tilingArray:::segment(y,
                                      maxseg = ifelse(is.null(maxseg), n/2, maxseg),
                                      maxk = ifelse(is.null(maxk), n, maxk))
        optk <- piccardsKO(obj1@logLik, n, s)
        if (verbose) {
            cat("\n Chromosome ", ic, ";  Optimal k ", optk, "\n")
        }
        finalsegm <- piccardsStretch01(obj1, optk, n, y)
        segmeans <- c(segmeans, finalsegm[, 1])
        segstates <- c(segstates, finalsegm[, 2])
    }
    if(domergeLevels) {
        tmp <- ourMerge(logratio, segmeans)
        segmeans <- tmp[, 2]
        segstates <- tmp[, 3]
    }
    return(cbind(Observed = logratio, SmoothedMean = segmeans,
                 Alteration = segstates))
}



#### Choosing a good s for Piccard's approach
####  and collapsing levels


### Verify with Piccard's paper, figure 1.

## coriel.data <- read.table("gm03563.txt", header = TRUE)
## cd3 <- coriel.data[coriel.data$Chromosome == 3, 3]
    

## out.lai <- CGHsegWrapper(cd3, optK = "L")  ## k = 14
## out.our <- CGHsegWrapper(cd3, optK = "O")  ## k = 2

## so our implementation seems correct. look at where the
## breakpoint is located, etc, and it is like figure 1 of
## Picard's paper.





#######################################################
#######################################################
#######################################################
###
###                 DNA copy
###                 Olshen and Venkatraman
###
#######################################################
#######################################################
#######################################################

mpiCBS <- function(wdir = getwd()) {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(cghMCR))
    mpi.remote.exec(library(cluster))
    mpi.remote.exec(library(DNAcopy))
    mpi.remote.exec(library(ADaCGH))
    mpi.remote.exec(library(aCGH))
    mpi.bcast.Robj2slave(wdir)
    mpi.remote.exec(setwd(wdir))
}


internalSmoothCNA <- function(genomdat,
                              smooth.region = 2, outlier.SD.scale = 4,
                              smooth.SD.scale = 2, trim = 0.025) {
    ## this is just the original smoothCNA funct. adapted to use
    ## a single array *chromosome and to be parallelized and fed to internalDNAcopy
   cat("\n      Starting smoothing \n")
    ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
   trimmed.SD <- sqrt(adacgh_trimmed.variance(genomdat[ina], trim))
   outlier.SD <- outlier.SD.scale * trimmed.SD
   smooth.SD <- smooth.SD.scale * trimmed.SD
   k <- smooth.region
   n <- length(genomdat[ina])
   
   smoothed.data <-
       sapply(1:n,
              function(i, x, n, nbhd, oSD, sSD) {
                  xi <- x[i]
                  nbhd <- i + nbhd
                  xnbhd <- x[nbhd[nbhd > 0 & nbhd <= n]]
                  if (xi > max(xnbhd) + oSD) 
                      xi <- median(c(xi, xnbhd)) + sSD
                  if (xi < min(xnbhd) - oSD) 
                      xi <- median(c(xi, xnbhd)) - sSD
                      xi
              },
              genomdat[ina], n, c(-k:-1, 1:k), outlier.SD, smooth.SD)
   genomdat[ina] <- smoothed.data
   cat("\n       Done smoothing \n")
   genomdat
}
    

    
internalDNAcopy <- function(acghdata,
                            alpha,
                            nperm,
                            kmax,
                            nmin,
                            overlap, 
                            trim,
                            undo.prune,
                            undo.SD,                            
                            sbdry,
                            sbn) {
    ## tries to follow the original "segment"
    cat("\n        Starting segmentation \n")

    data.type <- "logratio"
    p.method <- "hybrid"
    window.size <- NULL
    undo.splits <- "none"
    genomdati <- acghdata
    min.width <- 2
    ina <- which(!is.na(genomdati) & !(abs(genomdati)==Inf))

    ## The code allows for dealing with NA and Inf, but would need to
    ## adjust other functions (as different arrays would have different
    ## length of pos, genenames, etc. So for now stop:
    if (length(ina) != length(genomdati))
        stop("Either an NA or an infinite in the data")

    genomdati <- genomdati[ina]
    trimmed.SD <- sqrt(adacgh_trimmed.variance(genomdati, trim))
    segci <- adacgh_changepoints(genomdati, data.type = "logratio",
                          alpha = alpha, sbdry = sbdry, sbn = sbn,
                          nperm = nperm, p.method = p.method,
                          ##                               window.size = window.size, 
                          ##                               overlap = overlap,
                          kmax = kmax, nmin = nmin,
                          trimmed.SD = trimmed.SD,
                          undo.splits = undo.splits,
                          undo.prune = undo.prune,
                          undo.SD = undo.SD, verbose = 2,
                          min.width = min.width)

    sample.lsegs <- segci$lseg
    sample.segmeans <- segci$segmeans
    if(length(sample.lsegs) != length(sample.segmeans))
        stop("Something terribly wrong: length(sample.lsegs) != length(sample.segmeans).")
    stretched.segmeans <- rep(sample.segmeans, sample.lsegs)
    cat("\n        Done segmentation \n")

    return(stretched.segmeans)
}



    

wrapperDNAcopySmooth <- function(data, alpha, nperm, kmax, nmin, overlap, trim,
                                 undo.prune, undo.SD, sbdry, sbn)  {
    smooth.region <- 2
    outlier.SD.scale <- 4
    smooth.SD.scale <- 2
    data2 <- internalSmoothCNA(data, smooth.region, outlier.SD.scale,
                              smooth.SD.scale, trim)
    outseg <- internalDNAcopy(data2, alpha, nperm,  kmax,  nmin, overlap,   
                              trim, undo.prune, undo.SD, sbdry, sbn)
    outseg
}

wrapperDNAcopyNoSmooth <- function(data, alpha, nperm, kmax, nmin, overlap, trim,
                                 undo.prune, undo.SD, sbdry, sbn) {
    outseg <- internalDNAcopy(data, alpha, nperm,  kmax,  nmin, overlap,   
                              trim, undo.prune, undo.SD, sbdry, sbn)
    outseg
}






#######################################################
#######################################################
#######################################################
###
###            Wavelet approach
###            Hsu et al.
###
#######################################################
#######################################################
#######################################################

mpiWave <- function(wdir = getwd()) {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(cluster))
    mpi.remote.exec(library(waveslim))
    mpi.remote.exec(library(ADaCGH))
    mpi.bcast.Robj2slave(wdir)
    mpi.remote.exec(setwd(wdir))    
}



#### The first part is the code as provided by Hsu and Grove.
####  Below are my (R.D.-U.) modifications

####################################################################
## Thresholding functions (slightly diff. than ones in 'waveslim'
####################################################################


our.sure <- function (wc, max.level = 4, hard = TRUE)
{
    if (hard) { shrink <- function(w,s) w*(abs(w)>s)
    } else shrink <- function(w,s) sign(w)*(abs(w)-s)*(abs(w)>s)

    for (i in 1:max.level) {
        wci <- wc[[i]]
        ni <- length(wci)
        factor <- mad(wci)

        sxi <- sort(abs( wci/factor ))^2
        s <- cumsum(sxi) + ((ni - 1):0) * sxi
        risk <- (ni - (2 * (1:ni)) + s)/ni

        surethresh <- sqrt(sxi[order(risk)[1]])*factor
        wc[[i]] <- shrink(wci,surethresh)
    }
    wc
}


nominal.thresh <- function (wc, max.level = 4, hard = TRUE, sig.lvl=.05)
  ## If you want threshold all but the coefficients significant at 
  ## level .05 (two-sided, and assuming normality) then set sig.lvl=.05
{
     if (hard) { shrink <- function(w,t) w*(abs(w)>t)
     } else shrink <- function(w,t) sign(w)*(abs(w)-t)*(abs(w)>t)

     for (i in 1:max.level) {
         wci <- wc[[i]]
         ni <- length(wci)
         factor <- mad(wci)

         sd.thresh <- qnorm(1-sig.lvl/2)*factor
         wc[[i]] <- shrink(wci,sd.thresh)
     }
     wc
}


our.hybrid <- function (wc, max.level = 4, hard=TRUE)
{
     if (hard) { shrink <- function(w,t) w*(abs(w)>t)
     } else shrink <- function(w,t) sign(w)*(abs(w)-t)*(abs(w)>t)

     for (i in 1:max.level) {
         wci <- wc[[i]]
         ni <- length(wci)
         factor <- mad(wci)

         if ((sum((wci/factor)^2)-ni)/ni <= sqrt(log2(ni)^3/ni)) {
              ## If not enough spread in coefficients, use
              ## 'universal threshold'
              unithresh <- factor*sqrt(2*log(ni))
              wc[[i]] <- shrink(wci,unithresh)
         }
         else
         {
             ## otherwise use sure threshold
             sxi <- sort(abs( wci/factor ))^2
             s <- cumsum(sxi) + ((ni - 1):0) * sxi
             risk <- (ni - (2 * (1:ni)) + s)/ni

             surethresh <- sqrt(sxi[order(risk)[1]])*factor

             wc[[i]] <- shrink(wci,surethresh)
         }
     }
     wc
}
###################################################################
##               End of Thresholding functions                   ##
###################################################################




####################################################################
##  segment() is a function that:
##  (a) clusters the threshheld data, 
##  (b) collapses together clusters "closer" than 'minDiff',
##  (c) return the median value of the cluster to which each data 
##      point was assigned as its predicted value 
## 
####################################################################

## I (RDU) rename segment to segmentW to prevent confussion

segmentW <- function(obs.dat, rec.dat, minDiff=0.25, n.levels=10) {
    ## 'obs.dat' is OBServed DATa
    ## 'rec.dat' is "REConstructed DATa" following wavelet thresholding
    ## 'n.levels' is the initial number of clusters to form 
    ## 'minDiff' is the MINimum (absolute) DIFFerence between the medians
    ##           of two adjacent clusters for them to be considered truly
    ##           different.  Clusters "closer" together than this are
    ##           collapsed together to form a single cluster.

    pam.out <- pam(rec.dat, n.levels)
    clust.indx <- pam.out$clustering
    med <- as.vector(pam.out$medoids)

    ord  <- order(med)
    ord.med  <- med[ord]

    ## ord.lab contains the unique group labels (i.e. if 4 groups then 1:4)
    ## ordered according to the values in 'med'
    ord.lab <- (1:length(med))[ord]

    diff.med   <- diff(ord.med)

    done <- !(min(diff.med) < minDiff)
    while (!done) {
  
        ## get labels of groups that are closest together
        w <- which.min(diff.med)
        grp.1 <- ord.lab[w]
        grp.2 <- ord.lab[w+1]

        ## rename grp.2 to grp.1
        clust.indx[which(clust.indx == grp.2)] <- grp.1
        
        ## delete grp.2 from the vector ord.lab
        ord.lab <- ord.lab[which(ord.lab != grp.2)] 

        ## get new medians
        ord.med <- unlist(lapply(ord.lab, function(x,I,dat) median(dat[I==x]), 
                                 clust.indx, rec.dat))

        ## figure out if we are done or not
        if (length(ord.med)>1) {
            diff.med <- diff(ord.med)
            done <- !(min(diff.med) < minDiff)
        } else done <- TRUE
    }

    ## create and output vector of same length as original data
    new.medians <- integer(length=length(rec.dat))
    for (i in 1:length(ord.lab)){
        ind <- which(clust.indx == ord.lab[i])
        new.medians[ind] <- median(obs.dat[ind])
    }
    return(new.medians)
}
####################################################################






wave.aCGH <- function(dat, chrom, minDiff) {
## level to use for wavelet decomposition and thresholding
## The 'recommended' level is floor(log2(log(N))+1)), which
## equals 3 for:  21 <= N <= 1096
    thrLvl <- 3

    ncloneschrom <- tapply(dat[, 1], chrom, function(x) length(x))
    if((max(ncloneschrom) > 1096) | (min(ncloneschrom) < 21))
        warningsForUsers <-
            c(warningsForUsers,
              paste("The number of clones/genes is either",
                    "larger than 1096 or smaller than 21",
                    "in at least one chromosome. The wavelet",
                    "thresholding of 3 might not be appropriate."))
    
    Nsamps  <- ncol(dat)
    uniq.chrom <- unique(chrom)

## construct the list:
## the code below gives some partial support for missings.
    ##  but I need to carry that along, and since we are not dealing
    ##  with missings now, I just re-writte ignoring any NA,
    ##  since, by decree, we have no NAs.
##     datalist <- list()
##     klist <- 1
##     for(i in 1:Nsamps) {
##         ratio.i <- dat[,i]
##         noNA  <- !is.na(ratio.i)
##         for (j in uniq.chrom) {
##             chr.j <- (chrom == j)
##             use.ij <- which(noNA & chr.j)
##             datalist[klist] <- ratio.i[use.ij]
##             klist <- klist + 1
##         }
##     }

    datalist <- list()
    klist <- 1
    for(i in 1:Nsamps) {
        ratio.i <- dat[,i]
        for (j in uniq.chrom) {
            chr.j <- (chrom == j)
            use.ij <- which(chr.j)
            datalist[[klist]] <- ratio.i[use.ij]
            klist <- klist + 1
        }
    }
    
    funwv <- function(ratio, thrLvl, minDiff) {
        wc   <- modwt(ratio, "haar", n.levels=thrLvl)
        ## These are the three different thresholding functions used
        ##thH  <- our.sure(wc, max.level=thrLvl, hard=FALSE)
        thH  <- our.hybrid(wc, max.level=thrLvl, hard=FALSE)
        ##thH  <- nominal.thresh(wc, max.level=thrLvl, hard=FALSE, sig=.05)
        ## reconstruct the thresheld ('denoised') data
        recH <- imodwt(thH)
        ## Categorize the denoised data then combine ("merge") levels that
        ## have predicted values with an absolute difference < 'minDiff' 
        pred.ij <- segmentW(ratio, recH, minDiff=minDiff)
        labs <- as.character(1:length(unique(pred.ij)))
        state <- as.integer(factor(pred.ij, labels=labs))
        return(list(pred.ij = pred.ij, state = state))
    }

    papout <- papply(datalist,
                     function(z) funwv(z,
                                       thrLvl = slave_thrLvl,
                                       minDiff = slave_minDiff),
                      list(slave_thrLvl = thrLvl,
                           slave_minDiff = minDiff))
    pred <- matrix(unlist(lapply(papout, function(x) x$pred.ij)),
                   ncol = Nsamps)
    state <- matrix(unlist(lapply(papout, function(x) x$state)),
                   ncol = Nsamps)
                   
    out <- list(Predicted =pred, State = state)
    return(out)
}


    
plateau.wavelets <- function(res, xcenter,
                             by.array = TRUE, superimpose = FALSE,
                             ylim = NULL) {

    if(by.array) {
        for(i in 1:ncol(res$Predicted)) {
            op <- order(res$Predicted[, i])
            maint <- colnames(xcenter)[i]
            pcht <- "."
            if(superimpose) {
                if(i > 1) par(new = TRUE)
                maint <- "All, superimposed"
                pcht <- ""
            }
            plot(xcenter[op, i], ylab = "",
                 main = maint,
                 col = "orange", pch = pcht,
                 ylim = ylim)
         
            lines(res$Predicted[op, i], col = "black")
            abline(h = 0, lty = 2, col = "blue")
        }
    } else {
        stretched.data <- as.vector(as.matrix(xcenter))
        stretched.smooth <- as.vector(as.matrix(res$Predicted))
        op <- order(stretched.smooth)
        plot(stretched.data[op], ylab = "",
             main = "All arrays",
             col = "orange", pch = ".",
             ylim = ylim)
        lines(stretched.smooth[op], col = "black")
        abline(h = 0, lty = 2, col = "blue")
    }
}
             
        
#######################################################
#######################################################
#######################################################
###
###            Price-Smith-Waterman
###            
###
#######################################################
#######################################################
#######################################################

library(cgh)


mpiPSW <- function(wdir = getwd()) {
##    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(cluster))
    mpi.remote.exec(library(waveslim))
    mpi.remote.exec(library(cghMCR))
    mpi.remote.exec(library(DNAcopy))
    mpi.remote.exec(library(cgh))
    mpi.remote.exec(library(ADaCGH))
    mpi.bcast.Robj2slave(wdir)
    mpi.remote.exec(setwd(wdir))    

}





my.sw3b <- function(logratio, chrom, sign = -1, p.crit = 0.01,
                    nIter = 1000,
                    prec = 100,
                    name,
                    highest = FALSE, ## identifying highest scoring island can
                    ## cross chromosome boundaries
                    ...) {

    ## like my.sw3, but do not parallelize over chroms. (Will parall. over arrays)

    if(!is.numeric(chrom)) stop("Chrom not numeric; this will cause trouble")
    swt <- sw.threshold(logratio, sign = sign)
    
    swt.run.l <- swt.perm.l <- swt.rob.l <- list()
    chrom.nums <- unique(chrom)

    
    swtlist <- list()
    klist <- 1
    for(i in 1:length(chrom.nums)) {
        swtlist[[klist]] <- swt[chrom == i]
        klist <- klist + 1
    }

    funsw <- function(x) {
        swt.run <- sw(x, trace = FALSE)
        
        swt.perm <- try(sw.perm.test(x, max.nIslands = NULL,
                                        nIter = nIter))
        if(inherits(swt.perm, "try-error")) {
            return(list(swt.run = NA, swt.perm = "swt.perm.try-error",
                        swt.rob = NA))
        } else {
            swt.rob <- sw.rob(x, prec = prec)
            return(list(swt.run = swt.run, swt.perm = swt.perm,
                        swt.rob = swt.rob))
        }
    }

    papout <- lapply(swtlist, funsw) ## funny that this works!

    swt.perm <- unlist(lapply(papout, function(x) x$swt.perm))
    if(any(swt.perm == "swt.perm.try-error")) {
        return("swt.perm.try-error")
    } else {

        swt.rob <- unlist(lapply(papout, function(x) x$swt.rob))
        swt.run.l <- lapply(papout, function(x) x$swt.run)
        
        swt.run <- list()
        swt.run$length <- unlist(lapply(swt.run.l, function(x) x$length))
        swt.run$start <- unlist(lapply(swt.run.l, function(x) x$start))
        
        
        ## recall that now all genes are numbered stgarting at 1 for each chromos
        nsp <- lapply(swt.run.l, function(x) length(x$length))
        npc <- table(chrom)
        npca <- cumsum(npc[-length(npc)])
        to.add <- rep(c(0, npca), nsp)
        
        perm.p.values <- rep(NA, length(logratio))
        
        x0 <- swt.run$start + to.add
        x1 <- x0 + swt.run$length - 1
        
        for(jj in 1:length(swt.perm)) {
            perm.p.values[x0[jj]:x1[jj]] <- swt.perm[jj]
        }
        
        plotdat <- list(logratio = logratio,
                        sign = sign,
                        rob = swt.rob,
                        swt.run = swt.run,
                        swt.perm = swt.perm,
                        p.crit = p.crit,
                        chrom = chrom)
        
        out.values <- cbind(logratio, rep(sign, length(logratio)),
                            swt.rob, perm.p.values)
        colnames(out.values) <-
            paste(name, c(".Original", ".Sign", ".Robust", ".p.value"),
                  sep = "")
        
        out <- list(out=out.values,
                    plotdat = plotdat)
        class(out) <- c(class(out), "CGH.PSW")
        return(out)
    }     
}


    
sw.plot3 <- function (logratio, location = seq(length(logratio)),
                      threshold.func = function(x) median(x) + 
                      0.2 * mad(x), sign = -1, highest = TRUE, expected = NULL, 
                      rob = NULL, legend = TRUE, xlab = "Chromosomal location", 
                      ylab = "log ratio",
                      swt.run,
                      swt.perm,
                      p.crit,
                      chrom, main,
                      geneNames,
                      html = TRUE,
                      nameIm = NULL,
                      idtype = idtype,
                      organism = organism,
                      html_js = html_js,
                      imgheight = imgheight) {   
    ## this puts the chr call
    ## geneNames often = positions.merge1$name

    if(is.null(nameIm)) nameIm <- main
    if(html) {
##        imheight <- imgheight
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imgheight,
                         width = imwidth, ps = 12)
    }

    sw.plot2(logratio, sign = sign, rob = rob, main = main,
             highest = FALSE)

    sign.segments <- which(swt.perm < p.crit)

    if(sign == 1) {
        red.pos <- quantile(logratio[logratio > 0], p = 0.66)
    } else if(sign == -1) {
        red.pos <- quantile(logratio[logratio < 0], p = 0.33)
    }   
    
    if(length(sign.segments)) {
        f1A <- function(index) {
            si <- sign.segments[index]
            segments(swt.run$start[si] - 0.5, red.pos,
                     swt.run$start[si] - 0.5 + swt.run$length[si],
                     red.pos, col = "red", lwd = 2)
        }
        tmp <- mapply(f1A, 1:length(sign.segments))
    }

   
    ## Limits between chromosomes
    LimitChr <- tapply(1:length(logratio), chrom, max)
    abline(v=LimitChr, col="grey", lty=2)
    
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    chrom.nums <- unique(chrom)
    axis(1, at = pos.labels, labels = chrom.nums)

    if(html) {
        lxs <- c(1, LimitChr)
        maxlr <- max(logratio)
        minlr <- min(logratio)
        nd <- 1:length(LimitChr)
        xleft <- lxs[nd]
        names(xleft) <- 1:length(xleft)
        xright <- lxs[nd + 1]

        f1B <- function(xleft, xright, nd)
            imRect(xleft, maxlr, xright, minlr - 10,
                   title = paste("Chromosome", nd),
                   alt = paste("Chromosome", nd),
                   href= paste("Chr", nd, "@", nameIm, ".html", sep =""))
        rectslist <- mapply(f1B, xleft, xright, nd, SIMPLIFY=FALSE)
        for(ll in 1:length(rectslist))
            addRegion(im1) <- rectslist[[ll]]
        createIM2(im1, file = paste(nameIm, ".html", sep = ""))
        imClose3(im1)
    }

    if(html) { ## here is chromosome specific code
        pixels.point <- 3
##        imgheight <- imgheight
        for(cnum in 1:length(chrom.nums)) {
##            cat("\n .... doing chromosome ", cnum, ": ")
            indexchr <- which(chrom == chrom.nums[cnum])
            chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
            chrwidth <- max(chrwidth, 800)
            im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                             height = imgheight, width = chrwidth,
                             ps = 12)

            ## The following seems needed (also inside sw.plot2) for the coords.
            ## of points to work OK
            par(xaxs = "i")
            par(mar = c(5, 5, 5, 5))
            par(oma = c(0, 0, 0, 0))
    
            sw.plot2(logratio[indexchr], location = location[indexchr],
              sign = sign, rob = rob[indexchr],
              main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
              highest = FALSE, detail = TRUE)
           
            if(length(sign.segments)) {
                for(ii in 1:length(sign.segments)) {
                    if((swt.run$start[sign.segments[ii]] >= location[indexchr[1]])
                       & (swt.run$start[sign.segments[ii]] <= location[indexchr[length(indexchr)]])) {
                           x0 <- swt.run$start[sign.segments[ii]] - 0.5
                           x1 <- x0 + swt.run$length[sign.segments[ii]]
                           y0 <- red.pos
                           y1 <- red.pos
                           segments(x0, y0, x1, y1, col = "red", lwd = 2)
                       }
                }
            } ## end of segments part

            ## The within chromosome map for gene names
            usr2pngCircle <- function(x, y, rr = 2) {
                xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
                r <- abs(xyrc[2, 1] - xyrc[3, 1])
                return(c(xyrc[1, 1], xyrc[1, 2], r))
            }
            ccircle <- mapply(usr2pngCircle, location[indexchr],
                                logratio[indexchr])
            nameChrIm <- paste("Chr", chrom.nums[cnum], "@", nameIm, sep ="")
            write(ccircle, file = paste("pngCoordChr_", nameChrIm, sep = ""),
                  sep ="\t", ncolumns = 3)
            write(as.character(geneNames[indexchr]),
                  file = paste("geneNamesChr_", nameChrIm, sep = ""))
            imClose3(im2)
            ## call the Python function
            if(html_js)
                system(paste(.python.toMap.py,  nameChrIm, 
                             idtype, organism, sep = " "))
            

        } ## looping over chromosomes
    } ## if html
    cat("\n")
}

          





sw.plot2 <- function (logratio, location = seq(length(logratio)),
                      threshold.func = function(x) median(x) +
                      0.2 * mad(x), sign = -1, highest = TRUE, expected = NULL,
                      rob = NULL, legend = TRUE, xlab = "Chromosomal location",
                      ylab = "log ratio", detail = FALSE, ...)
{
  my.line <- function(x, y, ...) {
    len <- length(x)
    run <- rle(y)[[1]]
    run.len <- length(run)
    j <- 1
    m <- 2 * x[1] - x[2]
    if (run.len == 1)
      lines(x = c(3/2 * x[1] - 1/2 * x[2], 3/2 * x[len] -
              1/2 * x[len - 1]), y = c(y[1], y[len]), ...)
    else {
      for (i in 1:(run.len - 1)) {
        k <- run[i]
        lines(x = c((x[j] + m)/2, (x[j + k - 1] + x[j +
                k])/2), y = c(y[j], y[j]), ...)
        lines(x = rep((x[j + k - 1] + x[j + k])/2, 2),
              y = c(y[j], y[j + k]), ...)
        m <- x[j + k - 1]
        j <- j + k
      }
      lines(x = c((m + x[j])/2, 3/2 * x[len] - 1/2 * x[len -
              1]), y = c(y[j], y[j]), ...)
    }
  }
  island.line <- function(x, y, start = 1, len = length(x),
                          edge = 0, ...) {
    if (is.null(start) || is.null(len) || length(start) ==
        0 || length(len) == 0 || len <= 0)
      return
    lenx <- length(x)
                                x1 <- c(2 * x[1] - x[2], x, 2 * x[lenx] - x[lenx - 1])
    x2 <- x1[1:(lenx + 1)] + diff(x1)/2
    lines(x = c(rep(x2[start], 2), rep(x2[start + len], 2)),
          y = c(y - edge, y, y, y - edge), ...)
  }
  log2 <- log(2)
  len <- length(logratio)
  par <- par()
  
  if(detail) {
    par(xaxs = "i")
    par(mar = c(5, 5, 5, 5), yaxs = "r")
    par(oma = c(0, 0, 0, 0))
  }
  
  threshold <- threshold.func(sign * logratio)
  plot(y = logratio, x = location, type = "n", xlab = "Chromosomal location", ylab = "",
       ..., axes = FALSE)
  
  rug(location, ticksize = 0.01)
  maxlr <- max(logratio)
  minlr <- min(logratio)
  axis(side = 4, at = seq(minlr, maxlr, length = 5), labels = c("0",
                                                       ".25", " .50", ".75", "1"))
  box()
  axis(2)
  ##    mtext(xlab, side = 1, line = 3, cex = 0.8)
  mtext(ylab, side = 2, line = 3, cex = 0.8)
  ##    abline(h = maxlr, lty = 2)
  ##    abline(h = minlr, lty = 2)
  if (!is.null(rob)) {
    mtext("Robustness", side = 4, line = 3, cex = 0.8)
    my.line(y = (maxlr - minlr) * rob + minlr, x = location,
            col = "#99ffff", lwd = 2)
  }
  if (highest) {
    x <- sign * logratio - threshold
    swx <- sw(x)
    if (length(swx$score)) {
      island.line(x = location, y = maxlr + (maxlr - minlr) *
                  0.02, start = swx$start[1], len = swx$length[1],
                  edge = (maxlr - minlr) * 0.01, col = "green")
    }
  }
  ##    my.line(y = rep(sign * threshold, len), x = location, col = "#00ff00")
  if (!is.null(expected)) {
    my.line(y = log2(expected) - 1, x = location, col = "#ff0000")
  }
  ##    points(y = logratio, x = location, pch = 20, col = "orange", cex = 0.5)
  if(detail) cexp <- 1
  else cexp <- 0.5
  points(y = logratio, x = location, pch = 20, col = "orange", cex = cexp)
  if (legend) { ## I make changes here; essentially, I break a few things!
    ##         legend.str <- c("highest-scoring island", "robustness", "signif.")
    ##         legend.col <- c("green", "#99ffff", "red")
    ##         legend(location[1], minlr + (maxlr - minlr) * 0.25, legend.str,
    ##             lty = rep(1, 3), col = legend.col, cex = 0.8)
    legend.str <- c("robustness", "signif.")
    legend.col <- c("#99ffff", "red")
    legend(location[1], minlr + (maxlr - minlr) * 0.25, legend.str,
           lty = rep(1, 2), col = legend.col, cex = 0.8)
    
  }
  abline(h = 0, lty = 2, col = "blue")
  par(mar = par$mar, yaxs = par$yaxs)
}







#######################################################
#######################################################
#######################################################
###
###            ACE
###            
###
#######################################################
#######################################################
#######################################################



mpiACE <- function(wdir = getwd()) {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(ADaCGH))
    mpi.bcast.Robj2slave(wdir)
    mpi.remote.exec(setwd(wdir))    
}



##data("file.aux.RData")


ace.analysis.C <- function(x, coefs=file.aux, Sdev, array.names="x")
{
  ##Coefficients of file
####  coefs <- read.table(file=file.aux, header=TRUE, sep="\t")

  if(is.null(array.names)) array.names <- "x"
  
  ##Obtain the null distribution of the (L,H)-pairs
  alpha1 <- coefs$alpha1
  alpha2 <- coefs$alpha2
  beta1 <- coefs$beta1
  beta2 <- coefs$beta2
  ACEPgene <- coefs$Pgene
  Ngenes <- length(x)
  Nlevels <- nrow(coefs)
  
  obj.ace <- .C("aceAnalysis", x=as.double(x), sdev=as.double(Sdev),
                Ngenes = as.integer(Ngenes), Nlevels=as.integer(Nlevels),
                alpha1=as.double(alpha1), alpha2=as.double(alpha2),
                beta1=as.double(beta1), beta2=as.double(beta2),
                ACEPgene = as.double(ACEPgene), FDR=as.double(rep(0, Nlevels)),
                calledGenes=as.integer(rep(0, Nlevels)),
                called = as.integer(rep(-1, Ngenes* Nlevels)),
                first= as.integer(rep(-1, Ngenes)),
                last=as.integer(rep(-1, Ngenes)))
                
  ####Clean up the variables and add 1 to first and last
  ####for the difference in array enumeration in C and R
  obj.ace$called[obj.ace$called==-1] <- NA
  obj.ace$called <- obj.ace$called[!is.na(obj.ace$called)]
  obj.ace$called <- matrix(obj.ace$called, Nlevels)
  obj.ace$first[obj.ace$first==-1] <- NA
  obj.ace$first <- obj.ace$first[!is.na(obj.ace$first)]
  obj.ace$first <- obj.ace$first + 1
  obj.ace$last[obj.ace$last==-1] <- NA
  obj.ace$last <- obj.ace$last[!is.na(obj.ace$last)]
  obj.ace$last <- obj.ace$last + 1
  obj.ace <- list(obj.ace$x, obj.ace$FDR, obj.ace$calledGenes,
                  obj.ace$sdev, obj.ace$called,
                  obj.ace$ACEPgene, obj.ace$first, obj.ace$last)
  names(obj.ace) <- c(array.names, "FDR", "calledGenes", "Sdev",
                      "called", "ACEPgene", "first", "last")
  class(obj.ace) <- "ACE"

  obj.ace

}






ace.analysis <-function(x, coefs = file.aux, Sdev, echo=FALSE, array.names="x") {
    Ngenes = length(x)
    if(is.null(array.names)) array.names <- "x"
    
    ##Step 1		
    ##Segmentation
    ##Compute running mean
    
    yhat <- rep(NA, Ngenes)
    yhat[1:2] = x[1:2]
    for (i in 3:(Ngenes - 2)) {
        yhat[i] = 0;
        for (j in (i - 2):(i + 2)) {
            yhat[i] = yhat[i] + x[j]
        }
        yhat[i] = yhat[i]/5
    }
    yhat[Ngenes - 1] = x[Ngenes - 1]
    yhat[Ngenes] = x[Ngenes]
    
    ##Find change points
    
    ##nknots <- 0
    knots <- rep(NA, Ngenes)
    signo <- (yhat>=0)
    for (i in 3:(Ngenes-2)) {
        if(yhat[i-2]>=0 && yhat[i-1]>=0 && yhat[i+1]>=0 && yhat[i+2]>=0) {
            signo[i] <- TRUE
        }
		
        if(yhat[i-2]<0 && yhat[i-1]<0 && yhat[i+1]<0 && yhat[i+2]<0) {
            signo[i] <- FALSE
        }
    }
    knots <- (1:(Ngenes-1))*(signo[1:(Ngenes-1)]!=signo[2:Ngenes])
    knots <- knots[knots>0]
    knots <- c(0, knots, Ngenes)
    
    
    ##Step 2
    ##Feature extraction
    ##Compute size, height, first, last tuples
    
    Nclusters <- length(knots) - 1
    first <- knots[1:Nclusters] + 1
    last <- knots[2:(Nclusters+1)] 
    size <- last - first + 1
    height <- mapply(function(first, last) mean(x[first:last]), first=first, last=last)
    z1 <- size
    z2 <- abs(height)
    
    
    
    ##Step 3
    ##Obtain the null distribution of the (L,H)-pairs
    
    ##Step 6
    ##Report genes
    ##Adjust start/end positions of clusters
    ##Could be optimized
    for (j in 1:Nclusters) {
        fi <- first[j]
        la <- last[j]
        err.opt <- 999999999
        p.opt <- fi
        q.opt <- la
        r <- floor(max(0, min(16, (la-fi)/2)-1))
        for (p in fi:(fi+r)) {
            for(q in la:(la-r)) {
                ##Hay que vigilar que se cumplan las condiciones de los for
                e1 <- ifelse(fi<p,sum(x[fi:(p-1)]^2),0)
                ##Hay que vigilar que se cumplan las condiciones de los for
                e2 <- ifelse(p<=q,sum((x[p:q]-height[j])^2),0)
                ##Hay que vigilar que se cumplan las condiciones de los for
                e3 <- ifelse(q+1<=la,sum(x[(q+1):la]^2),0)
                err <- (e1 + e2 + e3) /(la-fi)
                ####In case that (la-fi)==0
                err <- ifelse(is.nan(err), 99999999999, err)
                if (err < err.opt) {
                    err.opt <- err
                    p.opt <- p
                    q.opt <- q
                }
            }
        }
        if(echo) {
            cat("\nUpdated:", p.opt, "-->", q.opt)
        }
        first[j] <- p.opt
        last[j] <- q.opt
        size[j] <- q.opt - p.opt + 1
    }
    alpha1 <- coefs$alpha1 * Sdev
    alpha2 <- coefs$alpha2 * Sdev
    beta1 <- coefs$beta1 * Sdev
    beta2 <- coefs$beta2 * Sdev
    ##Step 4
    ##Find significant genes
    ##Determine which clusters are inside/outside
    Nlevels <- nrow(coefs)
    called <- matrix(NA, Nlevels, Nclusters)
    
    v1 <- t(as.matrix(mapply(function(alpha1, beta1) {
        z2-(alpha1+beta1*z1) }, alpha1=alpha1, beta1=beta1)))
    v2 <- t(as.matrix(mapply(function(alpha2, beta2) {
        z2-(alpha2+beta2*z1) }, alpha2=alpha2, beta2=beta2)))
###Check that mapply left the matrix in good shape
    if (nrow(v1)!=Nlevels) {
        v1 <- t(v1)
    }
    if (nrow(v2)!=Nlevels) {
        v2 <- t(v2)
    }
    ##number of rejections
    called <- !(v1<=0 & v2>=0)
    
    ##Step 5
    ##Estimate the positive false discovery rate
    
    
    ACEPgene <- coefs$Pgene
    calledClusters <- apply(called, 1, sum)
    calledGenes <- rep(NA, Nlevels)
    for (i in 1:Nlevels) {
        calledGenes[i] <- sum(last[called[i,]] - 
                              first[called[i,]] + 1)
        
    }
    
    FDR <- Ngenes*(1-ACEPgene)/calledGenes	
    
    ace <- list(x, FDR, calledGenes, Sdev, called, ACEPgene, first, last)
    names(ace) <- c(array.names, "FDR", "calledGenes", "Sdev", "called", "ACEPgene", "first", "last")
    class(ace) <- "ACE"
    ace
}

sd.ACE.analysis<- function(obj.ACE.analysis) {
	######We only calculate desv. if the number of genes >= 10, 'minCount'
	x <- obj.ACE.analysis[[1]]
	if (length(x)>=10) {
		first <- obj.ACE.analysis$first
		last <- obj.ACE.analysis$last
		called <- obj.ACE.analysis$called
		##Estimation of the variance parameter: with all samples and chromosomes, 
		##so we save the number of genes to make weighted mean
		##Indexes have changed from 21 to 41 in this version
		indexes <- unlist(mapply(function(f,l)c(f:l), 
			f=first[!called[41,]], l=last[!called[41,]]))	
		n <- length(indexes)
		if (n>1) {
			Sdev <- sqrt(var(x[indexes])*(n-1)/n)
			}
		else {
			Sdev <- 0
			}
		Sdev
	}
	else {
		Sdev <- 0
		Sdev
		}
}

## ace.analysisP <-function(x, coefs, Sdev, array.names) {
##     ace.analysis.C(x, coefs, Sdev, array.names)
## }


firstACEestimate <- function(x, coefs, Sdev, array.names) {
    obj1 <- ace.analysis.C(x, coefs, Sdev, array.names)
    return(sd.ACE.analysis(obj1))
}


ACE <- function(x, chrom.numeric, coefs = file.aux, Sdev=0.2, echo=FALSE) {
    
####### x is log2.ratio
####### Chrom.numeric ---MUST BE NUMERIC-- 
    if (!is.numeric(chrom.numeric)) {
        stop("Chromosome variable must be numeric")
    }
     
    array.names <- colnames(x)
    nchrom <- length(unique(chrom.numeric))
    
    if(is.null(dim(x)) || (dim(x)[2]==1)) {
### this is only parallelization possible, since only one subject
      genes <- split(x, chrom.numeric)
      first.estimate <- papply(genes,
                               function(z) firstACEestimate(z,
                                                            coefs = slave_coefs,
                                                            Sdev = slave_Sdev,
                                                            array.names = slave_array.names),
                                list(slave_coefs = coefs,
                                     slave_Sdev = Sdev,
                                     slave_array.names = array.names))
      Sdevs <- unlist(first.estimate)
      Sdevs <- mean(Sdevs[Sdevs>0])
      if(is.nan(Sdevs)) Sdevs <- 0
      res <- papply(genes, function(z) ace.analysis.C(x = z,
                                                      coefs = slave_coefs,
                                                      Sdev = slave_Sdev,
                                                      array.names = slave_array.names),
                     list(slave_coefs = coefs,
                          slave_Sdev = Sdev,
                          slave_array.names = array.names))
      class(res) <- "ACE"
    } else {
      res <- list()
      nsample <- ncol(x)
      datalist <- list()
      klist <- 1
      for(i in 1:nsample) {
        for (j in 1:nchrom) {
          datalist[[klist]] <- list()
          datalist[[klist]]$logr <- x[chrom.numeric == j, i]
          datalist[[klist]]$arrayn <- array.names[i]
          klist <- klist + 1
        }
      }
      first.estimate <- papply(datalist,
                               function(z) firstACEestimate(z$logr,
                                                            coefs = slave_coefs,
                                                            Sdev = slave_Sdev,
                                                            array.names = z$arrayn),
                                list(slave_coefs = coefs,
                                     slave_Sdev = Sdev))
      Sdevs <- unlist(first.estimate)
      Sdevs <- mean(Sdevs[Sdevs>0])
      if(is.nan(Sdevs)) Sdevs <- 0
      
      res <- papply(datalist,
                    function(z) ace.analysis.C(z$logr,
                                               coefs = coefs,
                                               Sdev = Sdev,
                                               array.names = z$arrayn),
                                list(slave_coefs = coefs,
                                     slave_Sdev = Sdev))
      resout <- list()
      i <- 1
      k <- 0
      kall <- length(res)
      
      while(k < kall) { #if k == klist it will bomb, which we want
        ## as that will signal an error
        ##             cat("\n i is ", i, "\n")
        resout[[i]] <- list()
        for(j in 1:nchrom) {
          k <- k + 1
          ##                 cat("\n     inner k is ", k, "\n")
          resout[[i]][[j]] <- res[[k]]
        }
        class(resout[[i]]) <- "ACE"
        i <- i + 1
      }
      class(resout) <- "ACE.array"
    }
    resout$chrom.numeric <- chrom.numeric
    resout$array.names <- array.names
    invisible(resout)
  }

get.FDR <- function(object, nchrom) {
	 ### Table to select FDR. Only interesting columns
         ### Number of samples
	 if (class(object)=="ACE") {
	 	nsamples <- 1
		}
	else {
		nsamples <- length(object)
		}
	 ### First column
	 calledGenes <- matrix(unlist(lapply(object, "[", "calledGenes")), ncol=nchrom)
         calledGenes <- apply(calledGenes, 1, sum)
         ### Third column
         ###These values are from file 'coeftable.cgc'
         ACEPgene <- object[[1]]$ACEPgene
         FDR <- mapply(function(A, b) {(1-A) / b}, A=ACEPgene, b=calledGenes)
         ###Total of genes
         tot.genes <- sum(sapply(object, function(x) length(x[[1]])))
         FDR <- FDR*tot.genes
         FDR.table <- cbind(1:length(ACEPgene), calledGenes, round(FDR,4))
         colnames(FDR.table) <- c("Index", "Genes", "FDR")
	 invisible(FDR.table)

}


ace.fdr.html.table <- function(FDR.table, outhtml, index) {
    ## Place this at bottom of page
    ## what gets change is this itself and the figures
    cat("<TABLE border>\n", file = outhtml, append=FALSE)
    cat("<tr><td>Index</td><td>Number of Genes with gains/losses</td><td>FDR</td></tr>\n",
        file = outhtml, append = TRUE)
    for(ni in 1:dim(FDR.table)[1]) {
        if(ni == index) 
            cat(paste("<b><tr><td><b>", FDR.table[ni, 1], "</b></td><td><b>",
                      FDR.table[ni, 2], "</b></td><td><b>", FDR.table[ni, 3],
                      "</b></td></tr></b>\n"), file = outhtml, append = TRUE)
        else
            cat(paste("<tr><td>", FDR.table[ni, 1], "</td><td>",
                      FDR.table[ni, 2], "</td><td>", FDR.table[ni, 3],
                      "</td></tr>\n"), file = outhtml, append = TRUE)
    }
    cat("</TABLE>", file = outhtml, append= TRUE)
}

summary.ACE <- function(object, fdr=NULL, html = TRUE,
                        outhtml ="ace.fdrtable.html", ...) {

    chrom.numeric <- object$chrom.numeric
    object$chrom.numeric <- NULL
    array.names <- object$array.names
    object$array.names <- NULL
 

    
	nchrom <- length(object)
	FDR.table <- get.FDR(object, nchrom)

	if (is.null(fdr)) fdr <- 0.15

        index <- which.min(abs(FDR.table[, 3] - fdr))
	print(FDR.table)
        print(paste("Selected index", index))
        cat(FDR.table[index, 3], file ="aceFDR")

        aceFDR.for.output <- FDR.table[index, 3]
        
        if(html) ace.fdr.html.table(FDR.table, outhtml, index)
		
	#Recover altered genes at the FDR level
	start <- sapply(object, "[", "first")
	end <- sapply(object, "[", "last") 
	called <- sapply(object, "[", "called")
	#Recover observations to get the sign of their averages
	obs <- sapply(object, "[", 1)
	called <- sapply(called, function(x, index) x[index,], index=index)
	called <- lapply(called, function(x) x == 1)
	## Select starting and ending points of called genes
	start <- mapply(function(start, called) { start[called]}, start=start, called=called)
	end <- mapply(function(end, called) { end[called]}, end=end, called=called)
	####Index for every cluster of genes
    gene.clusters <- sapply(start, function(x) if(length(x)>0) 1:length(x))
    altered <- mapply(function(start, end){ mapply(function(x,y)x:y, x=start, y=end) } , start=start, end=end)
    gene.clusters <- mapply(function(start, end, gene.clusters) {
      mapply(function(x,y,z) rep(z,length(x:y)), x=start, y=end, z=gene.clusters) },
                            start=start, end=end, gene.clusters=gene.clusters)
    altered <- sapply(altered, unlist)

    gene.clusters <- sapply(gene.clusters, function(x) as.vector(unlist(x)))
    genes.altered <- mapply(function(x,y) {x<-rep(0, length(x)); x[y]<-1;x}, x=obs, y=altered)	
    gene.clusters <- mapply(function(x,y,z) {w<-rep(NA, length(x));if(length(y)>0) w[y]<-z;w}, x=obs,
                            y=altered, z=gene.clusters)
    
    ## browser() ## the next one is the one that breaks when a single chromosome
    ## because there is something weird with "gene.clusters".
    cluster.means <- mapply(function(x,y) sign(ave(x,y)), x=obs, y=gene.clusters)
    size <- lapply(sapply(object,"[", 1), length)
    Chrom <- mapply(function(x,y) rep(paste("Chrom", y),x), x=size, y=1:nchrom)
    res <- mapply(function(x,y,z,w) data.frame(Chromosome=w, x,Gain.Loss=y*z), 
                  x=obs, y=genes.altered, z=cluster.means, w=Chrom, SIMPLIFY=FALSE)
    
    res <- do.call("rbind", res)
    medians.gl <- tapply(res$x, res$Gain.Loss, median)
    medians.state <- rep(NA, length(res$x))
    medians.state[res$Gain.Loss == 1] <- medians.gl[3]
    medians.state[res$Gain.Loss == 0] <- medians.gl[2]
    medians.state[res$Gain.Loss == -1] <- medians.gl[1]
    out <- list()
    out$segm <- list()
    out$segm[[1]] <- cbind(Observed = res$x, Smoothed = medians.state,
                           State = res$Gain.Loss)
    out$chrom.numeric <- chrom.numeric
    out <- add.names.as.attr(out, array.names)

    class(out) <- c("adacgh.generic.out", "summaryACE")
    attr(out, "aceFDR.for.output") <- aceFDR.for.output
    return(out)
        ##         class(res) <- c("summary.ACE", "CGH.ACE.summary")
## 	rownames(res) <- 1:nrow(res)
##         ncr <- ncol(res) - 1
##         attr(res, "aceFDR.for.output") <- aceFDR.for.output
##         class(res) <- c("summary.ACE", "CGH.ACE.summary", "adacgh.generic.out")
## 	res
    }



summary.ACE.array <- function(object, fdr=NULL, html = TRUE,
                              outhtml ="ace.fdrtable.html", ...) {

    chrom.numeric <- object$chrom.numeric
    array.names <- object$array.names
    object$chrom.numeric <- NULL
    object$array.names <- NULL
    
    nchrom <- length(object[[1]])
                                        #Recalculate FDR for multiple arrays
    FDR.table <- lapply(object, get.FDR, nchrom=nchrom)
        
    genes <- apply(do.call("cbind",lapply(FDR.table, function(x)x[,2])), 1, sum)
    ACEPgene <- object[[1]][[1]]$ACEPgene
    FDR <- mapply(function(A, b) {(1-A) / b}, A=ACEPgene, b=genes)
###Total of genes. Calculated as in original java function. Could be missing values?	
    tot.genes <- sum(sapply(object[[1]], function(x) length(x[[1]])))*length(object)
    FDR <- FDR*tot.genes
    FDR.table <- cbind(1:length(ACEPgene), genes, round(FDR,4))
    colnames(FDR.table) <- c("Index", "Genes", "FDR")
    
    if (is.null(fdr)) fdr <- 0.15
    
    index <- which.min(abs(FDR - fdr))
    print(FDR.table)
    print(paste("Selected index", index))
    cat(FDR.table[index, 3], file ="aceFDR")
    aceFDR.for.output <- FDR.table[index, 3] ## xx: oh yes, ugly
    
    if(html) ace.fdr.html.table(FDR.table, outhtml, index)
    
    out <- list()
    out$segm <- list()
    res<-list()
    ## nextis no longer needed, as we pass array.names as component but
    ## all this is likely to be soon deprecated or, at least, little used
    array.names <- rep(NA, length(object))
    for (i in 1:length(object)) {
                                        #Recover altered genes at the FDR level
      start <- sapply(object[[i]], "[", "first")
      end <- sapply(object[[i]], "[", "last") 
      called <- sapply(object[[i]], "[", "called")
                                        #Recover observations to get the sign of their averages
      obs <- sapply(object[[i]], "[", 1)
      array.names[i] <- names(object[[i]][[1]][1])
      called <- sapply(called, function(x, index) x[index,], index=index)
      called <- lapply(called, function(x) x == 1)   
      ## Select starting and ending points of called genes
      start <- mapply(function(start, called) { start[called]}, start=start, called=called)
      end <- mapply(function(end, called) { end[called]}, end=end, called=called)
####Index for every cluster of genes
            gene.clusters <- sapply(start, function(x) if(length(x)>0) 1:length(x))
            altered <- mapply(function(start, end){ mapply(function(x,y)x:y, x=start, y=end) } , start=start, end=end)
            gene.clusters <- mapply(function(start, end, gene.clusters) {
                mapply(function(x,y,z) rep(z,length(x:y)), x=start, y=end, z=gene.clusters) },
                                    start=start, end=end, gene.clusters=gene.clusters)
            altered <- sapply(altered, unlist)
            gene.clusters <- sapply(gene.clusters, function(x) as.vector(unlist(x)))
            genes.altered <- mapply(function(x,y) {x<-rep(0, length(x)); x[y]<-1;x}, x=obs, y=altered)	
            gene.clusters <- mapply(function(x,y,z) {
                w<-rep(NA, length(x));if(length(y)>0) w[y]<-z;w}, x=obs, y=altered, z=gene.clusters)
            cluster.means <- mapply(function(x,y) sign(ave(x,y)), x=obs, y=gene.clusters)
            size <- lapply(sapply(object[[i]],"[", 1), length)
            Chrom <- mapply(function(x,y) rep(paste("Chrom", y),x), x=size, y=1:nchrom)
            res[[i]] <- mapply(function(x,y,z,w) data.frame(Chromosome=w, x,Gain.Loss=y*z), 
                               x=obs, y=genes.altered, z=cluster.means, w=Chrom, SIMPLIFY=FALSE)
            class(res[[i]]) <-"summary.ACE"
            res[[i]] <- do.call("rbind", res[[i]])
            medians.gl <- tapply(res[[i]]$x, res[[i]]$Gain.Loss, median)
            medians.state <- rep(NA, length(res$x))
            medians.state[res[[i]]$Gain.Loss == 1] <- medians.gl[3]
            medians.state[res[[i]]$Gain.Loss == 0] <- medians.gl[2]
            medians.state[res[[i]]$Gain.Loss == -1] <- medians.gl[1]
            out$segm[[i]] <- cbind(Observed = res[[i]]$x, Smoothed = medians.state,
                                   State = res[[i]]$Gain.Loss)
            
            ## 		names(res[[i]])[2] <- array.names[i]
            ##                 names(res[[i]])[3] <- paste("Gain.Loss.", array.names[i], sep = "")
        }
        ##         class(res) <- c("summary.ACE.array", "CGH.ACE.summary")
        ##         attr(res, "aceFDR.for.output") <- aceFDR.for.output
        ##         class(res) <- c("summary.ACE.array", "CGH.ACE.summary", "adacgh.generic.out")
        ##         res
       out$chrom.numeric <- chrom.numeric
    out <- add.names.as.attr(out, array.names)

        class(out) <- c("adacgh.generic.out", "summaryACE")
    attr(out, "aceFDR.for.output") <- aceFDR.for.output
    return(out)
        
        
  }


################  Parallel over arrays only

pSegmentHMM_A <- function(x, chrom.numeric, ...) {
    ## The original HMM functions chocke if chromosome
    ## is not a sequential integer starting at 1.
    ## Non-present chrom. numbers or not starting at 1 bombs.
    ## So recode
    chrom.numeric.seq <- as.numeric(as.factor(chrom.numeric))
    
    out <- papply(data.frame(x),
                  function(z) hmmWrapper_A(z, Chrom = slave_chrom.numeric.seq),
                   list(slave_chrom.numeric.seq = chrom.numeric.seq))
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    outl <- add.names.as.attr(outl, colnames(x))
    class(outl) <- c("adacgh.generic.out","mergedHMM")
    return(outl)
}

hmmWrapper_A <- function(logratio, Chrom, Pos = NULL) {
    ## Fit HMM, and do mergeLevels
    Clone <- 1:length(logratio)
    if(is.null(Pos)) Pos <- Clone
    obj.aCGH <- create.aCGH(data.frame(logratio),
                            data.frame(Clone = Clone,
                                       Chrom = Chrom,
                                       kb = Pos))
    res <- find.hmm.states(obj.aCGH, aic = TRUE, bic = FALSE)
    hmm(obj.aCGH) <- res
    out <- ourMerge(obj.aCGH$hmm$states.hmm[[1]][, 8],
                      obj.aCGH$hmm$states.hmm[[1]][, 6]) 
    return(out)
}



pSegmentBioHMM_A <- function(x, chrom.numeric, Pos, ...) {
    out <- papply(data.frame(x),
                  function(z) BioHMMWrapper_A(z,
                                              Chrom = slave_chrom.numeric,
                                              Pos    = slave_Pos),
                   list(slave_chrom.numeric = chrom.numeric,
                        slave_Pos = Pos)
                   )

    te <- unlist(unlist(lapply(out, function(x) inherits(x, "try-error"))))
    if(any(te)) {
        m1 <- "The BioHMM code occassionally crashes (don't blame us!)."
        m2 <- "You can try rerunning it a few times."
        m3 <- "You can also tell the original authors that you get the error(s): \n\n  "
        mm <- paste(m1, m2, m3, paste(out[which(te)], collapse = "    \n   "))
        caughtError(mm)
    }
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    outl$pos <- Pos
    outl <- add.names.as.attr(outl, colnames(x))
    class(outl) <- c("adacgh.generic.out","mergedBioHMM")
    return(outl)
}

BioHMMWrapper_A <- function(logratio, Chrom, Pos) {
  logratio <- matrix(logratio, ncol=1)
  uchrom <- unique(Chrom)
##  obs <- NULL
  smoothed <- NULL
  for (ic in uchrom) {
      ydat <- logratio[Chrom == ic]
      n <- length(ydat)
      res <- try(myfit.model(sample = 1, chrom = ic, dat = matrix(ydat, ncol = 1),
                       datainfo = data.frame(Name = 1:n, Chrom = rep(ic, n),
                       Position = Pos[Chrom == ic])))
      if(inherits(res, "try-error")) {
          return(res)
      } else {
          smoothed <- c(smoothed, res$out.list$mean)
      }
  }
  out <- ourMerge(logratio, smoothed)
  ##out <- ourMerge(obs, smoothed)
  return(out)
}




pSegmentDNAcopy_A <- function(x, chrom.numeric, mergeSegs = TRUE, smooth = TRUE,
                            alpha=0.01, nperm=10000, kmax=25, nmin=200,
                            eta = 0.05, overlap=0.25, trim = 0.025,
                            undo.prune=0.05, undo.SD=3, merge.pv.thresh =
                            1e-04, merge.ansari.sign = 0.05,
                            merge.thresMin = 0.05, merge.thresMax = 0.5, ...) {
    if (nperm == 10000 & alpha == 0.01 & eta == 0.05) {
        sbdry <- default.DNAcopy.bdry
    } else {
        max.ones <- floor(nperm * alpha) + 1
        sbdry <- getbdry(eta, nperm, max.ones)
    }
  sbn <- length(sbdry)
    datalist <- data.frame(x)
    papply_common <- list(slave_cnum         = chrom.numeric,
                          slave_alpha        = alpha,
                          slave_nperm        = nperm,
                          slave_kmax         = kmax,
                          slave_nmin         = nmin,
                          slave_overlap      = overlap,
                          slave_trim         = trim,
                          slave_undo.prune   = undo.prune,
                          slave_undo.SD      = undo.SD,
                          slave_sbdry        = sbdry,
                          slave_sbn          = sbn,
                          slave_merge        = mergeSegs,
                          slave_pv_th_merge  = merge.pv.thresh,
                          slave_ansari_sign  = merge.ansari.sign,
                          slave_merge_tmin   = merge.thresMin,
                          slave_merge_tmax   = merge.thresMax,
                          slave_smooth       = smooth)
  ## up to here, we are in the master
    papfunc <- function(data) {
        if(slave_smooth)
            data <- internalSmoothCNA_A(data,
                                        chrom.numeric = slave_cnum,
                                        smooth.region = 2, outlier.SD.scale = 4,
                                        smooth.SD.scale = 2, trim = 0.025)
        outseg <-
            internalDNAcopy_A(data,
                              chrom.numeric = slave_cnum,      
                              alpha =         slave_alpha,     
                              nperm =         slave_nperm,    
                              kmax =          slave_kmax,      
                              nmin =          slave_nmin,      
                              overlap =       slave_overlap,   
                              trim =          slave_trim,      
                              undo.prune =    slave_undo.prune,
                              undo.SD =       slave_undo.SD,   
                              sbdry =         slave_sbdry,     
                              sbn =           slave_sbn)
        if(!slave_merge) {
            return(outseg)
        } else {
            outmerge <- ourMerge(outseg[, 1], outseg[, 2],
                                   merge.pv.thresh = slave_pv_th_merge,
                                   merge.ansari.sign = slave_ansari_sign,
                                   merge.thresMin = slave_merge_tmin,
                                   merge.thresMax = slave_merge_tmax)
            return(outmerge)
        }
    }
    papout <- papply(datalist,
                     papfunc,
                     papply_common)

    outl <- list()
    outl$segm <- papout
    outl$chrom.numeric <- chrom.numeric
    outl <- add.names.as.attr(outl, colnames(x))
    class(outl) <- "DNAcopy" ## why not adacgh.generic.out if not mergeSegs?? FIXME!!!
    if(mergeSegs) class(outl) <- c(class(outl), "adacgh.generic.out")
    return(outl)
}


internalSmoothCNA_A <- function(acghdata, chrom.numeric,
                              smooth.region = 2, outlier.SD.scale = 4,
                              smooth.SD.scale = 2, trim = 0.025) {
    ## this is just the original smoothCNA funct. adapted to use
    ## a single array and to be parallelized and fed to internalDNAcopy
##     cat("\n DEBUG: STARTING internalSmoothCNA_A \n")
    
    chrom <- chrom.numeric
    uchrom <- unique(chrom)
    genomdat <- acghdata
    ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
    trimmed.SD <- sqrt(adacgh_trimmed.variance(genomdat[ina], trim))
    outlier.SD <- outlier.SD.scale * trimmed.SD
    smooth.SD <- smooth.SD.scale * trimmed.SD
    
    k <- smooth.region
    for (i in uchrom) {
        ina <-
            which(!is.na(genomdat) & !(abs(genomdat) == Inf) & chrom == i)
        n <- length(genomdat[ina])
        smoothed.data <-
            sapply(1:n,
                   function(i, x, n, nbhd, oSD, sSD) {
                       xi <- x[i]
                       nbhd <- i + nbhd
                       xnbhd <- x[nbhd[nbhd > 0 & nbhd <= n]]
                       if (xi > max(xnbhd) + oSD) 
                           xi <- median(c(xi, xnbhd)) + sSD
                       if (xi < min(xnbhd) - oSD) 
                           xi <- median(c(xi, xnbhd)) - sSD
                       xi
                   },
                   genomdat[ina], n, c(-k:-1, 1:k), outlier.SD, smooth.SD)
        acghdata[ina] <- smoothed.data
    }
    acghdata
}

internalDNAcopy_A <- function(acghdata,
                            chrom.numeric,
                            sbdry,
                            sbn,
                            alpha,
                            nperm,
                            kmax,
                            nmin,
                            overlap, 
                            trim,
                            undo.prune,
                            undo.SD) {
    ## tries to follow the original "segment"
    
    uchrom <- unique(chrom.numeric)
    data.type <- "logratio"
    p.method <- "hybrid"
    window.size <- NULL
    undo.splits <- "none"
    genomdati <- acghdata
    min.width <- 2
    ina <- which(!is.na(genomdati) & !(abs(genomdati)==Inf))
    
    ## The code allows for dealing with NA and Inf, but would need to
    ## adjust other functions (as different arrays would have different
    ## length of pos, genenames, etc. So for now stop:
    if (length(ina) != length(genomdati))
        stop("Either an NA or an infinite in the data")

    genomdati <- genomdati[ina]
    trimmed.SD <- sqrt(adacgh_trimmed.variance(genomdati, trim))
    chromi <- chrom.numeric[ina]
    sample.lsegs <- NULL
    sample.segmeans <- NULL
##     cat("\n DEBUG: internalDNAcopy_A: before loop uchrom \n")

    for (ic in uchrom) {
##         cat("\n DEBUG: internalDNAcopy_A: before changepoints \n")
##         browser()
        segci <- adacgh_changepoints(genomdati[chromi==ic],
                              data.type = "logratio",
                              alpha = alpha, sbdry = sbdry, sbn = sbn,
                              nperm = nperm, p.method = p.method,
##                               window.size = window.size, 
##                               overlap = overlap,
                              kmax = kmax, nmin = nmin,
                              trimmed.SD = trimmed.SD,
                              undo.splits = undo.splits,
                              undo.prune = undo.prune,
                              undo.SD = undo.SD, verbose = 2,
                              min.width = min.width)
##         cat("\n DEBUG: internalDNAcopy_A: end of changepoints \n")

        sample.lsegs <- c(sample.lsegs, segci$lseg)
        sample.segmeans <- c(sample.segmeans, segci$segmeans)

    }
##     cat("\n DEBUG: internalDNAcopy_A: after loop uchrom \n")

    if(length(sample.lsegs) != length(sample.segmeans))
        stop("Something terribly wrong: length(sample.lsegs) != length(sample.segmeans).")
    stretched.segmeans <- rep(sample.segmeans, sample.lsegs)
    stretched.state    <- rep(1:length(sample.lsegs), sample.lsegs)
    return(cbind(Observed = genomdati, Predicted = stretched.segmeans,
                 State = stretched.state))
}


##ace.analysisP_C <-function(x) {
## ace.analysisP_C <-function(x, coefs, Sdev, array.names) {
##   ace.analysis.C(x, coefs, Sdev, array.names)
## }




ACE_C <- function(x, Chrom, coefs = file.aux, Sdev=0.2, echo=FALSE) {
  
####### x is log2.ratio
####### Chrom ---MUST BE NUMERIC-- 
###### only for 1 array at a time
  if (!is.numeric(Chrom)) {
    stop("Chromosome variable must be numeric")
  }
  
  array.names <- colnames(x)
  
  nchrom <- length(unique(Chrom))
  if(is.null(dim(x)) || (dim(x)[2]==1)) {
    genes <- split(x, Chrom)
    first.estimate <- papply(genes,
                             function(z) ace.analysis.C(z,
                                                        coefs = slave_coefs,
                                                        Sdev = slave_Sdev,
                                                        array.names = slave_array.names),
                              list(slave_coefs = coefs,
                                   slave_Sdev = Sdev,
                                   slave_array.names = array.names))
    Sdevs.estimate <- papply(first.estimate, sd.ACE.analysis)
    Sdevs <- unlist(lapply(Sdevs.estimate, "[", 1))
    Sdevs <- mean(Sdevs[Sdevs>0])
    if(is.nan(Sdevs)) Sdevs <- 0
    
    res <- papply(genes, function(z) ace.analysis.C(z,
                                                    coefs = slave_coefs,
                                                    Sdev = slave_Sdev,
                                                    array.names = slave_array.names),
                              list(slave_coefs = coefs,
                                   slave_Sdev = Sdev,
                                   slave_array.names = array.names))
    class(res) <- "ACE"
  }
  else {
    
    Sdevs <- matrix(NA, ncol(x),nchrom)
    res <- list()
    for(i in 1:ncol(x)) {
      genes <- split(x[,i], Chrom)
      first.estimate <- papply(genes,
                               function(z) ace.analysis.C(z,
                                                          coefs = slave_coefs,
                                                          Sdev = slave_Sdev,
                                                          array.names = slave_array.name),
                                list(slave_coefs = coefs,
                                     slave_Sdev = Sdev,
                                     slave_array.name = array.names[i]))
      Sdevs.estimate <- papply(first.estimate, sd.ACE.analysis)
      Sdevs[i,] <- unlist(lapply(Sdevs.estimate, "[", 1))
    }
    Sdevs <- mean(Sdevs[Sdevs>0])
    if(is.nan(Sdevs)) Sdevs <- 0
    for (i in 1:ncol(x)) {
      genes <- split(x[,i], Chrom)
      
      res[[i]] <- papply(genes,
                               function(z) ace.analysis.C(z,
                                                          coefs = slave_coefs,
                                                          Sdev = slave_Sdev,
                                                          array.names = slave_array.name),
                                list(slave_coefs = coefs,
                                     slave_Sdev = Sdev,
                                     slave_array.name = array.names[i]))
      class(res[[i]]) <- "ACE"
    }
    class(res) <- "ACE.array"
  }
  res$chrom.numeric <- Chrom
  res$array.names <- array.names
  invisible(res)
}













#######################################################
#######################################################
#######################################################
###
###            Printing results and PaLS
###            
###
#######################################################
#######################################################
#######################################################


writeResults <- function(obj, ...) {
    UseMethod("writeResults")
}

## writeResults.CGH.PSW <- function(obj, acghdata, commondata, file = "PSW.output.txt", ...) {
writeResults.CGH.PSW <- function(obj, commondata, file = "PSW.output.txt", ...) {    
    write.table(cbind(commondata, obj$Data), file = file,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
}

writeResults.summaryACE <- function(obj, acghdata, commondata, file = NULL, ...) {
    if(is.null(file)) {
        file <-  paste("ACE.output.FDR=",
                       attr(obj, "aceFDR.for.output"), ".txt", sep ="")
    }
    print.adacgh.generic.results(obj, acghdata, commondata, output = file)
}

writeResults.CGH.wave <- function(obj, acghdata, commondata,
                                  file = "wavelet.output.txt", ...) {
    if(inherits(obj, "CGH.wave.merged")) {pals <- TRUE} else {pals <- FALSE}
    print.adacgh.generic.results(obj, acghdata, commondata, output = file,
                                 send_to_pals = pals)
}

writeResults.DNAcopy <- function(obj, acghdata, commondata, 
                                 file = "CBS.output.txt", ...) {
  print.adacgh.generic.results(obj, acghdata,
                               commondata, output = file)
}

writeResults.CGHseg <- function(obj, acghdata, commondata, 
                                 file = "CGHseg.output.txt", ...) {
    print.adacgh.generic.results(obj, acghdata,
                                commondata, output = file,
                                 send_to_pals = FALSE)
}

writeResults.adacghHaarSeg <- function(obj, acghdata, commondata, 
                                 file = "HaarSeg.output.txt", ...) {
    print.adacgh.generic.results(obj, acghdata,
                                 commondata, output = file)
}



writeResults.mergedHMM <- function(obj, acghdata, commondata, 
                                 file = "HMM.output.txt", ...) {
    print.adacgh.generic.results(obj, acghdata,
                                commondata, output = file)
}

writeResults.adacghGLAD <- function(obj, acghdata, commondata, 
                                 file = "GLAD.output.txt", ...) {
    print.adacgh.generic.results(obj, acghdata,
                                commondata, output = file)
}

writeResults.mergedBioHMM <- function(obj, acghdata, commondata, 
                                 file = "BioHMM.output.txt", ...) {
    print.adacgh.generic.results(obj, acghdata,
                                commondata, output = file)
}

print.adacgh.generic.results <- function(res, xcenter,
                                 commondata,
                                 output = "ADaCGH.results.txt",
                                 send_to_pals = TRUE){

    out <- data.frame(commondata)
    if(ncol(out) > 5) {
        stop("This sucks, but if your commondata has more than 5 columns, this function will blow up.")
    }

    for(i in 1:ncol(xcenter)) {
      out <- cbind(out, res$segm[[i]])
    }
    colnames(out)[(ncol(commondata) + 1):(ncol(out))] <-
        paste(rep(colnames(xcenter),rep(3, ncol(xcenter))),
              c(".Original", ".Smoothed", ".Status"),
              sep = "")

    write.table(out, file = output,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)

    if (exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv) & send_to_pals) {
        cols.look <- seq(from = 8, to = ncol(out), by = 3)

        Ids <- apply(out[, cols.look, drop = FALSE], 2,
                     function(z) commondata$ID[which( z == -1)])
        writeForPaLS(Ids, colnames(xcenter), "Lost_for_PaLS.txt")
        
        Ids <- apply(out[, cols.look, drop = FALSE], 2,
                     function(z) commondata$ID[which( z == 1)])
        writeForPaLS(Ids, colnames(xcenter), "Gained_for_PaLS.txt")

        Ids <- apply(out[, cols.look, drop = FALSE], 2,
                     function(z) commondata$ID[which( z != 0)])
        writeForPaLS(Ids, colnames(xcenter), "Gained_or_Lost_for_PaLS.txt")
    }

}

writeForPaLS <- function(alist, names, outfile) {
    ## alist: a list with as many lists as subjects; each sublist are the
    ##        genes of interest.
    ## names: subject or array names
    ## outfile: guess what? is the name of the output file

    
  if(is.array(alist) | is.matrix(alist) )
    if (dim(alist)[2] == 1) alist <- as.vector(alist)

    if(!is.list(alist) & is.vector(alist) & (length(names) == 1)) {
        ## we suppose we are dealing with a one-array data set
        alist <- list(alist)
    }
  if(length(alist) == 0) {
      write("", file = outfile)
  } else if(length(names) != length(alist)) {
      print("names are ")
      print(names)
      print("alist is ")
      print(alist)
      stop("ERROR in writeForPaLS: names and alist should have the same length")
  } else {
      write(
            unlist(
                   mapply(function(x, y) return(c(paste("#", y, sep = ""), as.character(x))),
                          alist, names)
                   ),
            file = outfile)
  }
}




##################################################
##################################################
##################################################
##################################################


##### Internal stuff, mostly for the web-based part

PSWtoPaLS <- function(x, y,
                      out = "Gained_or_Lost_for_PaLS.txt") {
  ## Create the gained or lost
  ooo <- function(u, v) mapply(function(x, y, nx)
                               return(c(paste("#", nx, sep = ""), x, y)),
                               u, v, names(u))
  write(unlist(ooo(x, y)), file = out)     
}


##### FIXME: this is all VERY ugly. We should not have code here that refers
##    to the web app, because it is a mess to change the web-app behavior, and
##    that should not be retrofitted here.

##    I do not want to break the old stuff for ADaCGH web app. So I call
##    another function for the new server. This will allow smoother changes
##    when we finally finish with the old ADaCGH web app.

caughtOtherError <- function(message) {
    png.height <- 400
    png.width  <- 400
    png.pointsize <- 10

    if(exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) {
        GDD("ErrorFigure.png", width = png.width,
               height = png.height, 
               ps = png.pointsize)
        plot(x = c(0, 1), y = c(0, 1),
             type = "n", axes = FALSE, xlab = "", ylab = "")
        box()
        text(0.5, 0.7, "There was a PROBLEM with this run.")
        dev.off()
        sink(file = "results.txt")
        cat(message)
        sink()
        sink(file = "exitStatus")
        cat("Error\n\n")
        cat(message)
        sink()
        quit(save = "no", status = 11, runLast = TRUE)
    } else if(exists(".__ADaCGH_SERVER_APPL", env = .GlobalEnv)) {
        caughtOtherError.Web(message)

    } else {
        message <- paste(message, " ", collapse = " ")
        message <- paste("There is a possible problem: ", message)
        stop(message)
    }
}

caughtError <- function(message) {
    png.height <- 400
    png.width  <- 400
    png.pointsize <- 10

    if(exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) {
        GDD("ErrorFigure.png", width = png.width,
               height = png.height, 
               ps = png.pointsize)
        plot(x = c(0, 1), y = c(0, 1),
             type = "n", axes = FALSE, xlab = "", ylab = "")
        box()
        text(0.5, 0.7, "There was a PROBLEM with this run.")
        dev.off()
        sink(file = "results.txt")
        cat(message)
        sink()
        sink(file = "exitStatus")
        cat("Error\n\n")
        cat(message)
        sink()
        quit(save = "no", status = 11, runLast = TRUE)
    } else if(exists(".__ADaCGH_SERVER_APPL", env = .GlobalEnv)) {
        caughtOtherPackageError.Web(message)
    } else {
        message <- paste("This is a known problem in a package we depend upon. ", message)
        stop(message)
    }
}



caughtOurError <- function(message) {
    png.height <- 400
    png.width  <- 400
    png.pointsize <- 10

    if(exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) {
        GDD("ErrorFigure.png", width = png.width,
               height = png.height, 
               ps = png.pointsize)
        plot(x = c(0, 1), y = c(0, 1),
             type = "n", axes = FALSE, xlab = "", ylab = "")
        box()
        text(0.5, 0.7, "There was a PROBLEM with the code.")
        text(0.5, 0.5,
             "Please let us know (send us the URL),")
        
        text(0.5, 0.3, "so that we can fix it.")
        dev.off()
        sink(file = "results.txt")
        cat(message)
        sink()
        sink(file = "exitStatus")
        cat("Error\n\n")
        cat(message)
        sink()
        quit(save = "no", status = 11, runLast = TRUE)
    } else if(exists(".__ADaCGH_SERVER_APPL", env = .GlobalEnv)) {
        caughtOurError.Web(message)
    } else {
        message <- paste("It looks like you found a bug. Please let us know. ", message)
        stop(message)
    }
}
    



caughtOtherError.Web <- function(message) {
    mpi.clean.quit.Web()
    sink(file = "R_Error_msg.txt")
    cat(message)
    cat("\n")
    sink()
    sink(file = "R_Status.txt")
    cat("Other Error\n\n")
    sink()
    quit(save = "no", status = 11, runLast = FALSE)
}

caughtOtherPackageError.Web <- function(message) {
    mpi.clean.quit.Web()
    message <- paste("This is a known problem in a package we depend upon. ",
                     message)
    sink(file = "R_Error_msg.txt")
    cat(message)
    cat("\n")
    sink()
    sink(file = "R_Status.txt")
    cat("Other Error\n\n")
    sink()
    quit(save = "no", status = 11, runLast = FALSE)
}



mpi.clean.quit.Web <- function() {
    if (is.loaded("mpi_initialize")){ 
        if (mpi.comm.size(1) > 0){ 
            try(mpi.close.Rslaves() , silent = TRUE)
        } 
    }
    try(mpi.exit(), silent = TRUE)
}


caughtOurError.Web <- function(message) {
    mpi.clean.quit.Web()
    message <- paste("There was a problem with our code. Please let us know.\n", 
                     message)
    sink(file = "R_Error_msg.txt")
    cat(message)
    cat("\n")
    sink()
    sink(file = "R_Status.txt")
    cat("Our Error\n\n")
    sink()
    quit(save = "no", status = 11, runLast = FALSE)

}


caughtUserError.Web <- function(message) {
    mpi.clean.quit.Web()
    message <- paste("There was a problem with something you did.\n",
                     "Check the error message, your data and options and try again.\n",
                     message, "\n")
    sink(file = "R_Error_msg.txt")
    cat(message)
    cat("\n")
    sink()
    sink(file = "R_Status.txt")
    cat("User Error\n\n")
    sink()
    quit(save = "no", status = 11, runLast = FALSE)
}


doCheckpoint <- function(num) {
##    checkpoint.num.new <- num
    save.image()
##    checkpoint.num <<- num
    sink("checkpoint.num")
    cat(num)
    sink()
    return(num)
}




my.html.data.frame <- function (object, first.col = "Name",
                             file = paste(first.word(deparse(substitute(object))), 
                             "html", sep = "."), append = FALSE, link = NULL, linkCol = 1, 
                             linkType = c("href", "name"), ...) 
{
    ## modifying html, from Hmisc: Their function always has first column
    ## named "Name". I allow to pass a name.
   
    linkType <- match.arg(linkType)
    x <- format.df(object, numeric.dollar = FALSE, ...)
    adj <- attr(x, "col.just")
    if (any(adj == "r")) 
        for (i in seq(along = adj)[adj == "r"]) x[, i] <- paste("<div align=right>", 
            x[, i], "</div>", sep = "")
    if (length(r <- dimnames(x)[[1]])) 
        x <- cbind(first.col = r, x)
    colnames(x)[1] <- first.col
    cat("<TABLE BORDER>\n", file = file, append = append)
    cat("<tr>", paste("<td>", dimnames(x)[[2]], "</td>", sep = ""), 
        "</tr>\n", sep = "", file = file, append = file != "")
    if (length(link)) 
        x[, linkCol] <- ifelse(link == "", x[, linkCol], paste("<a ", 
            linkType, "=\"", link, "\">", x[, linkCol], "</a>", 
            sep = ""))
    for (i in 1:nrow(x)) cat("<tr>", paste("<td>", x[i, ], "</td>", 
        sep = ""), "</tr>\n", sep = "", file = file, append = file != 
        "")
    cat("</TABLE>\n", file = file, append = file != "")
    structure(list(file = file), class = "html")
}






    
######################################################

#########  Imagemap stuff


imClose3 <- function (im) {
    ## prevent all the "Closing PNG device ..."
    dev.off(im$Device)
}

imagemap3 <- function(filename,width=480,height=480,
                      title='Imagemap from R', ps = 12){

    GDD(file = paste(filename,".png",sep=''),w=width, h=height,
        type = "png", ps = ps)	  
	  
    im <- list()
    im$Device <- dev.cur()
    im$Filename=filename
    im$Height=height
    im$Width=width
    im$Objects <- list()
    im$HTML <- list()
    im$title <- title
    
    class(im) <- "imagemap"
    im
}

createIM2 <- function(im, file = "", imgTags = list(),
                      title = "Genome View") {
    cat(paste("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n",
              "<html> <head> <title>", title, "</title></head><body>"),
        file = file)
    cat(buildIM(im, imgTags), sep = "\n", file = file, append = TRUE)
    cat("</body></html>", file = file, append = TRUE)
}


plot.adacgh.nonsuperimpose <- function(res, chrom,  main, colors,
                                       ylim, geneNames, idtype, organism,
                                       geneLoc, html_js, imgheight,
                                       genomewide_plot= FALSE,
                                       chromsplot = NULL) {
  cat("\n plot.adacgh.nonsuperimpose: Doing sample ", main, "\n")
  if(genomewide_plot) {
    plot.adacgh.genomewide(res, chrom, geneNames, imgheight,
                           main, colors,
                           ylim, geneLoc)
  }
  
  plot.adacgh.chromosomewide(res, chrom,
                             geneNames,
                             imgheight,
                             main, colors,
                             ylim, idtype, organism, geneLoc,
                             html_js = html_js,
                             chromsplot = chromsplot)
}

plot.adacgh.genomewide <- function(res, chrom,
                                   geneNames,
                                   imgheight,
                                   main = NULL,
                                   colors = c("orange", "red", "green", "blue", "black"),
                                   ylim = NULL,
                                   geneLoc = NULL) {
    warning(paste("This function is likely to be soon deprecated.",
            "With huge arrays, genomewide plots make little sense."))
    
    cat("\n        plot.adacgh.genomewide: Doing sample ", main, "\n")
    pch <- 20
    im1 <- mapGenomeWideOpen(main)
    nameIm <- main
    logr <- res[, 1]
    smoothdat <- res[, 2]
    if(is.null(geneLoc)) {
        simplepos <- (1:length(logr))
    } else {
        ## geneLoc is withing chromosome,
        ## thus, we need some absolute, increasing pos.
        lchr <- tapply(geneLoc, chrom, length)
        mchr <- cumsum(as.numeric(tapply(geneLoc, chrom, max))) ## not very efficient
        sumpos <- rep(c(0, mchr[-length(mchr)]),
                      lchr)
        simplepos <- geneLoc + sumpos
    }
                        
    res.dat <- res[, 3]
    col <- rep(colors[1],length(res.dat))
    col[which(res.dat == -1)] <- colors[3]
    col[which(res.dat == 1)] <- colors[2]
  
    plotGenomeWide(logr, simplepos, col, main, pch, ylim, colors)
    im1 <- mapLinkChrom(logr, simplepos, chrom, nameIm, im1)
   
    lines(smoothdat ~ simplepos, col=colors[4],
          lwd = 2)
  
    mapGenomeWideClose(nameIm, im1)
}



plot.adacgh.chromosomewide <- function(res, chrom,
                                       geneNames, imgheight,
                                       main = NULL,
                                       colors = c("orange", "red", "green", "blue", "black"),
                                       ylim = NULL,
                                       idtype = idtype,
                                       organism = organism,
                                       geneLoc = NULL,
                                       html_js = FALSE,
                                       chromsplot = NULL) {

  pixels.point <- 3
  pch <- 20
  col <- rep(colors[1],length(res[, 3]))
  col[which(res[, 3] == -1)] <- colors[3]
  col[which(res[, 3] == 1)] <- colors[2]
  nameIm <- main
  if(is.null(chromsplot)) chrom.nums <- unique(chrom)
  else chrom.nums <- chromsplot
  cat("\n  plot.adacgh.chromosomewide: doing sample ", main, "\n")
  
  for(cnum in 1:length(chrom.nums)) {
    cat("\n        plot.adacgh.chromosomewide: doing chromosome  ", cnum, "\n")
    indexchr <- which(chrom == chrom.nums[cnum])
    ccircle <- NULL
    simplepos <- if(is.null(geneLoc)) (1:length(indexchr)) else geneLoc[indexchr]
    ## Formerly mapChromOpen
    chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
    chrwidth <- max(chrwidth, 800)
    im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     height = imgheight, width = chrwidth,
                     ps = 12)
    ##
    ## Formerly plotChromWide()
    par(xaxs = "i")
    par(mar = c(5, 5, 5, 5))
    par(oma = c(0, 0, 0, 0))
    
    plot(res[indexchr,1] ~ simplepos, col=col[indexchr], cex = 1,
         xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE,
         main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
         pch = pch, ylim = ylim)
    box()
    axis(2)
    abline(h = 0, lty = 2, col = colors[5])
    rug(simplepos, ticksize = 0.01)
    lines(res[indexchr, 2] ~ simplepos,
          col = colors[4], lwd = 2, type = "l")
    dummy.coord <- usr2png(cbind(c(2, 0), c(0, 0)), im2)
    cc1.r <- max(abs(dummy.coord[1, 1]  - dummy.coord[2, 1]), 4)
    ccircle <- rbind(t(usr2png(cbind(simplepos, res[indexchr, 1]), im2)),
                     rep(cc1.r, length(simplepos)))
    nameChrIm <- paste("Chr", chrom.nums[cnum], "@", nameIm, sep ="")
    write(ccircle, file = paste("pngCoordChr_", nameChrIm, sep = ""),
          sep ="\t", ncolumns = 3)
    calcnarrays <- ncol(ccircle)/length(geneNames[indexchr])
    ## what are the next three lines for here?? FIXME
    if(!exists("arraynums")) arraynums <- 1
    if(calcnarrays != arraynums)
      stop("Serious problem: number of arrays does not match")
    write(rep(as.character(geneNames[indexchr]), arraynums), 
          file = paste("geneNamesChr_", nameChrIm, sep = ""))
    imClose3(im2)
    if(html_js) 
      system(paste(.python.toMap.py, nameChrIm, 
                   idtype, organism, sep = " "))
  }        ## looping over chromosomes
}


plot.gw.superimp <- function(res, chrom, ylim, geneNames, imgheight, main = NULL,
                             colors = c("orange", "red", "green", "blue", "black"),
                             geneLoc = NULL) {
  
  pch <- ""
  arraynums <- length(res)
  im1 <- mapGenomeWideOpen(main)
  nameIm <- main
  nfig <- 1
  
  for (arraynum in 1:arraynums) {
    cat("\n      plot.gw.superimp: doing array ", arraynum, "\n")
    logr <- res[[arraynum]][, 1]
    smoothdat <- res[[arraynum]][, 2]
    if(is.null(geneLoc)) {
      simplepos <- (1:length(logr))
    } else {
      ## geneLoc is withing chromosome,
      ## thus, we need some absolute, increasing pos.
      lchr <- tapply(geneLoc, chrom, length)
      mchr <- cumsum(as.numeric(tapply(geneLoc, chrom, max))) ## not very efficient
      sumpos <- rep(c(0, mchr[-length(mchr)]),
                    lchr)
      simplepos <- geneLoc + sumpos
    }
    res.dat <- res[[arraynum]][, 3]
    col <- rep(colors[1],length(res.dat))
    col[which(res.dat == -1)] <- colors[3]
    col[which(res.dat == 1)] <- colors[2]
    
    if(nfig == 1) {
      plotGenomeWide(logr, simplepos, col, main, pch, ylim, colors)
      im1 <- mapLinkChrom(logr, simplepos, chrom, nameIm, im1)
    }
    
    lines(smoothdat ~ simplepos,
          col = colors[4], lwd = 2, type = "l")
    nfig <- nfig + 1
    par(new = TRUE)
  }
  mapGenomeWideClose(nameIm, im1)
}

plot.cw.superimpA <- function(res, chrom,
                              geneNames,
                              main = "All_arrays",
                              colors = c("orange", "red", "green", "blue", "black"),
                              ylim =NULL, 
                              idtype = idtype, organism = organism,
                              geneLoc = NULL,
                              html_js = html_js,
                              imgheight = imgheight) {
  ## For superimposed: one plot per chr
  pch <- ""
  arraynums <- length(res)
  nameImage <- main
  pixels.point <- 3
  chrom.nums <- unique(chrom)
  datalist <- list()
  for(cnum in 1:length(chrom.nums)) {
    indexchr <- which(chrom == chrom.nums[cnum])
    datalist[[cnum]] <- list()
    datalist[[cnum]]$indexchr <- indexchr
    datalist[[cnum]]$cnum <- cnum
    datalist[[cnum]]$thiscn <- chrom.nums[cnum]
    datalist[[cnum]]$resl <- lapply(res, function(w) w[indexchr, , drop = FALSE])
    if(!is.null(geneLoc))
      datalist[[cnum]]$posn <- geneLoc[indexchr]
    else
      datalist[[cnum]]$posn <- NULL
  }
  pappl_common <- list(arraynums = arraynums, nameImage = main,
                       geneNames = geneNames, imgheight = imgheight,
                       pch = pch, ylim = ylim, html_js = html_js,
                       idtype = idtype, organism = organism)
  
  funp <- function(z) {
                   ## arraynums, nameImage,
                   ## geneNames, imgheight,
                   ## pch, ylim, html_js,
                   ## idtype, organism) {
    if(is.null(z)) return()
    indexchr <- z$indexchr
    ccircle <- NULL
    thiscn <- z$thiscn
    cnum <- z$cnum
    ##         environment(mapChromOpenA) <- environment()
    ##         im2 <- mapChromOpenA()
    pixels.point <- 3
    chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
    chrwidth <- max(chrwidth, 800)
    im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameImage, sep =""),
                     height = imgheight, width = chrwidth,
                     ps = 12)
    nfig <- 1
    for(arraynum in 1:arraynums) { ## first, plot the points
      cat("\n      plot.cw.superimpA: doing array ", arraynum, "\n")
      logr <- z$resl[[arraynum]][, 1]
      res.dat <- z$resl[[arraynum]][, 3]
      smoothdat <- z$resl[[arraynum]][, 2]
      col <- rep(colors[1],length(res.dat))
      col[which(res.dat == -1)] <- colors[3]
      col[which(res.dat == 1)] <- colors[2]
      simplepos <-
        if(is.null(z$posn)) (1:length(logr)) else z$posn
      
      
      if(nfig == 1) {
        ## Formerly plotChromWideA
        par(xaxs = "i")
        par(mar = c(5, 5, 5, 5))
        par(oma = c(0, 0, 0, 0))
        plot(logr ~ simplepos, col=col, cex = 1,
             xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE,
             main = paste("Chr", thiscn, "@", nameImage, sep =""),
             pch = pch, ylim = ylim)
        box()
        axis(2)
        abline(h = 0, lty = 2, col = colors[5])
        rug(simplepos, ticksize = 0.01)
        
      }
      ##             environment(pngCircleRegionA) <- environment()
      ##             ccircle <- pngCircleRegionA()
      
      dummy.coord <- usr2png(cbind(c(2, 0), c(0, 0)), im2)
      cc1.r <- max(abs(dummy.coord[1, 1]  - dummy.coord[2, 1]), 4)
      ccircle <- cbind(ccircle,
                       rbind(t(usr2png(cbind(simplepos, logr), im2)),
                             rep(cc1.r, length(simplepos))))
      ## we want all points, but only draw axes once
      points(logr ~ simplepos, col=col,
             cex = 1, pch = 20)
      lines(smoothdat ~ simplepos,
            col = colors[4], lwd = 2, type = "l")
      nfig <- nfig + 1
    }
    ##         environment(mapCloseAndPythonChromA) <- environment()
    ##         mapCloseAndPythonChromA()
    
    nameChrIm <- paste("Chr", thiscn, "@", nameImage, sep ="")
    write(ccircle, file = paste("pngCoordChr_", nameChrIm, sep = ""),
          sep ="\t", ncolumns = 3)
    calcnarrays <- ncol(ccircle)/length(geneNames[indexchr])
    if(!exists("arraynums")) arraynums <- 1
    if(calcnarrays != arraynums) {
      ##           cat("\n calcnarrays ", calcnarrays)
      ##           cat("\n arraynums ", arraynums)
      ##           cat("\n ncol(ccircle) ", ncol(ccircle))
      ##           cat("\n length((geneNames[indexchr]) ", length(geneNames[indexchr]))
      ##           cat("\n length(indexchr) ", length(indexchr), "\n")
      stop("Serious problem: number of arrays does not match")
    }
    
    write(rep(as.character(geneNames[indexchr]), arraynums), 
          file = paste("geneNamesChr_", nameChrIm, sep = ""))
    imClose3(im2)
    if(html_js)
      system(paste(.python.toMap.py, nameChrIm, 
                   idtype, organism, sep = " "))
  } ## END OF FUNP
  out <- papply(datalist, funp, pappl_common)
                ## function(z) funp(z,
                ##    arraynums = arraynums, nameImage = main,
                ##    geneNames = geneNames, imgheight = imgheight,
                ##    pch = pch, ylim = ylim, html_js = html_js,
                ##    idtype = idtype, organism = organism
                ##    ))
  ## papply_commondata = pappl_common)
}

mapGenomeWideOpen <- function(main) {
    nameIm <- main
    imgheight <- 500
    imwidth <- 1600
    im1 <- imagemap3(nameIm, height = imgheight,
                     width = imwidth, ps = 12)
    return(im1)
}

mapGenomeWideClose <- function(nameIm, im1) {
    createIM2(im1, file = paste(nameIm, ".html", sep = ""))
    imClose3(im1)
}

plotGenomeWide <- function(logr, simplepos, col, main, pch, ylim, colors) {
    plot(logr ~ simplepos, col= col, ylab = "log ratio",
         xlab ="Chromosome location", axes = FALSE, cex = 0.7, main = main,
         pch = pch, ylim = ylim)
    box()
    rug(simplepos, ticksize = 0.01)
    axis(2)
    abline(h = 0, lty = 2, col = colors[5])
}


mapLinkChrom <- function(logr, simplepos, chrom, nameIm, im1) {
    ## 1) add the vertical chromosome lines
    ## 2) html map: add links to chromosome-wide figures
    LimitChr <- tapply(simplepos,
                       chrom, max)
    abline(v=LimitChr, col="grey", lty=2)
    
    chrom.nums <- unique(chrom)
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)
    
    lxs <- c(1, LimitChr)
    maxlr <- max(logr)
    minlr <- min(logr)
    nd <- 1:length(LimitChr)
    xleft <- lxs[nd]
    names(xleft) <- 1:length(xleft)
    xright <- lxs[nd + 1]
    f1 <- function(xleft, xright, nd)
        imRect(xleft, maxlr, xright, minlr - 10,
               title = paste("Chromosome", nd),
               alt = paste("Chromosome", nd),
               href= paste("Chr", nd, "@", nameIm, ".html", sep =""))
    rectslist <- mapply(f1, xleft, xright, nd, SIMPLIFY=FALSE)
    for(ll in 1:length(rectslist))
        addRegion(im1) <- rectslist[[ll]]
    return(im1)
}



constructSegmObj <- function(x, chrom.numeric, data, Pos) {
    ## x: our generic output object

    l1 <- list()
    l1$data <- data.frame(chrom = chrom.numeric,
                          maploc = Pos,
                          data)
    l1$output <- list()
    cuniq <- unique(chrom.numeric)
    for(i in 1:length(x)) { ## over arrays
        running.start <- 0
        z <- x[[i]][, 2]
        for(chr in cuniq) {
            thispos <- Pos[chrom.numeric == chr]
            y <- z[chrom.numeric == chr]
            nelem <- length(y)
            poschange <- which(diff(c(NA, y)) != 0)
            this.start <- c(thispos[1], thispos[poschange])
            this.end <- c(thispos[poschange - 1], thispos[nelem])
            l1$output$ID <- c(l1$output$ID,
                rep(colnames(data)[i], length(this.start)))
            l1$output$chrom <- c(l1$output$chrom, rep(chr, length(this.start)))
            l1$output$loc.start <- c(l1$output$loc.start,
                                     this.start)
            l1$output$loc.end <- c(l1$output$loc.end,
                                   this.end)
            l1$output$seg.mean <- c(l1$output$seg.mean,
                                    y[c(1, poschange)])
        }
    }
    l1$output <- as.data.frame(l1$output)
    return(l1)
}



doMCR <- function(x, chrom.numeric, data,
                  MCR.gapAllowed = 500,
                  MCR.alteredLow = 0.03,
                  MCR.alteredHigh = 0.97,
                  MCR.recurrence = 75,
                  fsink = "results.txt",
                  hsink = "mcr.results.html",
                  Pos,  ...) {
    if(ncol(data) > 1) {
        delta <- 0.05
        res <- constructSegmObj(x, chrom.numeric, data, Pos)
        ## Hoping over the bugs in cghMCR
        tryms <- try({
            cghmcr <- cghMCR(res,
                             gapAllowed = MCR.gapAllowed,
                             alteredLow = MCR.alteredLow,
                             alteredHigh = MCR.alteredHigh,
                             recurrence = MCR.recurrence)
            mcrs <- MCR(cghmcr)
            if(is.null(dim(mcrs))) {
                mcrs <- NA
            } else {
                if(dim(mcrs)[1] == 0) mcrs <- NA
            }
        })
        if((class(tryms) == "try-error") || is.na(mcrs)) {
            sink(file = hsink)
            cat("\n<p> No common minimal regions found.</p>\n")
            sink()
            sink(file = fsink)
            cat("\n\n\nMinimal common regions\n")
            cat("\n No common minimal regions found.\n")
            sink()
        } else {
            msamples <- sapply(mcrs[, 9],
                               function(x) length(unlist(strsplit(x, ",")))) 
            mcrselect <- which(msamples > 1)
            mcrs <- mcrs[mcrselect, ,drop = FALSE]
            mcrsc <- data.frame(chromosome = mcrs[, 1],
                                samples = mcrs[, 9],
                                mcr.start  = as.numeric(mcrs[, 7]),
                                mcr.end  = as.numeric(mcrs[, 8]))
            mcrsc <- mcrsc[ order(as.numeric(as.character(mcrsc$chromosome)),
                                  mcrsc$mcr.start, mcrsc$mcr.end), ]
            sink(file = hsink)
            if (nrow(mcrsc) == 0) {
                cat("\n<p> No common minimal regions found.</p>\n")
            } else {
                my.html.data.frame(mcrsc, first.col = "Case",
                                file = hsink, append = TRUE)
            }
            sink()
            
            sink(file = fsink)
            cat("\n\n\nMinimal common regions\n")
            if (nrow(mcrsc) == 0)
                cat("\n No common minimal regions found.\n")
            else 
                print(mcrsc)
            sink()
        }
    }
}
   


papply0 <- function (arg_sets, papply_action, papply_commondata = list(), 
    show_errors = TRUE, do_trace = FALSE, also_trace = c()) 
{
### this is a crippled version of papply, for running sequentially.
    ### I'm getting lots of mpi problems in a few cases.
    ### instead of rewriting calls, I close mpi, and use papply0.
    ### If these issues ever get fixed, then use papply.
    
    if (!is.list(arg_sets)) {
        print("1st argument to papply must be a list")
        return(NULL)
    }
    if (!is.function(papply_action)) {
        print("2nd argument to papply must be a function")
        return(NULL)
    }
    if (!is.list(papply_commondata)) {
        print("3rd argument to papply must be a list")
        return(NULL)
    }
    papply_also_trace <- also_trace
    run_parallel <- 0
    if (run_parallel != 1) {
        print("Running serial version of papply\n")
        attach(papply_commondata)
        results <- lapply(arg_sets, papply_action)
        detach(papply_commondata)
        return(results)
    }
}


##### This is straight from aCGH, but fixing the problem of wilcox.text(exact = T)

mergeLevelsB <- function(vecObs, vecPred, pv.thres=0.0001, ansari.sign=0.05, thresMin=0.05,
                         thresMax=0.5,verbose=1,scale=TRUE){

# Check if supplied thresholds are valid
if(thresMin>thresMax){cat("Error, thresMax should be equal to or larger than thresMin\n");return()}


# Initializing threshold and threshold vector for keeping track of thresholds
thresAbs=thresMin
sq<-numeric()

#initializing threshold index (threshold count)
j=0

#initializing ansari p-values to keep track of ansari p-values for each threshold in sq
ansari=numeric()

# Initialize levels count
lv=numeric()

# Set backtracking flag. Start with flag=0 indicating significance not yet reached, backtracking not begun
flag=0

# If thresMin=thresMax, fixed threshold is used and we set flag=2, only one run of the algoritm with initial thresMin
if(thresMin==thresMax){flag=2}

# Evaluate optimum steps for algorithm
else {
 l.step <- signif((thresMax-thresMin)/10,1)
 s.step <- signif((thresMax-thresMin)/200,1)
}

while (1){

  # Print current threshold if verbose is 1 or larger
  if(verbose>=1){cat("\nCurrent thresAbs: ",thresAbs,"\n")}

  j=j+1

  # Save current threshold
  sq[j]<-thresAbs

  # temporary predicted values (to be updated)
  vecPredNow=vecPred

  #unmissing unique segment medians
  mnNow=unique(vecPred)
  mnNow=mnNow[!is.na(mnNow)]

  #continuing indicator otherwise get out of the loop
  cont=0

  while(cont==0 & length(mnNow)>1) {

        mnNow=sort(mnNow)  #currennt sorted vector of means
        n <- length(mnNow)  # number of means in mnNow

        # Print current number of levels (n) if verbose is 2 or larger
        if(verbose>=2){ cat("\r",n,":",length(unique(vecPred)),"\t")}

        # Get distances translated to copy number differences
        # Only distances to closest levels
        if(scale){d<-(2*2^mnNow)[-n]-(2*2^mnNow)[-1]}
        else{d<-(mnNow)[-n]-(mnNow)[-1]}

        #order distance between means with the closest on top and corresponding indices
        dst<-cbind(abs(d)[order(abs(d))],(2:n)[order(abs(d))],(1:(n-1))[order(abs(d))])

        #for each pair of means
        for (i in 1:nrow(dst))  {
                #set continuity index to "NOT continue" (=1)
                cont=1
                #test for combining of the two segment means
                out=combine.funcB(diff=dst[i,1],vecObs, vecPredNow, mnNow, mn1=mnNow[dst[i,2]], mn2=mnNow[dst[i,3]], pv.thres=pv.thres, thresAbs=if(scale){2*2^thresAbs-2}else{thresAbs})
                #if combine?
                if (out$pv > pv.thres) {

                       #set continuity index to "YES" (=0) and break out of the current pairs loop
                       cont=0

                       #update predicted values and segments
                       vecPredNow=out$vecPredNow
                       mnNow=out$mnNow
                       break
                 }                
          }               
 }

### When done merging for a given threshold, test for significance ####
        ansari[j]=my.ansari.test(sort(vecObs-vecPredNow), sort(vecObs-vecPred))$p.value
  if(is.na(ansari[j])){ansari[j]=0} # If too many numbers for test to be performed, a 0 is returned, resulting in no merging (please use fixed threshold to get any merging)
  lv[j]=length(mnNow) # get number of levels

### If backtracking flag=2, the merging is stopped at this thresMax (or fixed threshold) ###
  if(flag==2){ break }

  # If p.value is less than the significance threshold, set backtracking flag=1 (backtracking on)
  if(ansari[j]<ansari.sign){
                        flag=1
  }

        
### If backtracking is on, a smaller threshold is attempted ####
        if (flag){

        # Stop if backtracking is on and p.value is higher than sign threshold or threshold is less or equal to thresMin
        if (ansari[j]>ansari.sign | thresAbs == thresMin){

#        # Don't merge at all if all tested threshold including thresMin is significant
#                         if (ansari[j] <= ansari.sign) {
#                                 vecPredNow=vecPred
#                                 mnNow=unique(vecPred)
#                                 mnNow=mnNow[!is.na(mnNow)]
#                         }
                                                  
        break
        }

      # Attempt smaller threshold
        else {
        thresAbs=signif(thresAbs-s.step,3)

        # Set threshold to thresMin as a minimum
        if (thresAbs <= thresMin){ thresAbs = thresMin }
      }
        }
        

### Increase threshold if backtracking is not on ###
        else {thresAbs=thresAbs+l.step}

#### Control step so function won't keep running, max threshold = thresMax and if sign not reached, threshold = thresMax ###
          if (thresAbs >= thresMax){
        thresAbs=thresMax
                    flag=2
          }

} # End while


# Return list of results
return(list(vecMerged=vecPredNow,mnNow=mnNow,sq=sq,ansari=ansari))
}





my.ansari.test <- function (x, ...) {
UseMethod("ansari.test")
}

my.ansari.test.default <-
function (x, y, alternative = c("two.sided", "less", "greater"),
exact = NULL, conf.int = FALSE, conf.level = 0.95, ...)
{

  ##R_pansari and R_qansari are in stats namespace
  thisR_pansari <- stats:::R_pansari
  thisR_qansari <- stats:::R_qansari
  
 alternative <- match.arg(alternative)
 if (conf.int) {
   if (!((length(conf.level) == 1) && is.finite(conf.level) &&
     (conf.level > 0) && (conf.level < 1)))
     stop("'conf.level' must be a single number between 0 and 1")
 }
 DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
 x <- x[complete.cases(x)]
 y <- y[complete.cases(y)]
 m <- length(x)
 if (m < 1)
   stop("not enough 'x' observations")
 n <- length(y)
 if (n < 1)
   stop("not enough 'y' observations")
 N <- m + n
 r <- rank(c(x, y))
 STATISTIC <- sum(pmin(r, N - r + 1)[seq_along(x)])
 TIES <- (length(r) != length(unique(r)))
 if (is.null(exact))
   exact <- ((m < 50) && (n < 50))
 if (exact && !TIES) {
   pansari <- function(q, m, n) {
     .C(thisR_pansari, as.integer(length(q)), p = as.double(q),
       as.integer(m), as.integer(n))$p
   }
   PVAL <- switch(alternative, two.sided = {
     if (STATISTIC > ((m + 1)^2%/%4 + ((m * n)%/%2)/2))
       p <- 1 - pansari(STATISTIC - 1, m, n)
     else p <- pansari(STATISTIC, m, n)
     min(2 * p, 1)
   }, less = 1 - pansari(STATISTIC - 1, m, n), greater =
pansari(STATISTIC,
     m, n))
   if (conf.int) {
     qansari <- function(p, m, n) {
       .C(thisR_qansari, as.integer(length(p)), q = as.double(p),
        as.integer(m), as.integer(n))$q
     }
     alpha <- 1 - conf.level
     x <- sort(x)
     y <- sort(y)
     ab <- function(sig) {
       rab <- rank(c(x/sig, y))
       sum(pmin(rab, N - rab + 1)[seq_along(x)])
     }
     ratio <- outer(x, y, "/")
     aratio <- ratio[ratio >= 0]
     sigma <- sort(aratio)
     cci <- function(alpha) {
       u <- absigma - qansari(alpha/2, m, n)
       l <- absigma - qansari(1 - alpha/2, m, n)
       uci <- NULL
       lci <- NULL
       if (length(u[u >= 0]) == 0 || length(l[l > 0]) ==
        0) {
        warning("samples differ in location: cannot compute
confidence set, returning NA")
        return(c(NA, NA))
       }
       if (is.null(uci)) {
        u[u < 0] <- NA
        uci <- min(sigma[which(u == min(u, na.rm = TRUE))])
       }
       if (is.null(lci)) {
        l[l <= 0] <- NA
        lci <- max(sigma[which(l == min(l, na.rm = TRUE))])
       }
       if (uci > lci) {
        l <- absigma - qansari(alpha/2, m, n)
        u <- absigma - qansari(1 - alpha/2, m, n)
        u[u < 0] <- NA
        uci <- min(sigma[which(u == min(u, na.rm = TRUE))])
        l[l <= 0] <- NA
        lci <- max(sigma[which(l == min(l, na.rm = TRUE))])
       }
       c(uci, lci)
     }
     cint <- if (length(sigma) < 1) {
       warning("cannot compute confidence set, returning NA")
       c(NA, NA)
     }
     else {
       absigma <- sapply(sigma + c(diff(sigma)/2,
sigma[length(sigma)] *
        1.01), ab)
       switch(alternative, two.sided = {
        cci(alpha)
       }, greater = {
        c(cci(alpha * 2)[1], Inf)
       }, less = {
        c(0, cci(alpha * 2)[2])
       })
     }
     attr(cint, "conf.level") <- conf.level
     u <- absigma - qansari(0.5, m, n)
     sgr <- sigma[u <= 0]
     if (length(sgr) == 0)
       sgr <- NA
     else sgr <- max(sgr)
     sle <- sigma[u > 0]
     if (length(sle) == 0)
       sle <- NA
     else sle <- min(sle)
     ESTIMATE <- mean(c(sle, sgr))
   }
 }
 else {
   EVEN <- ((N%%2) == 0)
   normalize <- function(s, r, TIES, m = length(x), n = length(y)) {
##########################################
## Here is the problem: length(x) returns an integer
m <- as.double(m)
n <- as.double(n)
##########################################
     z <- if (EVEN)
       s - m * (N + 2)/4
     else s - m * (N + 1)^2/(4 * N)
     if (!TIES) {
       SIGMA <- if (EVEN)
        sqrt((m * n * (N + 2) * (N - 2))/(48 * (N -
         1)))
       else sqrt((m * n * (N + 1) * (3 + N^2))/(48 *
        N^2))
     }
     else {
       r <- rle(sort(pmin(r, N - r + 1)))
       SIGMA <- if (EVEN)
        sqrt(m * n * (16 * sum(r$lengths * r$values^2) -
         N * (N + 2)^2)/(16 * N * (N - 1)))
       else sqrt(m * n * (16 * N * sum(r$lengths * r$values^2) -
        (N + 1)^4)/(16 * N^2 * (N - 1)))
     }
     z/SIGMA
   }
   p <- pnorm(normalize(STATISTIC, r, TIES))
   PVAL <- switch(alternative, two.sided = 2 * min(p, 1 -
     p), less = 1 - p, greater = p)
   if (conf.int && !exact) {
     alpha <- 1 - conf.level
     ab2 <- function(sig, zq) {
       r <- rank(c(x/sig, y))
       s <- sum(pmin(r, N - r + 1)[seq_along(x)])
       TIES <- (length(r) != length(unique(r)))
       normalize(s, r, TIES, length(x), length(y)) -
        zq
     }
     srangepos <- NULL
     srangeneg <- NULL
     if (length(x[x > 0]) && length(y[y > 0]))
       srangepos <- c(min(x[x > 0], na.rm = TRUE)/max(y[y >
        0], na.rm = TRUE), max(x[x > 0], na.rm = TRUE)/min(y[y >
        0], na.rm = TRUE))
     if (length(x[x <= 0]) && length(y[y < 0]))
       srangeneg <- c(min(x[x <= 0], na.rm = TRUE)/max(y[y <
        0], na.rm = TRUE), max(x[x <= 0], na.rm = TRUE)/min(y[y <
        0], na.rm = TRUE))
     if (any(is.infinite(c(srangepos, srangeneg)))) {
       warning("cannot compute asymptotic confidence set or estimator")
       conf.int <- FALSE
     }
     else {
       ccia <- function(alpha) {
        statu <- ab2(srange[1], zq = qnorm(alpha/2))
        statl <- ab2(srange[2], zq = qnorm(alpha/2,
         lower.tail = FALSE))
        if (statu > 0 || statl < 0) {
         warning("samples differ in location: cannot compute confidence set, returning NA")
         return(c(NA, NA))
        }
        u <- uniroot(ab2, srange, tol = 1e-04, zq =
qnorm(alpha/2))$root
        l <- uniroot(ab2, srange, tol = 1e-04, zq = qnorm(alpha/2,
         lower.tail = FALSE))$root
        sort(c(u, l))
       }
       srange <- range(c(srangepos, srangeneg), na.rm = FALSE)
       cint <- switch(alternative, two.sided = {
        ccia(alpha)
       }, greater = {
        c(ccia(alpha * 2)[1], Inf)
       }, less = {
        c(0, ccia(alpha * 2)[2])
       })
       attr(cint, "conf.level") <- conf.level
       statu <- ab2(srange[1], zq = 0)
       statl <- ab2(srange[2], zq = 0)
       if (statu > 0 || statl < 0) {
        ESTIMATE <- NA
        warning("cannot compute estimate, returning NA")
       }
       else ESTIMATE <- uniroot(ab2, srange, tol = 1e-04,
        zq = 0)$root
     }
   }
   if (exact && TIES) {
     warning("cannot compute exact p-value with ties")
     if (conf.int)
       warning("cannot compute exact confidence intervals with ties")
   }
 }
 names(STATISTIC) <- "AB"
 RVAL <- list(statistic = STATISTIC, p.value = PVAL, null.value =
              c(`ratio of scales` = 1),
   alternative = alternative, method = "Ansari-Bradley test",
   data.name = DNAME)
 if (conf.int)
   RVAL <- c(RVAL, list(conf.int = cint, estimate = c(`ratio of scales` = ESTIMATE)))
 class(RVAL) <- "htest"
 return(RVAL)
}









#################################


combine.funcB <- function(diff,vecObs, vecPredNow, mnNow, mn1, mn2, pv.thres=0.0001, thresAbs=0)
{ 
  #observed values in the first segment
        vec1=vecObs[which(vecPredNow==mn1)]
  #observed values in the second segment
        vec2=vecObs[which(vecPredNow==mn2)]
        
  #if difference between segment medians does not exceed thresAbs, then set pv=1
        if (diff<=thresAbs) {
                pv=1
        }
  #otherwise test for difference in mean based on observed values
        else {
                if((length(vec1) > 10 & length(vec2) > 10) | sum(length(vec1),length(vec2))>100){
                        pv=wilcox.test(vec1,vec2)$p.value
                }
                else{pv=wilcox.test(vec1,vec2,exact=TRUE)$p.value  }       #/10^max(mn1,mn2)
                if(length(vec1) <= 3 | length(vec2) <= 3){pv=0}         
        }
        index.merged<-numeric()
  #if p-value exceeds pv.thres
        if (pv > pv.thres)      {
    #combine observed values
                vec=c(vec1,vec2)
    # Index values to be updated
                index.merged=which((vecPredNow==mn1) | (vecPredNow==mn2))               
    #update predicted values by median of the observed values
                vecPredNow[index.merged]=median(vec, na.rm=TRUE)
    #update segment medians  median of the observed values and remove one of the duplicates
                mnNow[which((mnNow==mn1) | (mnNow==mn2))]=median(vec, na.rm=TRUE)
                mnNow=unique(mnNow)
        }
        list(mnNow=mnNow, vecPredNow=vecPredNow, pv=pv)
}

#########################################



## in find.param.two.R
## Error in matrix(c(1 - p1, p1, p2, 1 - p2), ncol = 2, b = T) : 
##         T used instead of TRUE
## Error in matrix(c(1 - p1, p1, p2, 1 - p2), ncol = 2, b = T) : 
##         recursive default argument reference

## viterbvis
## Error in matrix(nrow = K, ncol = n, b = T) : 
##         T used instead of TRUE
## Error in matrix(nrow = K, ncol = n, b = T) : 
##         recursive default argument reference


## Error in max(-1e+05, log(BFGS.output$prior[i]), na.rm = T) :



















### Utility functions

add.names.as.attr <- function(x, anames) {
    if(is.null(anames))
        anames <- paste("sample.", 1:length(x[[1]]), sep = "")
    
    for(i in 1:length(x[[1]])) {
        attributes(x[[1]][[i]]) <- c(attributes(x[[1]][[i]]),
                                     "ArrayName" = anames[i])
    }
    return(x)
}


warn.too.few.in.chrom <- function(x, min.num.chrom = 20) {
  ## if too few samples, somethin funny is going on
  tt <- table(x)
  if(any(tt < min.num.chrom))
    warning("There are fewer than ", min.num.chrom, " observations",
            "in some group of ", deparse(substitute(x)),
            "!!!!!!!!!!. \n This is unlikely to make sense.\n",
            "Note that you can get weird errors, or no output at all,", 
            "if you are running parallelized with certaind methods.")
}

stop.na.inf <- function(x) {
  ## The code for many functions allows for dealing with NA and Inf, but
  ## would need to adjust other functions (as different arrays would have
  ## different length of pos, genenames, etc. So for now stop
  if(any(is.na(x)) | any(is.infinite(x)))
        stop("Either an NA or an infinite in the data: ",
             deparse(substitute(x)), ".\n",
             "   Eliminate those values or use imputation")
}


### Example of usage. Suppose we create missing values

## cghE1[1:10, 5:7] <- NA
## imputed.x <- my.impute.lowess(cghE1[1:40, 5:7], rep(1, 40))


my.impute.lowess <- function (x,
                              chrom.numeric,
                              Clone = NULL,
                              Pos = NULL,
                              chrominfo = human.chrom.info.Jul03,
                              maxChrom = 23,
                              smooth = 0.1)
{
  ## BEWARE: Pos MUST be in kilobases!!!
  if(is.null(Clone)) Clone <- 1:length(chrom.numeric)
  if(is.null(Pos)) Pos <- Clone
  aCGH.obj <- create.aCGH(data.frame(x),
                          data.frame(Clone = Clone,
                                     Chrom = chrom.numeric,
                                     kb = Pos))
  
    data.imp <- log2.ratios <- log2.ratios(aCGH.obj)
    clones.info <- clones.info(aCGH.obj)
    uniq.chrom <- unique(clones.info$Chrom)
    for (j in uniq.chrom[uniq.chrom <= maxChrom]) {
        cat("Processing chromosome ", j, "\n")
        centr <- chrominfo$centromere[j]
        indl <- which(clones.info$Chrom == j & clones.info$kb <
            centr)
        indr <- which(clones.info$Chrom == j & clones.info$kb >
            centr)
        kbl <- clones.info$kb[indl]
        kbr <- clones.info$kb[indr]
        for (i in 1:ncol(log2.ratios)) {
            if (length(indl) > 0) {
                vecl <- log2.ratios[indl, i]
                ind <- which(!is.na(vecl))
                if (length(ind) > 1)
                  data.imp[indl, i][-ind] <- approx(lowess(kbl[ind],
                    vecl[ind], f = smooth), xout = kbl[-ind])$y
            }
            if (length(indr) > 0) {
                vecr <- log2.ratios[indr, i]
                ind <- which(!is.na(vecr))
                if (length(ind) > 0)
                  data.imp[indr, i][-ind] <- approx(lowess(kbr[ind],
                    vecr[ind], f = smooth), xout = kbr[-ind])$y
            }
        }
    }
    prop.miss <- apply(data.imp, 2, prop.na)
    if (max(prop.miss, na.rm = TRUE) > 0) {
        for (i in 1:ncol(data.imp)) {
            vec <- data.imp[, i]
            ind <- which(is.na(vec))
            if (length(ind) > 0) {
                vec[ind] <- sapply(ind, function(i) {
                  chr <- clones.info$Chrom[i]
                  kb <- clones.info$kb[i]
                  if (kb >= chrominfo$centromere[chr])
                    median(vec[clones.info$Chrom == chr & clones.info$kb >=
                      chrominfo$centromere[chr]], na.rm = TRUE)
                  else median(vec[clones.info$Chrom == chr &
                    clones.info$kb < chrominfo$centromere[chr]],
                    na.rm = TRUE)
                })
                vec[is.na(vec)] <- 0
                data.imp[, i] <- vec
            }
        }
    }
    prop.miss <- apply(data.imp, 2, prop.na)
    if (max(prop.miss) > 0)
        print(paste("Missing values still remain in samples ",
            which(prop.miss > 0)))
    data.imp
}


the.time.with.ms <- function() {
    uu <- as.POSIXlt(Sys.time())
    return(paste(uu$hour, uu$min,
                 paste(unlist(strsplit(as.character(uu$sec), "\\.")),
                       collapse = ""), sep = ""))
}

tempdir2 <- function() {
    direxists <- TRUE
    while(direxists) {
        p1 <-  paste(round(runif(1, 1, 9999)),
                     the.time.with.ms(), sep = "_")
        p1 <- paste(tempfile(pattern = "tmpdir_ADaCGH_",
                             tmpdir = "."),
                    p1, sep = "_")
        if(!file.exists(p1)) direxists <- FALSE
    }
    dir.create(p1)
    return(p1)
}



### Simple examples of stochasticity of methods and working of SegmentPlotWrite

### library(ADaCGH)
### data(cghE1)
### tmpchr <- sub("chr", "", cghE1$Chromosome)
### chrom.numeric <- as.numeric(as.character(tmpchr))
### chrom.numeric[tmpchr == "X"] <- 23
### chrom.numeric[tmpchr == "Y"] <- 24
### rm(tmpchr)
### ### we need the data ordered
### reorder <- order(chrom.numeric,
###                  cghE1$UG.Start,
###                  cghE1$UG.End,
###                  cghE1$Name)
### cghE1 <- cghE1[reorder, ]
### chrom.numeric <- chrom.numeric[reorder]



### #### DNAcopy: differences between pSegmentDNA and
### ##   SegmentPlotWrite: because DNAcopy is variable
### set.seed(1)
### cna.obj <- CNA(as.matrix(cghE1[, 5:7]),
###               chrom = chrom.numeric,
###               maploc = cghE1$UG.Start,
###               data.type = "logratio")
### smoothed <- smooth.CNA(cna.obj)
### segmented1 <- segment(smoothed, undo.splits = "none", nperm = 10000)
### segmented2 <- segment(smoothed, undo.splits = "none", nperm = 10000)
### segmented3 <- segment(smoothed, undo.splits = "none", nperm = 10000)
### segmented3$output[32:34, ]
### segmented1$output[32:34, ]


### #### HMM: this is highly variable

### Clone <- 1:nrow(cghE1)
### obj <- create.aCGH(data.frame(cghE1[, 5:7]), 
###                    data.frame(Clone = Clone,
###                               Chrom = chrom.numeric,
###                               kb     = cghE1$UG.Start ))
### set.seed(1)
### res <- find.hmm.states(obj, aic = TRUE, bic = FALSE)
### res2 <- find.hmm.states(obj, aic = TRUE, bic = FALSE)

### res$states.hmm[[1]][1334:1342, 2 + 4]
### res2$states.hmm[[1]][1334:1342, 2 + 4]



### ### GLAD: no changes here between SegmentPlotWrite
### ### and pSegmentGLAD
### ## Note also that the testAdacghNum.py, with GALD, shows
### ## perfect 

### oGLAD <- pSegmentGLAD(cghE1[, 5:7],
###                     chrom.numeric)
### setwd("/tmp/mG")
### SegmentPlotWrite(cghE1[, 5:7], chrom.numeric,
###                  Pos = cghE1$UG.Start,
###                  mergeSegs = FALSE,
###                  idtype = "ug", organism = "Hs",
###                  method = "GLAD",
###                  geneNames = cghE1[, 1], commondata = cghE1[, 1:4])

### spw.out <- read.table("/tmp/mG/GLAD.output.txt")

### summary(oGLAD$segm$S1[, 2] - spw.out$S1.Smoothed)
### summary(oGLAD$segm$S2[, 2] - spw.out$S2.Smoothed)
### summary(oGLAD$segm$S3[, 2] - spw.out$S3.Smoothed)
### summary(oGLAD$segm$S1[, 3] - spw.out$S1.Status)
### summary(oGLAD$segm$S2[, 3] - spw.out$S2.Status)
### summary(oGLAD$segm$S3[, 3] - spw.out$S3.Status)


my.usr2png <- function(xy, imWidth, imHeight) {
    dev <- dev.cur()
    xy <- fig2dev(plt2fig(usr2plt(xy,dev),dev),dev)
    cbind(
          ceiling(xy[,1]*imWidth),
          ceiling((1-xy[,2])*imHeight)
          )
}





######################################

###### Diagnositc plots

#### For now, these are not documented nor publicly available

DNAcopyDiagnosticPlots <- function(CNA.object, CNA.smoothed.object,
                                   smoothed.CNA.object,
                                   array.chrom = FALSE, chrom.numeric = NULL) {
    numarrays <- ncol(CNA.object) - 2
    if(!array.chrom) {
        par(pty = "s")
        for(i in 1:numarrays)
            plot(CNA.object[, (i + 2)], smoothed.CNA.object[, (i + 2)],
                 main = colnames(CNA.object[i + 2]),
                 xlab = "Original data", ylab = "Smoothed data")
    } else {
        if(is.null(chrom.numeric))
            stop("To use array.chom = TRUE you need to provide the vector with chrom. number")
        par(pty = "s")
        for(i in 1:numarrays) { cat("\n Doing sample ", i, "\n")
        ##    par(mfrow = c(4, 6))
            ncr <- unique(chrom.numeric)
            for(j in ncr) {
                x <- CNA.object[chrom.numeric == j, (i + 2)]
                y <- smoothed.CNA.object[chrom.numeric == j, (i + 2)]
                plot(x, y,
                     main = paste(colnames(CNA.object[i+2]), "; Chr ", j,
                     sep = ""),
                     xlab = "Original data", ylab = "Smoothed data")
            }
        }
    }

#####  example
##### par(mfrow = c(1, 3))
##### DNAcopyDiagnosticPlots(CNA.object,
#####                        smoothed.CNA.object)

##### par(mfrow = c(3, 3))
##### par(ask = TRUE)
##### DNAcopyDiagnosticPlots(CNA.object,
#####                        smoothed.CNA.object, array.chrom = TRUE,
#####                        chrom.numeric = chrom.numeric)

}


WaveletsDiagnosticPlots <- function(acghdata, chrom.numeric) {
    numarrays <- ncol(acghdata)
    if(numarrays < 2) {
        stop("Only one array. No histogram possible")
    }
    cnu <- unique(chrom.numeric)
    ncrom <- length(cnu)
    dat <- as.matrix(acghdata)
    ar1s <- matrix(nrow = numarrays, ncol = ncrom)
    for(cn in 1:ncrom) { ## zz: parallelize this?
        index.dat <- which(chrom.numeric == cnu[cn])
        for(subject in 1:numarrays) {
            trythis <- try(
            ar1s[subject, cn] <-
                as.vector(acf(dat[index.dat, subject],
                              lag.max = 1, plot = FALSE)$acf)[2]
                           )
            if(inherits(trythis, "try-error"))
                caughtOurError(paste("acf bombed unexpectedly with error",
                                     trythis, ". \n Please let us know so we can fix the code."))

        }
    }
    rm(cn, subject, index.dat)
    ##         pdf("Autocorrelation.plots.pdf", width = 17.6, height = 12.5)
    ##        par(mfrow = c(6,4))
    ##        par(oma = c(2, 2, 2, 2))
    if (numarrays > 1) {
        for(i in 1:ncrom)
            hist(ar1s[, i], main = paste("Chr", i), xlab = "Autocorrelation, lag 1")
        mtext("Autocorrelation coefficient lag 1", side = 3,
              outer = TRUE, line = 0.5, font = 2)
        ##        dev.off()
    } else {
        plot(x = c(0, 1), y = c(0, 1),
             type = "n", axes = FALSE, xlab = "", ylab = "")
        box()
        text(0.5, 0.7, "Only one array.")
        text(0.5, 0.5,
             "No histogram possible.")
    }

##### example
##### data(cghMCRe)

##### chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
##### chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
##### chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

##### par(mfrow = c(5, 5))
##### WaveletsDiagnosticPlots(cghMCRe[, 5:7], chrom.numeric)
  }
