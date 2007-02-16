## Linking to Toronto DB for a chromosome
##http://projects.tcag.ca/variation/cgi-bin/tbrowse/tbrowse?source=hg18&table=Locus&show=table&keyword=&flop=AND&fcol=_C19&fcomp==&rnum=0&fkwd=chr13&cols=



## Analysis: must return observed, smoothed, state.
## Plot: must take the above plus geneNames and Position (where that can
##       be a vector of consecutive integers)



## .__ADaCGH_WEB_APPL <- TRUE in web appl!

if(exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) {
    warningsForUsers <- vector()
} else {
    warningsForUsers <- warning
}


## where do we live? to call the python script
.calltoMap.py <- function() {
    tmp <- library()[[2]]
    pathpy <- tmp[which(tmp[, 1] == "ADaCGH"), 2]
    .python.call <- paste(pathpy, "/ADaCGH/Python/toMap.py", sep = "")
    return(.python.call)
}

.python.toMap.py <- .calltoMap.py()
    
##############################################


###  Visible stuff

mpiInit <- function(wdir = getwd(), minUniverseSize = 15) {
    library(Rmpi)
    if(mpi.universe.size() < minUniverseSize) {
        stop("MPI problem: universe size < minUniverseSize")
    }
##    mpi.spawn.Rslaves(nslaves= mpi.universe.size())
    mpi.spawn.Rslaves(hosts = sample(lamhosts()))
    ## mpi.setup.rngstream() ## or 
    mpi.setup.sprng()
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    library(papply)
    mpi.remote.exec(library(ADaCGH))
    mpi.bcast.Robj2slave(wdir)
    mpi.remote.exec(setwd(wdir))    
}

pSegmentACE <- function(x, chrom.numeric, ...) {
    ## ACE is already parallelized
    ACE(x, chrom.numeric, echo = FALSE, coefs = file.aux, Sdev = 0.2)
}


pSegmentHMM <- function(x, chrom.numeric, ...) {
    out <- papply(data.frame(x),
                  function(z) hmmWrapper(z, Chrom = slave_chrom),
                  papply_commondata = list(
                  slave_chrom = chrom.numeric))
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    class(outl) <- c("adacgh.generic.out","mergedHMM")
    return(outl)
}



pSegmentGLAD <- function(x, chrom.numeric, ...) {
    out <- papply(data.frame(x),
                  function(z) gladWrapper(z,
                                          Chrom = slave_chrom),
                  papply_commondata = list(
                  slave_chrom = chrom.numeric))
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    class(outl) <- c("adacgh.generic.out", "adacghGLAD")
    return(outl)
}    

pSegmentBioHMM <- function(x, chrom.numeric, Pos, ...) {
    out <- papply(data.frame(x),
                  function(z) BioHMMWrapper(z,
                                         Chrom = slave_chrom,
                                         Pos    = slave_kb),
                  papply_commondata = list(
                  slave_chrom = chrom.numeric,
                  slave_kb    = Pos))

    te <- unlist(lapply(out, function(x) inherits(x, "try-error")))
    if(any(te)) {
        m1 <- "The BioHMM code occassionally crashes (don't blame us!)."
        m2 <- "You can try rerunning it a few times."
        m3 <- "You can also tell the original authors that you get the error "
        mm <- paste(m1, m2, m3, out[[which(te)]])
        caughtError(mm)
    }
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    outl$pos <- Pos
    class(outl) <- c("adacgh.generic.out","mergedBioHMM")
    return(outl)
}


pSegmentCGHseg <- function(x, chrom.numeric, ...) {
    ## Beware: we parallelize over chromosoes AND subjects
 
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
      
    out0 <- papply(datalist, function(z) CGHsegWrapper(z))
  
    out <- list()
    klist <- 1
    for(i in 1:nsample) {
        out[[i]] <- out0[[klist]]
        for(j in nchrom[-1]) {
            klist <- klist + 1
            out[[i]] <- rbind(out[[i]], out0[[klist]])
        }
    }
    outl <- list()
    outl$segm <- out
    outl$chrom.numeric <- chrom.numeric
    class(outl) <- c("adacgh.generic.out", "CGHseg")
    return(outl)
}


pSegmentPSW <- function(x, chrom.numeric, sign = -1,
                        nIter = 1000, prec = 100,  p.crit = 0.10,
                        name = NULL, ...) {
    numarrays <- ncol(x)
    ncrom <- length(unique(chrom.numeric))
    out <- list()
    ## out$Data
    out$plotData <- list()
    if (exists(".__ADaCGH_WEB_APPL", env = .GlobalEnv)) { ## send to PaLS
      print("testing value of .AD...")
      print(paste(".__ADaCGH_WEB_APPL  is", .__ADaCGH_WEB_APPL))
        palsVect <- vector()
        palsL <- list()
    }
 
    datalist <- list()
    for (i in 1:ncol(x)) {
        datalist[[i]] <- list()
        datalist[[i]]$data <- x[, i]
        datalist[[i]]$num  <- i
    }
    papout <- papply(datalist,
                  function(z) my.sw3b(z$data,
                                      chrom = sl.chrom.numeric,
                                      sign = sl.sign,
                                      p.crit = sl.p.crit,
                                      nIter = sl.nIter,
                                      prec = sl.prec,
                                      name = paste(sl.name, colnamesdata[z$num], sep = ""),
                                      highest = FALSE),
                  papply_commondata = list(
                  sl.chrom.numeric = chrom.numeric,
                  sl.sign = sign,
                  sl.nIter = nIter,
                  sl.prec = prec,
                  sl.name = name,
                  sl.p.crit = p.crit,
                  colnamesdata = colnames(x)))
    if(any(lapply(papout, function(z) z == "swt.perm.try-error"))) {
        m1 <- "There was a problem in the PSW routine; this is "
        m2 <- "probably related to the global thresholding + within "
        m3 <- "chromosome perm test with your data."
        m4 <- "You might want to try another method, or the original "
        m5 <- " (thresholding within chromosome) PSW "
        mm <- c(m1, m2, m3, m4, m5)
        caughtOtherError(mm)
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
    return(out)
}


pSegmentWavelets <- function(x, chrom.numeric, mergeSegs = TRUE,
                             minDiff = 0.25,
                             minMergeDiff = 0.05,
                             thrLvl = 3, initClusterLevels = 10, ...) {
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
    
    funwv <- function(ratio) {
        wc   <- modwt(ratio, "haar", n.levels=thrLvl)
        thH  <- our.hybrid(wc, max.level=thrLvl, hard=FALSE)
        recH <- imodwt(thH)
        ## cluster levels
        pred.ij <- segmentW(ratio, recH, minDiff=minDiff,
                            n.levels = initClusterLevels)
        labs <- as.character(1:length(unique(pred.ij)))
        state <- as.integer(factor(pred.ij, labels=labs))
        return(cbind(Observed = ratio,
                     Smoothed = pred.ij,
                     State = state))
    }
    out0 <- papply(datalist, funwv,
                   papply_commondata =list(thrLvl = thrLvl,
                   minDiff = if(mergeSegs) minMergeDiff else minDiff))

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
        class(outl) <- c(class(out), "CGH.wave", "CGH.wave.merged",
                         "adacgh.generic.out")
        return(outl)
    }
        
}


pSegmentDNAcopy <- function(x, chrom.numeric, mergeSegs = TRUE, smooth = TRUE,
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
       
    papfunc <- function(data) {
        if(slave_smooth)
            data <- internalSmoothCNA(data,
                                      chrom.numeric = slave_cnum,
                                      smooth.region = 2, outlier.SD.scale = 4,
                                      smooth.SD.scale = 2, trim = 0.025)
        outseg <-
            internalDNAcopy(data,
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
    class(outl) <- "DNAcopy"
    if(mergeSegs) class(outl) <- c(class(outl), "adacgh.generic.out")
    return(outl)
}






segmentPlot <- function(x, geneNames,
                        chrom.numeric = NULL,
                        cghdata = NULL,
                        arraynames = NULL,
                        idtype = "ug",
                        organism = "Hs",
                        yminmax = NULL,
                        numarrays = NULL,
                        ...) {
    if(is.null(numarrays)) {
        if(!is.null(arraynames)) numarrays <- length(arraynames)
        if(!is.null(cghdata)) numarrays <- ncol(cghdata)
    }
    if(is.null(yminmax)) {
        if(is.null(cghdata))
            stop("At least one of yminmax or cghdata has to be specified")
        yminmax <- c(min(as.matrix(cghdata)),
                     max(as.matrix(cghdata)))
    }
    
    if(is.null(arraynames)) arraynames <- colnames(cghdata)
    if(is.null(arraynames)) arraynames <- paste("sample.", 1:numarrays, sep = "")

    if(inherits(x, c("adacgh.generic.out"))) {
        geneLoc <- if(inherits(x, "mergedBioHMM")) x$pos else NULL
        if(inherits(x, "CGH.ACE.summary")) {
            original.pos <- 2
            segment.pos  <- NULL
        } else {
            original.pos <- 1
            segment.pos <- 2
        }
        if (inherits(x, "CGHseg") |
            (inherits(x, "CGH.wave") & (!inherits(x, "CGH.wave.merged")))){
            colors <- c(rep("orange", 3), "blue")
        } else {
            colors <- c("orange", "red", "green", "blue")
        }
        tmp_papout <-
            papply(as.list(1:numarrays),
                   function(z) {
                       cat("\n Doing sample ", z, "\n")
                       plot.adacgh.nonsuperimpose(res = res,
                                                  chrom = cnum_slave,
                                                  arraynum = z,
                                                  main = arraynames[z],
                                                  colors = colors,
                                                  ylim = yminmax,
                                                  geneNames = geneNames,
                                                  idtype = idtype,
                                                  organism = organism,
                                                  geneLoc = pos_slave)
                   },
                   papply_commondata =
                   list(res = x$segm,
                        cnum_slave= x$chrom.numeric,
                        arraynames = arraynames,
                        geneNames = geneNames,
                        idtype = idtype,
                        organism = organism,
                        colors = colors,
                        pos_slave = geneLoc))
        plot.adacgh.superimp(x$segm, x$chrom.numeric,  geneNames = geneNames,
                             main = "All_arrays",
                             colors = colors,
                             ylim= yminmax,
                             idtype = idtype, organism = organism,
                             geneLoc = geneLoc)
    }  else if(inherits(x, "CGH.PSW")) {
        if(x$plotData[[1]]$sign < 0) {
            main <- "Losses."
        } else {
            main <- "Gains."
        }
        tmp_papout <-
            papply(as.list(1:numarrays),
                   function(z) {
                       cat("\n Doing sample ", z, "\n")
                       sw.plot3(logratio = data_slave[[z]]$logratio,
                                         sign = data_slave[[z]]$sign,
                                         swt.perm = data_slave[[z]]$swt.perm,
                                         rob = data_slave[[z]]$rob,
                                         swt.run = data_slave[[z]]$swt.run,
                                         p.crit = data_slave[[z]]$p.crit.bonferroni,
                                         chrom = data_slave[[z]]$chrom,
                                         main = paste(main_slave, arraynames[z], sep=""),
                                         geneNames = geneNames,
                                         idtype = idtype,
                                         organism = organism)
                   },
                   papply_commondata = list(
                   data_slave = x$plotData,
                   main_slave = main,
                   arraynames = arraynames,
                   geneNames = geneNames,
                   idtype = idtype,
                   organism = organism))
    } else {
        stop("No plotting for this class of objects")
    }

##     }  else if (inherits(x, "DNAcopy")) {
##         ## FIXME: this is really obsolete stuff
##         ## should not be used in the web-based app
##         ## we leave it commented out; should be functional
##         ## but we ain't using it anymore.
##         tmp_papout <- papply(as.list(1:numarrays),
##                              function(z) {
##                                  cat("\n Doing sample ", z, "\n")
##                                  plot.olshen2(res = res,
##                                               arraynum = z,
##                                               main = arraynames[z],
##                                               html = TRUE,
##                                               geneNames = geneNames,
##                                               idtype = idtype,
##                                               organism = organism)
##                              },
##                              papply_commondata = list(res = x,
##                              arraynames = arraynames,
##                              geneNames = geneNames,
##                              idtype = idtype,
##                              organism = organism))
##         plot.olshen3(x, geneNames = geneNames,
##                      main = "All_arrays", ylim = yminmax,
##                      html = html,
##                      arraynums = 1:numarrays,
##                      idtype = idtype, organism = organism) ## all genome
        
##         plot.olshen4(x,  geneNames = geneNames,
##                      main = "All_arrays", ylim = yminmax,
##                      html = html,
##                      arraynums = 1:numarrays,
##                      idtype = idtype, organism = organism) ## by chromosome
##         }
}



plateauPlot <- function(obj, ...) {
    UseMethod("plateauPlot")
}

plateauPlot.DNAcopy <- function(obj, ...) {
    ## For each sample:
    plot.DNAcopy2(obj, ...)
    ## superimposing all:
    plot.DNAcopy3(obj, pt.pch = "", pt.cex = 0, ...)
    ## One for all
    plot.DNAcopy4(obj, ...)
    }
plateauPlot.CGH.wave <- function(obj, cghdata, ...) {
    numarrays <- ncol(cghdata)
    yminmax <- c(min(as.matrix(cghdata)),
                 max(as.matrix(cghdata)))
    ## by array:
    plateau.wavelets(obj, cghdata, by.array = TRUE)
    ## arrays, superimposed
    plateau.wavelets(obj, cghdata, by.array = TRUE, superimpose = TRUE,
                     ylim = yminmax)
    ## all, collapsed
    plateau.wavelets(obj, cghdata, by.array = FALSE)
}






SegmentPlotWrite <- function(data, chrom,
                             mergeSegs, Pos,
                             idtype, organism,
                             method,
                             geneNames,
                             commondata, ...) {
    ymax <- max(data)
    ymin <- min(data)
    numarrays <- ncol(data)
    
    fseg <- get(paste("pSegment", method, sep = ""))
    trythis <- try(
                   segmres <- fseg(data, chrom,
                                   Pos = Pos,
                                   mergeSegs = mergeSegs, ...)
                   )
    if(inherits(trythis, "try-error"))
        caughtOurError(trythis)
    cat("\n\n Segmentation done \n\n")

    save.image()
    save(segmres, file = "segmres.RData")

    trythis <- try(doMCR(segmres$segm, chrom = chrom, data = data,
                         Pos = Pos, ...))
    if(inherits(trythis, "try-error"))
        caughtOurError(trythis)
    
    trythis <- try(segmentPlot(segmres,
                               geneNames = geneNames,
                               chrom.numeric = chrom,
                               cghdata = data,
                               idtype = idtype,
                               organism = organism))
    if(inherits(trythis, "try-error"))
        caughtOurError(trythis)
    cat("\n\n Plotting done \n\n")

    trythis <- try(writeResults(segmres,
                                data, commondata))
    if(inherits(trythis, "try-error"))
        caughtOurError(trythis)
}                                


######################################

###### Diagnositc plots

#### For now, these are not documented nor publicly available

DNAcopyDiagnosticPlots <- function(CNA.object, CNA.smoothed.object,
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
    ncrom <- length(unique(chrom.numeric))
    dat <- as.matrix(acghdata)
    ar1s <- matrix(nrow = numarrays, ncol = ncrom)
    for(cn in 1:ncrom) { ## zz: parallelize this?
        index.dat <- which(chrom.numeric == cn)
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

writeResults.CGH.PSW <- function(obj, acghdata, commondata, file = "PSW.output.txt", ...) {
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
    if(inherits(obj, "adacgh.generic.out")) {
        print.adacgh.generic.results(obj, acghdata,
                                     commondata, output = file)
    } else {
        print.olshen.results(obj, acghdata, commondata,
                             merged = merged, output = file)
    }
}

writeResults.CGHseg <- function(obj, acghdata, commondata, 
                                 file = "CGHseg.output.txt", ...) {
    print.adacgh.generic.results(obj, acghdata,
                                commondata, output = file,
                                 send_to_pals = FALSE)
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
##     out <- data.frame(ID = commondata$name,
##                       Chromosome = commondata$chromosome,
##                       Start = commondata$start,
##                       End = commondata$end,
##                       MidPoint = commondata$MidPoint)

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
                     function(z) commondata$name[which( z == -1)])
        writeForPaLS(Ids, colnames(xcenter), "Lost_for_PaLS.txt")
        
        Ids <- apply(out[, cols.look, drop = FALSE], 2,
                     function(z) commondata$name[which( z == 1)])
        writeForPaLS(Ids, colnames(xcenter), "Gained_for_PaLS.txt")

        Ids <- apply(out[, cols.look, drop = FALSE], 2,
                     function(z) commondata$name[which( z != 0)])
        writeForPaLS(Ids, colnames(xcenter), "Gained_or_Lost_for_PaLS.txt")
    }

}



print.adacgh.cghseg.results <- function(res, xcenter,
                                 commondata,
                                 output = "CGHseg.results.txt"){
    ## This function "stretches out" the output and creates a table
    ## that matches the original names, etc.
    out <- data.frame(commondata)
##     out <- data.frame(ID = commondata$name,
##                       Chromosome = commondata$chromosome,
##                       Start = commondata$start,
##                       End = commondata$end,
##                       MidPoint = commondata$MidPoint)

    for(i in 1:ncol(xcenter)) {
            out <- cbind(out, res$segm[[i]][, c(1, 2)])
    }
    
    colnames(out)[(ncol(commondata) + 1):(ncol(out))] <-
        paste(rep(colnames(xcenter),rep(2, ncol(xcenter))),
              c(".Original", ".Smoothed"),
              sep = "")

    write.table(out, file = output,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
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

hmmWrapper <- function(logratio, Chrom, Pos = NULL) {
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


BioHMMWrapper <- function(logratio, Chrom, Pos) {
  logratio <- matrix(logratio, ncol=1)
  uchrom <- unique(Chrom)
##  obs <- NULL
  smoothed <- NULL
  for (ic in uchrom) {
      ydat <- logratio[Chrom == ic]
      n <- length(ydat)
      res <- try(fit.model(sample = 1, chrom = ic, dat = matrix(ydat, ncol = 1),
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
    segmentus2 <-
        mergeLevels(vecObs  = observed,
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


piccardsK <- function(loglik, n, s = -0.5) {
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

piccardsStretch <- function(obj, k, n, logratio) {
    if(k > 1) {
        poss <- obj@breakpoints[[k]]
        start <- c(1, poss)
        end <- c(poss - 1, n)
        smoothed <- mapply(function(start, end) mean(logratio[start: end]), start, end)
        reps <- diff(c(start, n + 1))
        smoothed <- rep(smoothed, reps)
        state <- rep(1:k, reps)
    } else {
        smoothed <- rep(mean(logratio), n)
        state <- rep(1, n)
    }
    return(cbind(smoothed = smoothed, state = state))
}
    
CGHsegWrapper <- function(logratio, maxseg = NULL, maxk = NULL) {
    ## beware, this is by chrom!!
    ## Just in case:
    n <- length(logratio)
    if(is.null(maxseg)) maxseg <- n/2
    if(is.null(maxk)) maxk   <- n
    obj1 <- tilingArray:::segment(logratio, maxseg = maxseg,
                                  maxk = maxk)
    optk <- piccardsK(obj1@logLik, n)
    finalsegm <- piccardsStretch(obj1, optk, n, logratio)
    return(cbind(Observed = logratio, SmoothedMean = finalsegm[, 1],
                 Alteration = finalsegm[, 2]))
}






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


internalSmoothCNA <- function(acghdata, chrom.numeric,
                              smooth.region = 2, outlier.SD.scale = 4,
                              smooth.SD.scale = 2, trim = 0.025) {
    ## this is just the original smoothCNA funct. adapted to use
    ## a single array and to be parallelized and fed to internalDNAcopy
   chrom <- chrom.numeric
   uchrom <- unique(chrom)
   genomdat <- acghdata
   ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
   trimmed.SD <- sqrt(trimmed.variance(genomdat[ina], trim))
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




internalDNAcopy <- function(acghdata,
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
    undo.splits <- "prune"
    genomdati <- acghdata
    ina <- which(!is.na(genomdati) & !(abs(genomdati)==Inf))

    ## The code allows for dealing with NA and Inf, but would need to
    ## adjust other functions (as different arrays would have different
    ## length of pos, genenames, etc. So for now stop:
    if (length(ina) != length(genomdati))
        stop("Either an NA or an infinite in the data")

    genomdati <- genomdati[ina]
    trimmed.SD <- sqrt(trimmed.variance(genomdati, trim))
    chromi <- chrom.numeric[ina]
    sample.lsegs <- NULL
    sample.segmeans <- NULL
    for (ic in uchrom) {
        segci <- changepoints(genomdati[chromi==ic], data.type = "logratio", alpha, 
                              sbdry, sbn, nperm, p.method, window.size, overlap, kmax,
                              nmin, trimmed.SD, undo.splits, undo.prune,
                              undo.SD, verbose = 2)
        sample.lsegs <- c(sample.lsegs, segci$lseg)
        sample.segmeans <- c(sample.segmeans, segci$segmeans)
    }
    if(length(sample.lsegs) != length(sample.segmeans))
        stop("Something terribly wrong: length(sample.lsegs) != length(sample.segmeans).")
    stretched.segmeans <- rep(sample.segmeans, sample.lsegs)
    stretched.state    <- rep(1:length(sample.lsegs), sample.lsegs)
    return(cbind(Observed = genomdati, Predicted = stretched.segmeans,
                 State = stretched.state))
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
    
    funwv <- function(ratio) {
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

    papout <- papply(datalist, funwv,
                     papply_commondata =list(thrLvl = thrLvl,
                     minDiff = minDiff))
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





my.sw3b <- function(logratio, chrom, sign = -1, p.crit = PSW.p.crit,
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

    papout <- lapply(swtlist, funsw)

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
                      idtype = idtype, organism = organism) {   
    ## this puts the chr call
    ## geneNames often = positions.merge1$name

    if(is.null(nameIm)) nameIm <- main
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
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
        f1 <- function(index) {
            si <- sign.segments[index]
            segments(swt.run$start[si] - 0.5, red.pos,
                     swt.run$start[si] - 0.5 + swt.run$length[si],
                     red.pos, col = "red", lwd = 2)
        }
        tmp <- mapply(f1, 1:length(sign.segments))
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

        f1 <- function(xleft, xright, nd)
            imRect(xleft, maxlr, xright, minlr - 10,
                   title = paste("Chromosome", nd),
                   alt = paste("Chromosome", nd),
                   href= paste("Chr", nd, "@", nameIm, ".html", sep =""))
        rectslist <- mapply(f1, xleft, xright, nd, SIMPLIFY=FALSE)
        for(ll in 1:length(rectslist))
            addRegion(im1) <- rectslist[[ll]]
        createIM2(im1, file = paste(nameIm, ".html", sep = ""))
        imClose(im1)
    }

    if(html) { ## here is chromosome specific code
        pixels.point <- 3
        chrheight <- 500
        for(cnum in 1:length(chrom.nums)) {
##            cat("\n .... doing chromosome ", cnum, ": ")
            indexchr <- which(chrom == chrom.nums[cnum])
            chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
            chrwidth <- max(chrwidth, 800)
            im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                             height = chrheight, width = chrwidth,
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
##             last.point.pixel.should.be <- chrwidth - ccircle[1, 1]
##             lichr <- length(indexchr)
##             ## recall ccircle is transposed
##             correction.factor <- (ccircle[1, lichr] - last.point.pixel.should.be)/lichr
##             richr <- (1:lichr) * correction.factor
##             ccircle[1, ] <- round(ccircle[1 ,] - richr)
            nameChrIm <- paste("Chr", chrom.nums[cnum], "@", nameIm, sep ="")
            write(ccircle, file = paste("pngCoordChr_", nameChrIm, sep = ""),
                  sep ="\t", ncolumns = 3)
            write(as.character(geneNames[indexchr]),
                  file = paste("geneNamesChr_", nameChrIm, sep = ""))
            imClose(im2)
            ## call the Python function
            system(paste(.python.toMap.py,  nameChrIm, 
                 idtype, organism, sep = " "))
            

        } ## looping over chromosomes
    } ## if html
    cat("\n")
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
  class(obj.ace) <- "ACE.analysis"

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
    
    v1 <- t(as.matrix(mapply(function(alpha1, beta1) {								                        z2-(alpha1+beta1*z1) }
                             , alpha1=alpha1, beta1=beta1)))
    v2 <- t(as.matrix(mapply(function(alpha2, beta2) {
        z2-(alpha2+beta2*z1) }
                             , alpha2=alpha2, beta2=beta2)))
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
    class(ace) <- "ACE.analysis"
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
			Sdev <- list(Sdev=sqrt(var(x[indexes])*(n-1)/n), n=n)
			}
		else {
			Sdev <- list(Sdev=0, n=n)
			}
		Sdev
	}
	else {
		Sdev <- 0
		Sdev
		}
}


ace.analysisP <-function(x) {
    ace.analysis.C(x, coefs, Sdev, array.names)
}




 ACE <- function(x, Chrom, coefs = file.aux, Sdev=0.2, echo=FALSE) {

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
 		first.estimate <- papply(genes, ace.analysisP,
                                          papply_commondata =list(coefs = coefs,
                                          Sdev = Sdev,
                                          echo = FALSE, array.names = array.names))
 		Sdevs.estimate <- papply(first.estimate, sd.ACE.analysis)
 		Sdevs <- unlist(lapply(Sdevs.estimate, "[", 1))
 		Sdevs <- mean(Sdevs[Sdevs>0])
                if(is.nan(Sdevs)) Sdevs <- 0
 		res <- papply(genes, ace.analysisP,
                               papply_commondata =list(coefs=coefs,
                               Sdev = Sdevs,
                               echo = FALSE, array.names = array.names))
 		class(res) <- "ACE"
 		}
 	else {
 		Sdevs <- matrix(NA, ncol(x),nchrom)
 		res <- list()
 		for(i in 1:ncol(x)) {
 			genes <- split(x[,i], Chrom)
                         first.estimate <- papply(genes, ace.analysisP,
                                                  papply_commondata =list(coefs= coefs,
                                                  Sdev = Sdev,
                                                  echo = FALSE, array.names = array.names[i]))
                         Sdevs.estimate <- papply(first.estimate, sd.ACE.analysis)
 			Sdevs[i,] <- unlist(lapply(Sdevs.estimate, "[", 1))
 		}
 		Sdevs <- mean(Sdevs[Sdevs>0])
                if(is.nan(Sdevs)) Sdevs <- 0
 		for (i in 1:ncol(x)) {
 			genes <- split(x[,i], Chrom)
                         res[[i]] <- papply(genes, ace.analysisP,
                                       papply_commondata =list(coefs= coefs,
                                       Sdev = Sdevs,
                                       echo = FALSE, array.names = array.names[i]))
 			class(res[[i]]) <- "ACE"
 		}
 		class(res) <- "ACE.array"
 	}
        res$chrom.numeric <- Chrom
 	invisible(res)
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
	gene.clusters <- mapply(function(x,y,z) {w<-rep(NA, length(x));if(length(y)>0) w[y]<-z;w}, x=obs, y=altered, z=gene.clusters)
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
        out$segm <- cbind(Observed = res$x, Smoothed = medians.state,
                          State = res$Gain.Loss)
        out$chrom.numeric <- chrom.numeric
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
    object$chrom.numeric <- NULL
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
        class(out) <- c("adacgh.generic.out", "summaryACE")
    attr(out, "aceFDR.for.output") <- aceFDR.for.output
    return(out)
        
        
    }




##################################################
##################################################
##################################################
##################################################


##### Internal stuff, mostly for the web-based part



PSWtoPaLS <- function(x = .__PSW_PALS.Lost_for_PaLS.txt,
                      y = .__PSW_PALS.Gained_for_PaLS.txt,
                      out = "Gained_or_Lost_for_PaLS.txt") {
    ## Creat the gained or lost
    ooo <- function(u, v) mapply(function(x, y, nx)
                                 return(c(paste("#", nx, sep = ""), x, y)),
                                 u, v, names(u))
    write(unlist(ooo(x, y)), file = out)     
}




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
    } else {
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
    } else {
        message <- paste("It looks like you found a bug. Please let us know. ", message)
        stop(message)
    }
}
    

my.html.data.frame <- function (object, first.col = "Name",
                             file = paste(first.word(deparse(substitute(object))), 
                             "html", sep = "."), append = FALSE, link = NULL, linkCol = 1, 
                             linkType = c("href", "name"), ...) 
{
    ## modifying html, from Hmisc: Their function always has first column
    ## named "Name". I allow to pass a name.
   
    linkType <- match.arg(linkType)
    x <- format.df(object, ...)
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


imClose <- function (im) {
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

linkGene2 <- function(id) {
    ## Returns a string for use in a web page with a call
    ## to IDClight.
    if ((idtype == "None") | (organism == "None"))
        return(id)
    else 
        paste("http://idclight.bioinfo.cnio.es/IDClight.prog",
              "?idtype=", idtype, "&id=", id, "&internal=0&org=",
              organism, sep = "")
}

linkGene <- function(id) {
    ## Returns a string for use in a web page with a call
    ## to IDClight.
    if ((idtype == "None") | (organism == "None"))
        return(id)
    else 
        paste("<a href=\"http://idclight.bioinfo.cnio.es/IDClight.prog",
              "?idtype=", idtype, "&id=", id, "&internal=0&org=",
              organism," target=\"icl_window\"\">",id,"</a>", sep = "")
}


createIM2 <- function(im, file = "", imgTags = list(),
                      title = "Genome View") {
    cat(paste("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n",
              "<html> <head> <title>", title, "</title></head><body>"),
        file = file)
    cat(buildIM(im, imgTags), sep = "\n", file = file, append = TRUE)
    cat("</body></html>", file = file, append = TRUE)
}






                 
    




### There is a lot of repetition in the plotting code. Could place a lot
##  into a function

plot.adacgh.nonsuperimpose <- function(res, chrom,  arraynum, main, colors,
                                       ylim, geneNames, idtype, organism,
                                       geneLoc) {
    plot.adacgh.genomewide(res, chrom,  arraynum, main, colors,
                           ylim, geneNames, geneLoc)
    plot.adacgh.chromosomewide(res, chrom,  arraynum, main, colors,
                               ylim, geneNames, idtype, organism, geneLoc)
}

plot.adacgh.genomewide <- function(res, chrom,
                                   arraynum, main = NULL,
                                   colors = c("orange", "red", "green", "blue"),
                                   ylim = NULL,
                                   geneNames = positions.merge1$name,
                                   geneLoc = NULL) {
    pch <- 20
    im1 <- mapGenomeWideOpen(main)
    nameIm <- main
    logr <- res[[arraynum]][, 1]
    smoothdat <- res[[arraynum]][, 2]
    if(is.null(geneLoc)) {
        simplepos <- (1:length(logr))
    } else {
        ## geneLoc is withing chromosome,
        ## thus, we need some absolute, increasing pos.
        lchr <- tapply(geneLoc, chrom, length)
        mchr <- cumsum(tapply(geneLoc, chrom, max)) ## not very efficient
        sumpos <- rep(c(0, mchr[-length(mchr)]),
                      lchr)
        simplepos <- geneLoc + sumpos
    }
                        
    res.dat <- res[[arraynum]][, 3]
    col <- rep(colors[1],length(res.dat))
    col[which(res.dat == -1)] <- colors[3]
    col[which(res.dat == 1)] <- colors[2]
  
    environment(plotGenomeWide) <- environment(mapLinkChrom) <- environment()
    plotGenomeWide()
    im1 <- mapLinkChrom()
   
    lines(smoothdat ~ simplepos, col="black",
          lwd = 2)
  
    mapGenomeWideClose(nameIm, im1)
}



plot.adacgh.chromosomewide <- function(res, chrom,
                                    arraynum, main = NULL,
                                    colors = c("orange", "red", "green", "blue"),
                                    ylim = NULL,
                                    geneNames = positions.merge1$name,
                                    idtype = idtype,
                                    organism = organism,
                                    geneLoc = NULL) {

    pch <- 20
    logr <- res[[arraynum]][, 1]
    smoothdat <- res[[arraynum]][, 2]
    res.dat <- res[[arraynum]][, 3]
    simplepos <- if(is.null(geneLoc)) (1:length(logr)) else geneLoc
   
    col <- rep(colors[1],length(res.dat))
    col[which(res.dat == -1)] <- colors[3]
    col[which(res.dat == 1)] <- colors[2]
    nameIm <- main
    chrom.nums <- unique(chrom)
    for(cnum in chrom.nums) {
        indexchr <- which(chrom == cnum)
        ccircle <- NULL
        environment(mapChromOpen) <- environment(plotChromWide) <- environment()
        im2 <- mapChromOpen()
        plotChromWide()

        lines(smoothdat[indexchr] ~ simplepos[indexchr],
              col = "black", lwd = 2, type = "l")

        environment(pngCircleRegion) <- environment()
        ccircle <- pngCircleRegion()
        
        environment(mapCloseAndPythonChrom) <- environment()
        mapCloseAndPythonChrom()
    }        ## looping over chromosomes
}



plot.gw.superimp <- function(res, chrom, main = NULL,
                             colors = c("orange", "red", "green", "blue"),
                             ylim =c(ymin, ymax), 
                             geneNames = positions.merge1$name,
                             geneLoc = NULL) {

    pch <- ""
    arraynums <- length(res)
    im1 <- mapGenomeWideOpen(main)
    nameIm <- main
    nfig <- 1
     
    for (arraynum in 1:arraynums) {
        logr <- res[[arraynum]][, 1]
        smoothdat <- res[[arraynum]][, 2]
        if(is.null(geneLoc)) {
            simplepos <- (1:length(logr))
        } else {
            ## geneLoc is withing chromosome,
            ## thus, we need some absolute, increasing pos.
            lchr <- tapply(geneLoc, chrom, length)
            mchr <- cumsum(tapply(geneLoc, chrom, max)) ## not very efficient
            sumpos <- rep(c(0, mchr[-length(mchr)]),
                          lchr)
            simplepos <- geneLoc + sumpos
        }
        res.dat <- res[[arraynum]][, 3]
        col <- rep(colors[1],length(res.dat))
        col[which(res.dat == -1)] <- colors[3]
        col[which(res.dat == 1)] <- colors[2]
 
        if(nfig == 1) {
            environment(plotGenomeWide) <- environment(mapLinkChrom) <- environment()
            plotGenomeWide()
            im1 <- mapLinkChrom()
        }

        lines(smoothdat ~ simplepos,
              col = "black", lwd = 2, type = "l")
        nfig <- nfig + 1
        par(new = TRUE)
    }
    mapGenomeWideClose(nameIm, im1)
}


plot.adacgh.superimp <- function(res, chrom, main,  colors, ylim, geneNames,
                                 idtype, organism, geneLoc) {
    plot.gw.superimp(res, chrom, main,  colors, ylim, geneNames,
                     geneLoc)
    plot.cw.superimp(res, chrom, main,  colors, ylim, geneNames,
                     idtype, organism, geneLoc)
}



plot.cw.superimp <- function(res, chrom, 
                                main = "All_arrays",
                                colors = c("orange", "red", "green", "blue"),
                                ylim =NULL, 
                                geneNames = positions.merge1$name,
                                idtype = idtype, organism = organism,
                                geneLoc = NULL) {
    
    ## For superimposed: one plot per chr
    pch <- ""
    arraynums <- length(res)
    nameIm <- main
    pixels.point <- 3
    chrheight <- 500
    chrom.nums <- unique(chrom)
    for(cnum in chrom.nums) {
        ccircle <- NULL
        environment(mapChromOpen) <- environment()
        im2 <- mapChromOpen()

        indexchr <- which(chrom == cnum)
        
        nfig <- 1
        for(arraynum in 1:arraynums) { ## first, plot the points
##            cat(" ........ for points doing arraynum ", arraynum, "\n")

            logr <- res[[arraynum]][, 1]
            res.dat <- res[[arraynum]][, 3]
            smoothdat <- res[[arraynum]][, 2]
            col <- rep(colors[1],length(res.dat))
            col[which(res.dat == -1)] <- colors[3]
            col[which(res.dat == 1)] <- colors[2]
            simplepos <- if(is.null(geneLoc)) (1:length(logr)) else geneLoc
            
            if(nfig == 1) {
                environment(plotChromWide) <- environment()
                plotChromWide()
            }
            
            environment(pngCircleRegion) <- environment()
            ccircle <- pngCircleRegion()

            ## we want all points, but only draw axes once
            points(logr[indexchr] ~ simplepos[indexchr], col=col[indexchr],
                   cex = 1, pch = 20)
            
##            cat(" ........ for segments doing arraynum ", arraynum, "\n")
            lines(smoothdat[indexchr] ~ simplepos[indexchr],
                  col = "black", lwd = 2, type = "l")
    
            nfig <- nfig + 1
        }
        environment(mapCloseAndPythonChrom) <- environment()
        mapCloseAndPythonChrom()
    }
}





mapCloseAndPythonChrom <- function() {
    nameChrIm <- paste("Chr", chrom.nums[cnum], "@", nameIm, sep ="")
    write(ccircle, file = paste("pngCoordChr_", nameChrIm, sep = ""),
          sep ="\t", ncolumns = 3)
    calcnarrays <- ncol(ccircle)/length(geneNames[indexchr])
    if(!exists("arraynums")) arraynums <- 1
    if(calcnarrays != arraynums)
        stop("Serious problem: number of arrays does not match")
    
    write(rep(as.character(geneNames[indexchr]), arraynums), 
          file = paste("geneNamesChr_", nameChrIm, sep = ""))
    imClose(im2)
    system(paste(.python.toMap.py, nameChrIm, 
                 idtype, organism, sep = " "))
}

    

plotChromWide <- function() {
    par(xaxs = "i")
    par(mar = c(5, 5, 5, 5))
    par(oma = c(0, 0, 0, 0))
    plot(logr[indexchr] ~ simplepos[indexchr], col=col[indexchr], cex = 1,
         xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE,
         main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
         pch = pch, ylim = ylim)
    box()
    axis(2)
    abline(h = 0, lty = 2, col = colors[4])
    rug(simplepos[indexchr], ticksize = 0.01)
}

pngCircleRegion <- function() {
    usr2pngCircle <- function(x, y, rr = 2, rmin = 4) {
        xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
        r <- min(abs(xyrc[2, 1] - xyrc[3, 1]), rmin)
        return(c(xyrc[1, 1], xyrc[1, 2], r))
    } 
    ccircle <- cbind(ccircle,
                     mapply(usr2pngCircle, simplepos[indexchr],
                            logr[indexchr]))
    return(ccircle)
}


mapChromOpen <- function() {
##    cat(" .... doing chromosome ", cnum, "\n")
    nameIm <- main
    pixels.point <- 3
    chrheight <- 500
    indexchr <- which(chrom == chrom.nums[cnum])
    chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
    chrwidth <- max(chrwidth, 800)
    im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     height = chrheight, width = chrwidth,
                     ps = 12)
    return(im2)
}

mapGenomeWideOpen <- function(main) {
    nameIm <- main
    imheight <- 500
    imwidth <- 1600
    im1 <- imagemap3(nameIm, height = imheight,
                     width = imwidth, ps = 12)
    return(im1)
}

mapGenomeWideClose <- function(nameIm, im1) {
    createIM2(im1, file = paste(nameIm, ".html", sep = ""))
    imClose(im1)
}

plotGenomeWide <- function() {
    ## BEWARE!!!: we need to set the environment properly!!
    plot(logr ~ simplepos, col= col, ylab = "log ratio",
         xlab ="Chromosome location", axes = FALSE, cex = 0.7, main = main,
         pch = pch, ylim = ylim)
    box()
    rug(simplepos, ticksize = 0.01)
    axis(2)
    abline(h = 0, lty = 2, col = colors[4])
}


mapLinkChrom <- function() {
    ## 1) add the vertical chromosome lines
    ## 2) html map: add links to chromosome-wide figures
    ## BEWARE!!!: we need to set the environment properly!!
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
        })
        if(class(tryms) == "try-error") {
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

    


old.doMCR <- function(x, chrom.numeric, data,
                  MCR.gapAllowed = 500,
                  MCR.alteredLow = 0.03,
                  MCR.alteredHigh = 0.97,
                  MCR.recurrence = 75,
                  fsink = "results.txt",
                  hsink = "mcr.results.html",
                  ...) {
    if(ncol(data) > 1) {

        delta <- 0.05
        res <- constructSegmObj(x, chrom.numeric, data)
        MCR.unsucc <- FALSE
        keep.running <- TRUE
        while(keep.running) {
            ## Hoping over the bugs in cghMCR
            tryms <- try({
                cghmcr <- cghMCR(res,
                                 gapAllowed = MCR.gapAllowed,
                                 alteredLow = MCR.alteredLow,
                                 alteredHigh = MCR.alteredHigh,
                                 recurrence = MCR.recurrence)
                mcrs <- MCR(cghmcr)
            })
            if(class(tryms) == "try-error") {
                MCR.unsucc <- TRUE
                MCR.alteredLow <- MCR.alteredLow + delta
                MCR.alteredHigh <- MCR.alteredHigh - delta
                warningsForUsers <- c(warningsForUsers,
                                      paste("MCR parameters changed to ",
                                            "alteredLow = ", MCR.alteredLow,
                                            "alteredHigh = ", MCR.alteredHigh,
                                            " to avoid a bug in a BioC pack."))
                
            } else {
                keep.running <- FALSE
            }
        }
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
##        return(mcrsc)

        sink(file = hsink)
        if (nrow(mcrsc) == 0) {
          cat("\n<p> No common minimal regions found.</p>\n")
        } else {
          html.data.frame(mcrsc, first.col = "Case",
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

    


