.__ADaCGH_WEB_APPL <- FALSE ## set to TRUE in web appl!


## .First.lib <- function(lib, pkg) {
## 	   library.dynam("ADaCGH2", pkg, lib)
## 	   }

## Olshen
#library(DNAcopy)
## Price-Smith-Waterman
#library(cgh)
## wavelets
#library("waveslim")
#library("cluster")  
##dyn.load("ADaCGH2.so")

if(.__ADaCGH_WEB_APPL) {
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

mpiInit <- function() {
    library(Rmpi)
    mpi.spawn.Rslaves(nslaves= mpi.universe.size())
    mpi.setup.rngstream() ## or mpi.setup.sprng()
##    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    library(papply)
    mpi.remote.exec(library(cluster))
    mpi.remote.exec(library(waveslim))
    mpi.remote.exec(library(cghMCR))
    mpi.remote.exec(library(DNAcopy))
    mpi.remote.exec(library(cgh))
    mpi.remote.exec(library(ADaCGH))
    

}


pSegmentDNAcopy <- function(x, alpha=0.01, nperm=10000, p.method = c("hybrid","perm"),
                    kmax=25, nmin=200, eta = 0.05, window.size=NULL, overlap=0.25, 
                    trim = 0.025, undo.splits=c("none","prune","sdundo"),
                    undo.prune=0.05, undo.SD=3)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be a copy number array object")
    call <- match.call()
    nsample <- ncol(x)-2
    sampleid <- colnames(x)[-(1:2)]
    uchrom <- unique(x$chrom)
    data.type <- attr(x, "data.type")
    p.method <- match.arg(p.method)
    if (p.method=="hybrid") window.size <- NULL
    undo.splits <- match.arg(undo.splits)
    segres <- list()
    segres$data <- x
    allsegs <- list()
    allsegs$ID <- NULL
    allsegs$chrom <- NULL
    allsegs$loc.start <- NULL
    allsegs$loc.end <- NULL
    allsegs$num.mark <- NULL
    allsegs$seg.mean <- NULL
    ## we parallelize over subjects


    datalist <- list()
    klist <- 1
    for(i in 1:nsample) {
        datalist[[klist]] <- x[, i + 2]
            klist <- klist + 1
    }

    
    funcbs <- function(genomdati) {
        ina <- which(!is.na(genomdati) & !(abs(genomdati)==Inf))
        genomdati <- genomdati[ina]
        trimmed.SD <- sqrt(trimmed.variance(genomdati, trim))
        chromi <- chrom[ina]
                                        #      maploci <- x$maploc[ina]
        sample.lsegs <- NULL
        sample.segmeans <- NULL
        if (nperm == 10000 & alpha == 0.01 & eta == 0.05) {
            sbdry <- default.DNAcopy.bdry
        } else {
            max.ones <- floor(nperm * alpha) + 1
            sbdry <- getbdry(eta, nperm, max.ones)
        }
        sbn <- length(sbdry)
        for (ic in uchrom) {
            segci <- changepoints(genomdati[chromi==ic], data.type, alpha, 
                                  sbdry, sbn, nperm, p.method, window.size, overlap, kmax,
                                  nmin, trimmed.SD, undo.splits, undo.prune,
                                  undo.SD, verbose)
            sample.lsegs <- c(sample.lsegs, segci$lseg)
            sample.segmeans <- c(sample.segmeans, segci$segmeans)
        }
        sample.nseg <- length(sample.lsegs)
        sample.segs.start <- ina[cumsum(c(1,sample.lsegs[-sample.nseg]))]
        sample.segs.end <- ina[cumsum(sample.lsegs)]
        
        chrom.o <- chrom[sample.segs.end]
        loc.start <- maploc[sample.segs.start]
        loc.end <- maploc[sample.segs.end]
        num.mark <- sample.lsegs
        seg.mean <- sample.segmeans

        return(list(chrom.o = chrom.o,
                    loc.start = loc.start,
                    loc.end = loc.end,
                    num.mark = num.mark,
                    seg.mean = seg.mean,
                    nseg = sample.nseg))
    }

    papout <- papply(datalist, funcbs,
                     papply_commondata =list(chrom = x$chrom,
                     uchrom=uchrom,
                     data.type = data.type,
                     alpha = alpha, 
                     nperm = nperm,
                     p.method = p.method,
                     window.size = window.size,
                     overlap = overlap,
                     kmax = kmax,
                     nmin = nmin,
                     trim = trim,
                     eta = eta,
                     undo.splits = undo.splits,
                     undo.prune = undo.prune,
                     undo.SD = undo.SD,
                     verbose = 2,
                     maploc = x$maploc),
                     do_trace = TRUE)                     

    allsegs$chrom <- unlist(lapply(papout, function(x) x$chrom.o))
    allsegs$loc.start <- unlist(lapply(papout, function(x) x$loc.start))
    allsegs$loc.end <- unlist(lapply(papout, function(x) x$loc.end))
    allsegs$num.mark <- unlist(lapply(papout, function(x) x$num.mark))
    allsegs$seg.mean <- unlist(lapply(papout, function(x) x$seg.mean))
    ## I could have just one statement for ID, but lets try to follow original
    ## code closely
    allsegs$ID <- rep(1:nsample, unlist(lapply(papout, function(x) x$nseg)))

    allsegs$ID <- sampleid[allsegs$ID]
    allsegs$seg.mean <- round(allsegs$seg.mean, 4)
    allsegs <- as.data.frame(allsegs)
    allsegs$ID <- as.character(allsegs$ID)
    segres$output <- allsegs
    segres$call <- call    
    class(segres) <- "DNAcopy"
    segres
}


###### We will use:
## a) segment plots which are same name for all, but call diff. functs.
## b) plateau plots (only for xx)


## a function that will merge the CBS?

mergeDNAcopy <- function(object) {
    numarrays <- ncol(object$data) - 2
    ## zz: we must get numarray from object
    
    if(!(inherits(object, "DNAcopy")))
        stop("This function can only be applied to DNAcopy objects")
    merged_segments <- list()
    merged_segments$chrom.numeric <- object$data$chrom ## verify its numeric zz!!
    merged_segments$segm <- list()
    for(arraynum in 1:numarrays) {
        obs <- object$data[, 2 + arraynum]
        segmented <-
            object$output[object$output$ID ==
                          colnames(object$data)[2 + arraynum], ]
        segmentus <- object$data$maploc
        for(i in 1:nrow(segmented)) {
            segmentus[(segmented[i,'loc.end'] >= segmentus) &
                      (segmented[i,'loc.start'] <= segmentus)] <-
                          segmented[i,'seg.mean']
        }
        segmentus2 <- mergeLevels(obs, segmentus)$vecMerged
        classes.ref <- which.min(abs(unique(segmentus2)))
        classes.ref <- unique(segmentus2)[classes.ref]
        ref <- rep(0, length(segmentus))
        ref[segmentus2 > classes.ref] <- 1
        ref[segmentus2 < classes.ref] <- -1
        merged_segments$segm[[arraynum]] <- cbind(merged.mean = segmentus, obs, alteration = ref)
    }
    class(merged_segments) <- c(class(merged_segments),
                                "mergedDNAcopy")
    return(merged_segments)
}    


pSegmentPSW <- function(common.dat,
                        acghdata,
                        chrom.numeric,
                        sign = -1,
                        nIter = 1000,
                        prec = 100,
                        p.crit = 0.10,
                        name = NULL) {
    numarrays <- ncol(acghdata)
    ncrom <- length(unique(chrom.numeric))
    out <- list()
    out$Data <- common.dat
    out$plotData <- list()
    for(i in 1:numarrays) {
        tmp <- my.sw3(acghdata[, i], chrom.numeric,
                      sign = sign, p.crit = p.crit,
                      main = colnames(acghdata)[i],
                      nIter = nIter,
                      prec = prec,
                      name = paste(name, colnames(acghdata)[i], sep = ""))
        out$Data <- cbind(out$Data, tmp$out)
        p.crit.bonferroni <- tmp$plotdat$p.crit / ncrom
        out$plotData[[i]] <- c(tmp$plotdat,
                              p.crit.bonferroni = p.crit.bonferroni)
    }
    class(out) <- c(class(out), "CGH.PSW")
    return(out)
}
        


segmentPlot <- function(x, geneNames,
                        chrom.numeric = NULL,
                        cghdata = NULL,
                        arraynames = NULL,
                        idtype = "ug",
                        organism = "Hs",
                        superimposed = FALSE,
                        html = TRUE,
                        yminmax = NULL,
                        numarrays = NULL,
                        ...) {
    if(is.null(numarrays)) {
        if(!is.null(arraynames)) numarrays <- length(arraynames)
        if(!is.null(cghdata)) numarrays <- ncol(cghdata)
    }
    if(is.null(yminmax)) {
        yminmax <- c(min(as.matrix(cghdata)),
                     max(as.matrix(cghdata)))
    }
    if(is.null(arraynames)) arraynames <- colnames(cghdata)
    if(is.null(arraynames)) arraynames <- paste("sample.", 1:numarrays, sep = "")
    
    classx <- inherits(x, c("DNAcopy", "CGH.wave", "CGH.PSW",
                            "CGH.ACE.summary", "mergedDNAcopy"))
    if (inherits(x, "DNAcopy")) {
        if(!superimposed) {
            for(i in 1:numarrays) { cat("\n Doing sample ", i, "\n")
                plot.olshen2(x,
                             i, main = arraynames[i],
                             html = TRUE,
                             geneNames = geneNames,
                             idtype = idtype, organism = organism)
            }
        } else {
            plot.olshen3(x, html = html, geneNames = geneNames,
                         main = "All_arrays", ylim = yminmax,
                         html = html,
                         idtype = idtype, organism = organism) ## all genome
            
            plot.olshen4(x,  geneNames = geneNames,
                         main = "All_arrays", ylim = yminmax,
                         html = html,
                         idtype = idtype, organism = organism) ## by chromosome
        }
    } else if(inherits(x, "mergedDNAcopy")) {
        if(!superimposed) {
            for(i in 1:numarrays) { cat("\n Doing sample ", i, "\n")
                plot.ace2(x$segm, x$chrom.numeric, arraynum = i, geneNames = geneNames,
                          main = arraynames[i],
                          idtype = idtype, organism = organism)
            }
        } else {
            plot.ace3(x$segm, x$chrom.numeric,  geneNames = geneNames,
                      main = "All_arrays",
                      ylim = yminmax,
                      pch = "",
                      idtype = idtype, organism = organism)
            plot.ace4(x, chrom.numeric,  geneNames = geneNames,
                      main = "All_arrays",
                      ylim = yminmax,
                      idtype = idtype, organism = organism)
        }
    } else if(inherits(x, "CGH.ACE.summary")) {
        if(!superimposed) {
            for(i in 1:numarrays) { cat("\n Doing sample ", i, "\n")
                plot.ace2(x, chrom.numeric, arraynum = i, geneNames = geneNames,
                          main = arraynames[i],
                          idtype = idtype, organism = organism)
            }
        } else {
            plot.ace3(x, chrom.numeric,  geneNames = geneNames,
                      main = "All_arrays",
                      ylim = yminmax,
                      pch = "",
                      idtype = idtype, organism = organism,
                      arraynums = 1:numarrays)
            plot.ace4(x, chrom.numeric,  geneNames = geneNames,
                      main = "All_arrays",
                      ylim = yminmax,
                      idtype = idtype, organism = organism,
                      arraynums = 1:numarrays)
        }
    } else if(inherits(x, "CGH.wave")) {
        if(!superimposed) {
            for(i in 1:numarrays) { cat("\n Doing sample ", i, "\n")
                plot.wavelets2(x, cghdata, chrom.numeric, arraynum = i,
                               geneNames = geneNames,
                               main = arraynames[i],
                               idtype = idtype, organism = organism)
            }
        } else {
            plot.wavelets3(x, cghdata, chrom.numeric,  geneNames = geneNames,
                           main = "All_arrays",
                           ylim = yminmax,
                           pch = "",
                           idtype = idtype, organism = organism)
            plot.wavelets4(x, cghdata, chrom.numeric,  geneNames = geneNames,
                           main = "All_arrays",
                           ylim = yminmax,
                           idtype = idtype, organism = organism)
        }
    } else if(inherits(x, "CGH.PSW")) {
        if(x$plotData[[i]]$sign < 0) {
            main <- "Losses."
        } else {
            main <- "Gains."
        }
        for(i in 1:numarrays) { cat("\n Doing sample ", i, "\n")
            sw.plot3(x$plotData[[i]]$logratio, sign=x$plotData[[i]]$sign,
                     swt.perm = x$plotData[[i]]$swt.perm, rob = x$plotData[[i]]$rob,
                     swt.run = x$plotData[[i]]$swt.run,
                     p.crit = x$plotData[[i]]$p.crit.bonferroni,
                     chrom = x$plotData[[i]]$chrom,
                     main = paste(main, arraynames[i], sep = ""),
                     geneNames = geneNames,
                     idtype = idtype, organism = organism)
        }
    } else {
        stop("No plotting for this class of objects")
    }
    
}
  

plateauPlot <- function(obj, ...) {
    UseMethod("plateauplot", ...)
}

plateauPlot.DNAcopy <- function() {}
plateauPlot.CGH.wave <- function() {}

writeResults <- function(obj, ...) {
    UseMethod("writeResults", ...)
}

writeResults.CGH.PSW <- function(obj, file = "PSW.output.txt") {
    write.table(obj, file = file,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
}

writeResults.CGH.ACE.summary <- function(obj, file = "ACE.output.txt") {
    print.ACE.results(obj, output = file)
}

writeResults.CGH.wave <- function(obj, acghdata, commondata,
                                  file = "wavelet.output.txt") {
    print.wavelets.results(obj, acghdata, commondata, output = file)
}

writeResults.DNAcopy <- function(obj, acghdata, commondata, merged = NULL,
                                 output = "CBS.output.txt") {
    print.olshen.results(obj, acghdata, commondata,
                         merged = merged, output = output) 
}


pSegmentWavelets <- function(acghdata, chrom.numeric, minDiff = 0.25, thrLvl = 3,
                             initClusterLevels = 10) {
## level to use for wavelet decomposition and thresholding
## The 'recommended' level is floor(log2(log(N)+1)), which
## equals 3 for:  21 <= N <= 1096
##    thrLvl <- 3

    ncloneschrom <- tapply(acghdata[, 1], chrom, function(x) length(x))
    if((thrLvl == 3) & ((max(ncloneschrom) > 1096) | (min(ncloneschrom) < 21)))
        warningsForUsers <-
            c(warningsForUsers,
              paste("The number of clones/genes is either",
                    "larger than 1096 or smaller than 21",
                    "in at least one chromosome. The wavelet",
                    "thresholding of 3 might not be appropriate."))
    
    Nsamps  <- ncol(acghdata)
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
        ratio.i <- acghdata[,i]
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
        pred.ij <- segmentW(ratio, recH, minDiff=minDiff, n.levels = initClusterLevels)
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
    class(out) <- c(class(out), "waveCGH")
    return(out)
}

pSegmentACE <- function(acghdata, chrom.numeric,  echo=FALSE) {
    ACE(acghdata, chrom.numeric, echo, coefs = file.aux, Sdev = 0.2)
}





######################################

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
                     main = paste(colnames(CNA.object[i+2]), "; Chr ", j, sep = ""),
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
            if(class(trythis) == "try-error")
                caughtOurError(paste("acf bombed unexpectedly with error",
                                     trythis, ". \n Please let us know so we can fix the code."))

        }
    }
    rm(cn, subject, index.dat)
    ##         pdf("Autocorrelation.plots.pdf", width = 17.6, height = 12.5)
    ##        par(mfrow = c(6,4))
    ##        par(oma = c(2, 2, 2, 2))
    for(i in 1:ncrom)
        hist(ar1s[, i], main = paste("Chr", i), xlab = "Autocorrelation, lag 1")
    mtext("Autocorrelation coefficient lag 1", side = 3,
          outer = TRUE, line = 0.5, font = 2)
    ##        dev.off()



##### example
##### data(cghMCRe)

##### chrom.numeric <- as.numeric(as.character(cghMCRe$Chromosome))
##### chrom.numeric[cghMCRe$Chromosome == "X"] <- 23
##### chrom.numeric[cghMCRe$Chromosome == "Y"] <- 24

##### par(mfrow = c(5, 5))
##### WaveletsDiagnosticPlots(cghMCRe[, 5:7], chrom.numeric)
}







caughtOurError <- function(message) {
    if(.__ADaCGH_WEB_APPL) {
        webPNG("ErrorFigure.png", width = png.width,
               height = png.height, 
               pointsize = png.pointsize,
               family = png.family)
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
    


    
######################################################

#########  Imagemap stuff



##require(GDD)
##require(imagemap)

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






mpiCBS <- function() {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(cghMCR))
    mpi.remote.exec(library(cluster))
    mpi.remote.exec(library(DNAcopy))
    mpi.remote.exec(library(ADaCGH2))
    mpi.remote.exec(library(aCGH))
}

plot.olshen2 <- function(res, arraynum, main = NULL,
                         colors = c("orange", "red", "green", "blue"),
                         pch = 20, ylim =NULL, html = TRUE,
                         superimpose = FALSE,
                         nsupimp = 0,
                         geneNames = positions.merge1$name,
                         idtype = idtype,
                         organism = organism) {
                                        #res is the results
                                        # color code for region status
## Like plot.olshen2, but with identifiers and imagemap
    logr <- res$data[, 2 + arraynum]
    segmented <-
        res$output[res$output$ID == colnames(res$data)[2 + arraynum], ]
    col <- rep(colors[1],length(logr))
##      col[which(res$profileValues$ZoneGNL==-1)] <- colors[2]
##      col[which(res$profileValues$ZoneGNL==1)] <- colors[3]
##      # outliers
##      outliers <- rep(20,length(res$profileValues$PosOrder))
##      outliers[which(res$profileValues$OutliersTot!=0)] <- 13



    nameIm <- main
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
                         width = imwidth, ps = 12)
    }

    plot(logr ~ res$data$maploc, col="orange", 
          xlab ="Chromosomal location", axes = FALSE, cex = 0.7, main = main,
          pch = pch, ylim = ylim)
     box()
     axis(2)
    abline(h = 0, lty = 2, col = "blue")
    rug(res$data$maploc, ticksize = 0.01)
    
    ## Limit between chromosomes
    LimitChr <- tapply(res$data$maploc,
                       res$data$chrom, max)
    abline(v=LimitChr, col="grey", lty=2)

    chrom.nums <- unique(res$data$chrom)
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)

    ## segments
    for(j in 1:nrow(segmented)) {
        segments(x0 = segmented$loc.start[j],
                 y0 = segmented$seg.mean[j],
                 x1 = segmented$loc.end[j],
                 y1 = segmented$seg.mean[j],
                 col = "black", lwd = 2)
    }

    if(html) {
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
        createIM2(im1, file = paste(nameIm, ".html", sep = ""))
        imClose(im1)
    }


    if(html) { ## here is chromosome specific code
        pixels.point <- 3
        chrheight <- 500
        chrom <- res$data$chrom
        for(cnum in 1:length(chrom.nums)) {
            cat(" .... doing chromosome ", cnum, "\n")
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
            plot(logr[indexchr] ~ res$data$maploc[indexchr], col="orange", 
                 xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE, cex = 1,
                 main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                 pch = pch, ylim = ylim)
            box()
            axis(2)
            abline(h = 0, lty = 2, col = "blue")
            rug(res$data$maploc[indexchr], ticksize = 0.01)            
            ## segments
            for(j in 1:nrow(segmented)) {
                if( (segmented$loc.start[j] >= res$data$maploc[indexchr[1]])
                   & ( segmented$loc.start[j] <= res$data$maploc[indexchr[length(indexchr)]])) {
                    segments(x0 = segmented$loc.start[j],
                             y0 = segmented$seg.mean[j],
                             x1 = segmented$loc.end[j],
                             y1 = segmented$seg.mean[j],
                             col = "black", lwd = 2)
                }
            }

            ## The within chromosome map for gene names
            ## in superimposed cases only in first map
            usr2pngCircle <- function(x, y, rr = 2) {
                xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
                r <- abs(xyrc[2, 1] - xyrc[3, 1])
                return(c(xyrc[1, 1], xyrc[1, 2], r))
            }
            
            ccircle <- mapply(usr2pngCircle, res$data$maploc[indexchr],
                              logr[indexchr])
            write(ccircle, file = "pngCoordChr",
                  sep ="\t", ncolumns = 3)
            write(as.character(geneNames[indexchr]), file = "geneNamesChr")
            imClose(im2)
            system(paste(.python.toMap.py,
                         paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         idtype, organism, sep = " "))
        }        ## looping over chromosomes
    } ## if html
}

   
plot.olshen3 <- function(res, arraynums = 1:numarrays, main = NULL,
                         colors = c("orange", "red", "green", "blue"),
                         pch = 20, ylim =NULL, html = TRUE,
                         geneNames = positions.merge1$name,
                         idtype = idtype, organism = organism) {
    ## For superimposed: only all genome plot with map to chromosomes

    nameIm <- main
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
                         width = imwidth, ps = 12)
    }
    nfig <- 1
    for (arraynum in arraynums) {
        logr <- res$data[, 2 + arraynum]
        segmented <-
            res$output[res$output$ID == colnames(res$data)[2 + arraynum], ]
        col <- rep(colors[1],length(logr))
        if(nfig == 1) {
            plot(logr ~ res$data$maploc, col="orange", 
                 xlab ="Chromosomal location", axes = FALSE, cex = 0.7, main = main,
                 pch = "", ylim = ylim)
            box()
            axis(2)
            abline(h = 0, lty = 2, col = "blue")
            rug(res$data$maploc, ticksize = 0.01)
        }

        for(j in 1:nrow(segmented)) {
            segments(x0 = segmented$loc.start[j],
                     y0 = segmented$seg.mean[j],
                     x1 = segmented$loc.end[j],
                     y1 = segmented$seg.mean[j],
                     col = "black", lwd = 2)
        }

        if(nfig == 1) {
            ## Limit between chromosomes
            LimitChr <- tapply(res$data$maploc,
                               res$data$chrom, max)
            abline(v=LimitChr, col="grey", lty=2)
            
            chrom.nums <- unique(res$data$chrom)
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
        }
        nfig <- nfig + 1
        par(new = TRUE)
    }
    createIM2(im1, file = paste(nameIm, ".html", sep = ""))
    imClose(im1)
}








plot.olshen4 <- function(res, arraynums = 1:numarrays, main = NULL,
                         colors = c("orange", "red", "green", "blue"),
                         pch = 20, ylim =NULL, html = TRUE,
                         geneNames = positions.merge1$name,
                         idtype = idtype,
                         organism = organism) {
    ## For superimposed: one plot per chr

    nameIm <- main
    
    pixels.point <- 3
    chrheight <- 500
    chrom <- res$data$chrom
    chrom.nums <- unique(chrom)
    for(cnum in 1:length(chrom.nums)) {
        cat(" .... doing chromosome ", cnum, "\n")
        indexchr <- which(chrom == chrom.nums[cnum])
        chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
        chrwidth <- max(chrwidth, 800)
        im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         height = chrheight, width = chrwidth,
                         ps = 12)
        ## The following seems needed (also inside sw.plot2) for the coords.
        ## of points to work OK
        nfig <- 1
        for(arraynum in arraynums) { ## first, plot the points
            cat(" ........ for points doing arraynum ", arraynum, "\n")
            logr <- res$data[, 2 + arraynum]
            if(nfig == 1) {
                par(xaxs = "i")
                par(mar = c(5, 5, 5, 5))
                par(oma = c(0, 0, 0, 0))
                plot(logr[indexchr] ~ res$data$maploc[indexchr], col="orange", 
                     xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE, cex = 1,
                     main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     pch = "", ylim = ylim)
                box()
                axis(2)
                abline(h = 0, lty = 2, col = "blue")
                rug(res$data$maploc[indexchr], ticksize = 0.01)
            }
            points(logr[indexchr] ~ res$data$maploc[indexchr], col="orange",
                   cex = 1, pch = 20)
            nfig <- nfig + 1
        }
        for(arraynum in arraynums) { ## now, do the segments
            cat(" ........ for segments doing arraynum ", arraynum, "\n")
            segmented <-
                res$output[res$output$ID == colnames(res$data)[2 + arraynum], ]
            
            
            ## segments
            for(j in 1:nrow(segmented)) {
                if( (segmented$loc.start[j] >= res$data$maploc[indexchr[1]])
                   & ( segmented$loc.start[j] <= res$data$maploc[indexchr[length(indexchr)]])) {
                    segments(x0 = segmented$loc.start[j],
                             y0 = segmented$seg.mean[j],
                             x1 = segmented$loc.end[j],
                             y1 = segmented$seg.mean[j],
                             col = "black", lwd = 2)
                }
            }
        }

        usr2pngCircle <- function(x, y, rr = 2) {
            xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
            r <- abs(xyrc[2, 1] - xyrc[3, 1])
            return(c(xyrc[1, 1], xyrc[1, 2], r))
        }
        ccircle <- mapply(usr2pngCircle, res$data$maploc[indexchr],
                          0)
        write(ccircle, file = "pngCoordChr",
              sep ="\t", ncolumns = 3)
        write(as.character(geneNames[indexchr]), file = "geneNamesChr")
        imClose(im2)
        system(paste(.python.toMap.py,
                     paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     idtype, organism, sep = " "))
    }
}



### Plateau plots by Olshen 

print.olshen.results <- function(res, xcenter,
                                 commondata,
                                 merged = NULL,
                                 output = "CBS.results.txt"){
    ## This function "stretches out" the output and creates a table
    ## that matches the original names, etc.

    out <- data.frame(ID = commondata$name,
                      Chromosome = commondata$chromosome,
                      Start = commondata$start,
                      End = commondata$end,
                      MidPoint = commondata$MidPoint)

    for(i in 1:ncol(xcenter)) {
        tmp <- res$output[res$output$ID == colnames(res$data)[2 + i], ]
        t1 <- rep(tmp$seg.mean, tmp$num.mark)
        t2 <- xcenter[, i]
        out <- cbind(out, t2, t1)
        if(!is.null(merged)) {
            out <- cbind(out, merged[[i]][ , 3])
        }
    }
    if(is.null(merged)) {
        colnames(out)[6:(ncol(out))] <-
            paste(rep(colnames(xcenter),rep(2, ncol(xcenter))),
                  c(".Original", ".Smoothed"), sep = "")
    } else {
        colnames(out)[6:(ncol(out))] <-
            paste(rep(colnames(xcenter),rep(3, ncol(xcenter))),
                  c(".Original", ".Smoothed", ".Merged"), sep = "")
    }
    write.table(out, file = output,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
}

    
    



plot.DNAcopy2 <- function (x, plot.type = "plateau", xmaploc = FALSE, altcol = TRUE, sbyc.layout = NULL, 
    cbys.nchrom = 1, cbys.layout = NULL, include.means = TRUE, 
    zeroline = TRUE, pt.pch = ".", pt.cex = NULL, pt.cols = NULL, 
    segcol = NULL, zlcol = NULL, ylim = NULL, lwd = NULL, ...) 
{
    if (!inherits(x, "DNAcopy")) 
        stop("First arg must be the result of segment")
    xdat <- x$data
    nsample <- ncol(xdat) - 2
    if (missing(ylim)) {
        uylim <- max(abs(xdat[, -(1:2)]), na.rm = TRUE)
        ylim <- c(-uylim, uylim)
    }
    xres <- x$output
    if (dev.cur() <= 1) 
        get(getOption("device"))()
    int.dev <- dev.interactive()
    plot.type <- match.arg(plot.type)
    op <- par(no.readonly = TRUE)
    parask <- par("ask")
##     if (int.dev & !parask & nsample > 1) 
##         par(ask = TRUE)
    sampleid <- colnames(xdat)[-(1:2)]
    chrom0 <- xdat$chrom
    uchrom <- unique(chrom0)
    nchrom <- length(uchrom)
    if (xmaploc) {
        maploc0 <- as.numeric(xdat$maploc)
        if (max(maploc0[chrom0 == uchrom[1]]) > min(maploc0[chrom0 == 
            uchrom[2]])) {
            plen <- max(maploc0[chrom0 == uchrom[1]])
            for (i in 2:nchrom) {
                maploc0[chrom0 == uchrom[i]] <- plen + maploc0[chrom0 == 
                  uchrom[i]]
                plen <- max(maploc0[chrom0 == uchrom[i]])
            }
        }
    }
    if (missing(pt.pch)) 
        pt.pch <- "."
    if (missing(pt.cex)) {
        if (pt.pch == ".") {
            pt.cex <- 3
        }
        else {
            pt.cex <- 1
        }
    }
    wcol0 <- rep(1, length(chrom0))
    if (altcol) {
        j <- 0
        for (i in uchrom) {
            j <- (j + 1)%%2
            wcol0[chrom0 == i] <- 1 + j
        }
    }
    if (missing(pt.cols)) 
        pt.cols <- c("black", "green")
    if (missing(segcol)) 
        segcol <- "red"
    if (missing(zlcol)) 
        zlcol <- "grey"
    if (missing(lwd)) 
        lwd <- 3
    if (plot.type == "chrombysample") {
        cat("Setting multi-figure configuration\n")
        par(mar = c(0, 4, 0, 2), oma = c(4, 0, 4, 0), mgp = c(2, 
            0.7, 0))
        if (missing(cbys.layout)) {
            nrow <- ncol <- ceiling(sqrt(nsample))
            if (nrow * ncol - nsample > 0) {
                nrow <- nrow - 1
                ncol <- ncol + 1
            }
            if (nrow * ncol - nsample >= nrow) 
                ncol <- ncol - 1
            cbys.layout <- c(nrow, ncol)
        }
        lmat0 <- lmat1 <- c(1:nsample, rep(-cbys.nchrom * nsample, 
            prod(cbys.layout) - nsample))
        for (i in 1:(cbys.nchrom - 1)) {
            lmat1 <- c(lmat1, lmat0 + nsample * i)
        }
        lmat1[lmat1 < 0] <- 0
        lmat <- matrix(lmat1, nrow = cbys.layout[1], ncol = cbys.nchrom * 
            cbys.layout[2], byrow = FALSE)
        layout(lmat)
    }
    if (plot.type == "samplebychrom") {
        cat("Setting multi-figure configuration\n")
        par(mar = c(4, 4, 4, 2), oma = c(0, 0, 2, 0), mgp = c(2, 
            0.7, 0))
        if (missing(sbyc.layout)) {
            nrow <- ncol <- ceiling(sqrt(nchrom))
            if (nrow * ncol - nchrom > 0) {
                nrow <- nrow - 1
                ncol <- ncol + 1
            }
            if (nrow * ncol - nchrom > ncol) 
                ncol <- ncol - 1
            sbyc.layout <- c(nrow, ncol)
        }
        lmat <- matrix(c(1:nchrom, rep(0, prod(sbyc.layout) - 
            nchrom)), nrow = sbyc.layout[1], ncol = sbyc.layout[2], 
            byrow = TRUE)
        layout(lmat)
    }
    if (plot.type == "chrombysample") {
        atchrom <- 0.5/cbys.nchrom
        for (ichrom in uchrom) {
            for (isamp in 1:nsample) {
                genomdat <- xdat[chrom0 == ichrom, isamp + 2]
                ina <- which(!is.na(genomdat) & !(abs(genomdat) == 
                  Inf))
                genomdat <- genomdat[ina]
                ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp] & 
                  xres$chrom == ichrom]))
                mm <- xres$seg.mean[xres$ID == sampleid[isamp] & 
                  xres$chrom == ichrom]
                kk <- length(ii)
                zz <- cbind(ii[-kk] + 1, ii[-1])
                plot(genomdat, pch = pt.pch, cex = pt.cex, xaxt = "n", 
                  ylim = ylim, ylab = sampleid[isamp])
                if (zeroline) 
                  abline(h = 0, col = zlcol, lwd = lwd)
                if (isamp%%cbys.layout[1] == 0) {
                  axis(1, outer = TRUE)
                  title(xlab = "Index")
                }
                if (include.means) {
                  for (i in 1:(kk - 1)) {
                    lines(zz[i, ], rep(mm[i], 2), col = segcol, 
                      lwd = lwd)
                  }
                }
            }
            mtext(paste("Chromosome", ichrom), side = 3, line = 1, 
                at = atchrom, outer = TRUE, font = 2)
            atchrom <- atchrom + 1/cbys.nchrom
            atchrom <- atchrom - floor(atchrom)
        }
    }
    else {
        for (isamp in 1:nsample) {
            genomdat <- xdat[, isamp + 2]
            ina <- which(!is.na(genomdat) & !(abs(genomdat) == 
                Inf))
            genomdat <- genomdat[ina]
            wcol <- wcol0[ina]
            chrom <- chrom0[ina]
            if (xmaploc) 
                maploc <- maploc0[ina]
            ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp]]))
            mm <- xres$seg.mean[xres$ID == sampleid[isamp]]
            kk <- length(ii)
            zz <- cbind(ii[-kk] + 1, ii[-1])
            if (missing(ylim)) 
                ylim <- range(c(genomdat, -genomdat))
            if (plot.type == "whole") {
                if (xmaploc) {
                  plot(maploc, genomdat, pch = pt.pch, cex = pt.cex, 
                    col = pt.cols[wcol], main = sampleid[isamp], 
                    ylab = "", ylim = ylim)
                  if (zeroline) 
                    abline(h = 0, col = zlcol, lwd = lwd)
                }
                else {
                  plot(genomdat, pch = pt.pch, cex = pt.cex, 
                    col = pt.cols[wcol], main = sampleid[isamp], 
                    ylab = "", ylim = ylim)
                  if (zeroline) 
                    abline(h = 0, col = zlcol, lwd = lwd)
                }
                if (include.means) {
                  for (i in 1:(kk - 1)) {
                    if (xmaploc) {
                      lines(maploc[zz[i, ]], rep(mm[i], 2), col = segcol, 
                        lwd = lwd)
                    }
                    else {
                      lines(zz[i, ], rep(mm[i], 2), col = segcol, 
                        lwd = lwd)
                    }
                  }
                }
            }
            if (plot.type == "samplebychrom") {
                cc <- xres$chrom[xres$ID == sampleid[isamp]]
                for (ichrom in uchrom) {
                  plot(genomdat[chrom == ichrom], pch = pt.pch, 
                    cex = pt.cex, ylab = "", main = paste("Chromosome", 
                      ichrom), ylim = ylim)
                  if (zeroline) 
                    abline(h = 0, col = zlcol, lwd = lwd)
                  if (include.means) {
                    jj <- which(cc == ichrom)
                    jj0 <- min(jj)
                    for (i in jj) {
                      lines(1 + zz[i, ] - zz[jj0, 1], rep(mm[i], 
                        2), col = segcol, lwd = lwd)
                    }
                  }
                }
                mtext(sampleid[isamp], side = 3, line = 0, outer = TRUE, 
                  font = 2)
            }
            if (plot.type == "plateau") {
                omm <- order(mm)
                ozz <- zz[omm, ]
                ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
                plot(genomdat[ina], pch = ".", cex = pt.cex, 
                  main = sampleid[isamp], ylab = "", ylim = ylim, col = "orange")
                if (zeroline) 
                  abline(h = 0, col = zlcol, lwd = lwd)
                if (include.means) {
                  ii <- cumsum(c(0, xres$num.mark[xres$ID == 
                    sampleid[isamp]][omm]))
                  smm <- mm[omm]
                  zz <- cbind(ii[-kk] + 1, ii[-1])
                  for (i in 1:(kk - 1)) lines(zz[i, ], rep(smm[i], 
                    2), col = segcol, lwd = lwd)
                }
            }
        }
    }
    on.exit(if (plot.type == "chrombysample" | plot.type == "samplebychrom") {
        par(op)
    } else {
        if (int.dev & !parask & nsample > 1) par(ask = parask)
    })
}





## one single plateau plot
plot.DNAcopy4 <- function (x, plot.type = "plateau",
                            xmaploc = FALSE, altcol = TRUE, sbyc.layout = NULL, 
    cbys.nchrom = 1, cbys.layout = NULL, include.means = TRUE, 
    zeroline = TRUE, pt.pch = NULL, pt.cex = NULL, pt.cols = NULL, 
    segcol = NULL, zlcol = NULL, ylim = NULL, lwd = NULL, ...) 
{
    if (!inherits(x, "DNAcopy")) 
        stop("First arg must be the result of segment")
    xdat <- x$data
    nsample <- ncol(xdat) - 2
    if (missing(ylim)) {
        uylim <- max(abs(xdat[, -(1:2)]), na.rm = TRUE)
        ylim <- c(-uylim, uylim)
    }
    xres <- x$output
    if (dev.cur() <= 1) 
        get(getOption("device"))()
    int.dev <- dev.interactive()
    plot.type <- match.arg(plot.type)
    op <- par(no.readonly = TRUE)
    parask <- par("ask")
##     if (int.dev & !parask & nsample > 1) 
##         par(ask = TRUE)
    sampleid <- colnames(xdat)[-(1:2)]
    chrom0 <- xdat$chrom
    uchrom <- unique(chrom0)
    nchrom <- length(uchrom)
    if (xmaploc) {
        maploc0 <- as.numeric(xdat$maploc)
        if (max(maploc0[chrom0 == uchrom[1]]) > min(maploc0[chrom0 == 
            uchrom[2]])) {
            plen <- max(maploc0[chrom0 == uchrom[1]])
            for (i in 2:nchrom) {
                maploc0[chrom0 == uchrom[i]] <- plen + maploc0[chrom0 == 
                  uchrom[i]]
                plen <- max(maploc0[chrom0 == uchrom[i]])
            }
        }
    }
    if (missing(pt.pch)) 
        pt.pch <- "."
    if (missing(pt.cex)) {
        if (pt.pch == ".") {
            pt.cex <- 3
        }
        else {
            pt.cex <- 1
        }
    }
    wcol0 <- rep(1, length(chrom0))
    if (altcol) {
        j <- 0
        for (i in uchrom) {
            j <- (j + 1)%%2
            wcol0[chrom0 == i] <- 1 + j
        }
    }
    if (missing(pt.cols)) 
        pt.cols <- c("black", "green")
    if (missing(segcol)) 
        segcol <- "red"
    if (missing(zlcol)) 
        zlcol <- "grey"
    if (missing(lwd)) 
        lwd <- 3
    genomdat <- as.vector(as.matrix(xdat[ , -c(1, 2)]))
    ii <- cumsum(c(0, xres$num.mark))
    mm <- xres$seg.mean
    kk <- length(ii)
    zz <- cbind(ii[-kk] + 1, ii[-1])
    if (missing(ylim)) 
        ylim <- range(c(genomdat, -genomdat))
    omm <- order(mm)
    ozz <- zz[omm, ]
    ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
    plot(genomdat[ina], pch = ".", cex = 0.7, 
         main = "All arrays", ylab = "", ylim = ylim, col = "orange")
    if (zeroline) 
        abline(h = 0, col = zlcol, lwd = lwd)


##    mm <- xres$seg.mean
##    omm <- order(mm)
    ii <- cumsum(c(0, xres$num.mark[omm]))
##    kk <- length(ii)
    smm <- mm[omm]
    zz <- cbind(ii[-kk] + 1, ii[-1])
    ozz <- zz[omm, ]
##    ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
    
    for (i in 1:(kk - 1))
        lines(zz[i, ], rep(smm[i], 2),
              col = segcol, lwd = lwd)
    
    
    
    on.exit(if (plot.type == "chrombysample" | plot.type == "samplebychrom") {
        par(op)
    } else {
        if (int.dev & !parask & nsample > 1) par(ask = parask)
    })
}




## Superimposes all chromosome plots
plot.DNAcopy3 <- function (x, plot.type = c("whole", "plateau", "samplebychrom", 
    "chrombysample"), xmaploc = FALSE, altcol = TRUE, sbyc.layout = NULL, 
    cbys.nchrom = 1, cbys.layout = NULL, include.means = TRUE, 
    zeroline = TRUE, pt.pch = NULL, pt.cex = NULL, pt.cols = NULL, 
    segcol = NULL, zlcol = NULL, ylim = NULL, lwd = NULL, ...) 
{
    if (!inherits(x, "DNAcopy")) 
        stop("First arg must be the result of segment")
    xdat <- x$data
    nsample <- ncol(xdat) - 2
    if (missing(ylim)) {
        uylim <- max(abs(xdat[, -(1:2)]), na.rm = TRUE)
        ylim <- c(-uylim, uylim)
    }
    xres <- x$output
    if (dev.cur() <= 1) 
        get(getOption("device"))()
    int.dev <- dev.interactive()
    plot.type <- match.arg(plot.type)
    op <- par(no.readonly = TRUE)
    parask <- par("ask")
##     if (int.dev & !parask & nsample > 1) 
##         par(ask = TRUE)
    sampleid <- colnames(xdat)[-(1:2)]
    chrom0 <- xdat$chrom
    uchrom <- unique(chrom0)
    nchrom <- length(uchrom)
    if (xmaploc) {
        maploc0 <- as.numeric(xdat$maploc)
        if (max(maploc0[chrom0 == uchrom[1]]) > min(maploc0[chrom0 == 
            uchrom[2]])) {
            plen <- max(maploc0[chrom0 == uchrom[1]])
            for (i in 2:nchrom) {
                maploc0[chrom0 == uchrom[i]] <- plen + maploc0[chrom0 == 
                  uchrom[i]]
                plen <- max(maploc0[chrom0 == uchrom[i]])
            }
        }
    }
    if (missing(pt.pch)) 
        pt.pch <- "."
    if (missing(pt.cex)) {
        if (pt.pch == ".") {
            pt.cex <- 3
        }
        else {
            pt.cex <- 1
        }
    }
    wcol0 <- rep(1, length(chrom0))
    if (altcol) {
        j <- 0
        for (i in uchrom) {
            j <- (j + 1)%%2
            wcol0[chrom0 == i] <- 1 + j
        }
    }
    if (missing(pt.cols)) 
        pt.cols <- c("black", "green")
    if (missing(segcol)) 
        segcol <- "red"
    if (missing(zlcol)) 
        zlcol <- "grey"
    if (missing(lwd)) 
        lwd <- 3
    
    for (isamp in 1:nsample) {
        genomdat <- xdat[, isamp + 2]
        ina <- which(!is.na(genomdat) & !(abs(genomdat) == 
                                          Inf))
        genomdat <- genomdat[ina]
        wcol <- wcol0[ina]
        chrom <- chrom0[ina]
        if (xmaploc) 
            maploc <- maploc0[ina]
        ii <- cumsum(c(0, xres$num.mark[xres$ID == sampleid[isamp]]))
        mm <- xres$seg.mean[xres$ID == sampleid[isamp]]
        kk <- length(ii)
        zz <- cbind(ii[-kk] + 1, ii[-1])
        if (missing(ylim)) 
            ylim <- range(c(genomdat, -genomdat))
        omm <- order(mm)
        ozz <- zz[omm, ]
        ina <- unlist(apply(ozz, 1, function(ii) ii[1]:ii[2]))
        if(isamp > 1) par(new = TRUE)
        plot(genomdat[ina], pch = pt.pch, cex = pt.cex, 
             main = NULL, ylab = "", ylim = ylim, col = "orange")
        if(isamp == nsample)
            title("All arrays, superimposed")
        if (zeroline) 
            abline(h = 0, col = zlcol, lwd = lwd)
        if (include.means) {
            ii <- cumsum(c(0, xres$num.mark[xres$ID == 
                                            sampleid[isamp]][omm]))
            smm <- mm[omm]
            zz <- cbind(ii[-kk] + 1, ii[-1])
            for (i in 1:(kk - 1)) lines(zz[i, ], rep(smm[i], 
                                                     2), col = segcol, lwd = lwd)
            
        }
    }
    on.exit(if (plot.type == "chrombysample" | plot.type == "samplebychrom") {
        par(op)
    } else {
        if (int.dev & !parask & nsample > 1) par(ask = parask)
    })
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



mpiWave <- function() {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(cluster))
    mpi.remote.exec(library(waveslim))
    mpi.remote.exec(library(ADaCGH2))
}



wave.aCGH <- function(dat, chrom, minDiff) {
## level to use for wavelet decomposition and thresholding
## The 'recommended' level is floor(log2(log(N)+1)), which
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
    class(out) <- c(class(out), "waveCGH")
    return(out)
}



plot.wavelets2 <- function(res, xdata, chrom,
                           arraynum, main = NULL,
                           colors = c("orange", "red", "green", "blue"),
                           pch = 20, ylim = NULL, html = TRUE,
                           geneNames = positions.merge1$name,
                           idtype = idtype, organism = organism) {
                                        #res is the results
                                        # color code for region status

    logr <- xdata[, arraynum]
    segmented <- res$Predicted[, arraynum]
    col <- rep(colors[1],length(logr))
    simplepos <- 1:length(logr)

    nameIm <- main
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
                         width = imwidth, ps = 12)
    }
    plot(logr ~ simplepos, col="orange", ylab = "log ratio",
         xlab ="Chromosome location", axes = FALSE, cex = 0.7, main = main,
         pch = pch, ylim = ylim,)
    box()
    rug(simplepos, ticksize = 0.01)
    axis(2)
    abline(h = 0, lty = 2, col = "blue")
    
    ## Limit between chromosomes
    LimitChr <- tapply(simplepos,
                       chrom, max)
    abline(v=LimitChr, col="grey", lty=2)

    chrom.nums <- unique(chrom)
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)


    ## segments
    lines(segmented ~ simplepos,
         col = "black", lwd = 2, type = "l")

    if(html) {
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
        createIM2(im1, file = paste(nameIm, ".html", sep = ""))
        imClose(im1)
    }

    if(html) { ## here is chromosome specific code
        pixels.point <- 3
        chrheight <- 500
        for(cnum in 1:length(chrom.nums)) {
            cat(" .... doing chromosome ", cnum, "\n")
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
            plot(logr[indexchr] ~ simplepos[indexchr], col="orange", cex = 1,
                 xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE,
                 main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                 pch = pch, ylim = ylim)
            box()
            axis(2)
            abline(h = 0, lty = 2, col = "blue")
            rug(simplepos[indexchr], ticksize = 0.01)
            ## segments
            lines(segmented[indexchr] ~ simplepos[indexchr],
                  col = "black", lwd = 2, type = "l")

            ## The within chromosome map for gene names
            ## in superimposed cases only in first map
            usr2pngCircle <- function(x, y, rr = 2) {
                xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
                r <- abs(xyrc[2, 1] - xyrc[3, 1])
                return(c(xyrc[1, 1], xyrc[1, 2], r))
            }
            
            ccircle <- mapply(usr2pngCircle, simplepos[indexchr],
                              logr[indexchr])
            write(ccircle, file = "pngCoordChr",
                  sep ="\t", ncolumns = 3)
            write(as.character(geneNames[indexchr]), file = "geneNamesChr")
            imClose(im2)
            system(paste(.python.toMap.py,
                         paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         idtype, organism, sep = " "))
        }        ## looping over chromosomes
    } ## if html
}


plot.wavelets3 <- function(res, xdata, chrom, arraynums = 1:numarrays, main = NULL,
                         colors = c("orange", "red", "green", "blue"),
                         pch = 20, ylim =NULL, html = TRUE,
                         geneNames = positions.merge1$name,
                           idtype = idtype, organism = organism) {
    ## For superimposed: only all genome plot with map to chromosomes

    nameIm <- main
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
                         width = imwidth, ps = 12)
    }
    nfig <- 1
    for (arraynum in arraynums) {
        logr <- xdata[, arraynum]
        segmented <- res$Predicted[, arraynum]
        col <- rep(colors[1],length(logr))
        simplepos <- 1:length(logr)
        
        if(nfig == 1) {
            plot(logr ~ simplepos, col="orange", ylab = "log ratio", 
                 xlab ="Chromosome location", axes = FALSE, main = main,
                 pch = "", ylim = ylim)
            box()
            axis(2)
            abline(h = 0, lty = 2, col = "blue")
            rug(simplepos, ticksize = 0.01)
        }
        lines(segmented ~ simplepos,
              col = "black", lwd = 2, type = "l")
    

        if(nfig == 1) {
            ## Limit between chromosomes
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
        }
        nfig <- nfig + 1
        par(new = TRUE)
    }
    createIM2(im1, file = paste(nameIm, ".html", sep = ""))
    imClose(im1)
}

plot.wavelets4 <- function(res, xdata, chrom, arraynums = 1:numarrays, main = NULL,
                           colors = c("orange", "red", "green", "blue"),
                           pch = 20, ylim =NULL, html = TRUE,
                           geneNames = positions.merge1$name,
                           idtype = idtype, organism = organism) {
    ## For superimposed: one plot per chr
    
    nameIm <- main
    
    pixels.point <- 3
    chrheight <- 500
    chrom.nums <- unique(chrom)
    for(cnum in 1:length(chrom.nums)) {
        cat(" .... doing chromosome ", cnum, "\n")
        indexchr <- which(chrom == chrom.nums[cnum])
        chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
        chrwidth <- max(chrwidth, 800)
        im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         height = chrheight, width = chrwidth,
                         ps = 12)
        ## The following seems needed (also inside sw.plot2) for the coords.
        ## of points to work OK
        nfig <- 1
        for(arraynum in arraynums) { ## first, plot the points
            cat(" ........ for points doing arraynum ", arraynum, "\n")
            logr <- xdata[, arraynum]
            simplepos <- 1:length(logr)
            if(nfig == 1) {
                par(xaxs = "i")
                par(mar = c(5, 5, 5, 5))
                par(oma = c(0, 0, 0, 0))
                plot(logr[indexchr] ~ simplepos[indexchr], col="orange", 
                     xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE, cex = 1,
                     main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     pch = "", ylim = ylim)
                box()
                axis(2)
                abline(h = 0, lty = 2, col = "blue")
                rug(simplepos[indexchr], ticksize = 0.01)
            }
            points(logr[indexchr] ~ simplepos[indexchr], col="orange",
                   cex = 1, pch = 20)
            nfig <- nfig + 1
        }
        for(arraynum in arraynums) { ## now, do the segments
            cat(" ........ for segments doing arraynum ", arraynum, "\n")
            segmented <- res$Predicted[indexchr, arraynum]
            lines(segmented ~ simplepos[indexchr], col = "black", lwd = 2, type = "l")
        }

        usr2pngCircle <- function(x, y, rr = 2) {
            xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
            r <- abs(xyrc[2, 1] - xyrc[3, 1])
            return(c(xyrc[1, 1], xyrc[1, 2], r))
        }
        ccircle <- mapply(usr2pngCircle, simplepos[indexchr],
                          0)
        write(ccircle, file = "pngCoordChr",
              sep ="\t", ncolumns = 3)
        write(as.character(geneNames[indexchr]), file = "geneNamesChr")
        imClose(im2)
        system(paste(.python.toMap.py,
                     paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     idtype, organism, sep = " "))
    }
}

plot.wavelets <- function(res, xdata, chrom,
                          arraynum, main = NULL,
                          colors = c("orange", "red", "green", "blue"),
                          pch = 20, ylim = NULL) {
                                        #res is the results
                                        # color code for region status

    logr <- xdata[, arraynum]
    segmented <- res$Predicted[, arraynum]
    col <- rep(colors[1],length(logr))
    simplepos <- 1:length(logr)
    
    plot(logr ~ simplepos, col="orange", ylab = "log ratio",
         xlab ="", axes = FALSE, cex = 0.7, main = main,
         pch = pch, ylim = ylim)
    box()
    axis(2)
    abline(h = 0, lty = 2, col = "blue")
    
    ## Limit between chromosomes
    LimitChr <- tapply(simplepos,
                       chrom, max)
    abline(v=LimitChr, col="grey", lty=2)

    chrom.nums <- unique(chrom)
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)


    ## segments
    lines(segmented ~ simplepos,
         col = "black", lwd = 2, type = "l")
}



print.wavelets.results <- function(res, xcenter, commondata, output =
                                   "Wavelets.results.txt"){
    ## This function "stretches out" the output and creates a table
    ## that matches the original names, etc.

    out <- data.frame(ID = commondata$name,
                      Chromosome = commondata$chromosome,
                      Start = commondata$start,
                      End = commondata$end,
                      MidPoint = commondata$MidPoint)

    for(i in 1:ncol(xcenter)) {
        t2 <- xcenter[, i]
        t1 <- res$Predicted[, i]
        t3 <- res$State[, i]
        out <- cbind(out, t2, t1, t3)
    }
    colnames(out)[6:(ncol(out))] <-
        paste(rep(colnames(xcenter),rep(3, ncol(xcenter))),
              c(".Original", ".Smoothed", ".State"), sep = "")
    write.table(out, file = output,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
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


mpiPSW <- function() {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(cgh))
    mpi.remote.exec(library(ADaCGH2))
}


## I modify a few things from original sw.plot

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
    mtext(xlab, side = 1, line = 3, cex = 0.8)
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










######## This is the one to use!!!!
my.sw2 <- function(logratio, chrom, sign = -1, p.crit = PSW.p.crit,
                   main = NULL,
                   nIter = 1000,
                   prec = 100,
                   name,
                   highest = FALSE, ## identifying highest scoring island can
                   ## cross chromosome boundaries
                   ...) {

    ## the thresholding is common to all genome,
    ## the anal. is by chromosome
    
    ## all parameters as in the corresponding sw functions.
    ## except p.crit. p.crit is the largest p-value for
    ## which we want a region to be shown, in red,
    ## in the plot.

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
        swt.perm <- sw.perm.test(x, max.nIslands = NULL,
                                        nIter = nIter)
        swt.rob <- sw.rob(x, prec = prec)
        list(swt.run = swt.run, swt.perm = swt.perm, swt.rob = swt.rob)
    }

    papout <- papply(swtlist, funsw,
                     papply_commondata =list(nIter = nIter,
                     prec = prec))

    swt.rob <- unlist(lapply(papout, function(x) x$swt.rob))
    swt.perm <- unlist(lapply(papout, function(x) x$swt.perm))
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
    
##     pdf(file = paste(name, ".pdf", sep = ""),
##         width = png.width,
##         height = png.height)

    sw.plot2(logratio, sign = sign, rob = swt.rob, main = main,
            ...)
#    axis(2)
#    box()
#    axis(3)
#    axis(4)
    ## the kludge: plotting "significant" segments
    sign.segments <- which(swt.perm < p.crit)

    if(sign == 1) {
        red.pos <- quantile(logratio[logratio > 0], p = 0.66)
    } else if(sign == -1) {
        red.pos <- quantile(logratio[logratio < 0], p = 0.33)
    }   
    if(length(sign.segments)) {
        for(i in 1:length(sign.segments)) {
            ## could index by i; but this prevents mistakes
            ## if code changes and signif. no longer in order.
            x0 <- swt.run$start[sign.segments[i]] - 0.5
            x1 <- x0 + swt.run$length[sign.segments[i]]
            y0 <- red.pos
            y1 <- red.pos
            segments(x0, y0, x1, y1, col = "red", lwd = 2)
        }
    }

    ## Limits between chromosomes
    LimitChr <- tapply(1:length(logratio), chrom, max)
    abline(v=LimitChr, col="grey", lty=2)
    
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)
##     dev.off()

    out.values <- cbind(rep(sign, length(logratio)),
                        swt.rob, perm.p.values)
    colnames(out.values) <-
        paste(name, c(".Sign", ".Robust", ".p.value"),
              sep = "")
    return(out.values)

}     


my.sw3 <- function(logratio, chrom, sign = -1, p.crit = PSW.p.crit,
                   main = NULL,
                   nIter = 1000,
                   prec = 100,
                   name,
                   highest = FALSE, ## identifying highest scoring island can
                   ## cross chromosome boundaries
                   ...) {

    ## like my.sw2 but take plotting outside and redefine return output

    
    ## the thresholding is common to all genome,
    ## the anal. is by chromosome
    
    ## all parameters as in the corresponding sw functions.
    ## except p.crit. p.crit is the largest p-value for
    ## which we want a region to be shown, in red,
    ## in the plot.

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
        swt.perm <- sw.perm.test(x, max.nIslands = NULL,
                                        nIter = nIter)
        swt.rob <- sw.rob(x, prec = prec)
        list(swt.run = swt.run, swt.perm = swt.perm, swt.rob = swt.rob)
    }

    papout <- papply(swtlist, funsw,
                     papply_commondata =list(nIter = nIter,
                     prec = prec))

    swt.rob <- unlist(lapply(papout, function(x) x$swt.rob))
    swt.perm <- unlist(lapply(papout, function(x) x$swt.perm))
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
                    

    
    out.values <- cbind(rep(sign, length(logratio)),
                        swt.rob, perm.p.values)
    colnames(out.values) <-
        paste(name, c(".Sign", ".Robust", ".p.value"),
              sep = "")

    out <- list(out=out.values,
                plotdat = plotdat)
    class(out) <- c(class(out), "CGH.PSW")
    return(out)
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
            cat("\n .... doing chromosome ", cnum)
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
            write(ccircle, file = "pngCoordChr",
                  sep ="\t", ncolumns = 3)
            write(as.character(geneNames[indexchr]), file = "geneNamesChr")
            imClose(im2)
            ## call the Python function
            system(paste(.python.toMap.py,
                         paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         idtype, organism, sep = " "))
            

        } ## looping over chromosomes
    } ## if html
}



##             usr2pngRect <- function(x, y, side = 2) {
##                 xleft <- x - side
##                 xright <- x + side
##                 ytop <- y + side
##                 ybottom <- y - side
##                 xycorns <- usr2png(rbind(c(xleft,ytop),c(xright,ybottom)),im2)
##                 return(c(xycorns[1,1],xycorns[2,2],xycorns[2,1],xycorns[1,2]))
##             }
##             write(mapply(usr2pngCircle, location[indexchr],
##                                 logratio[indexchr]),
##                   file = "pngCoordChr",
##                   sep ="\t", ncolumns = 3)
##             write(mapply(usr2pngRect, location[indexchr],
##                                 logratio[indexchr]),
##                   file = "pngCoordChrRect",
##                   sep ="\t", ncolumns = 4)

##             usr2png <- function (xy, im) {
##                  xy <- usr2dev(xy, dev.cur())
##                  cbind(ceiling(xy[, 1] * im$Width), ceiling((1 - xy[, 2]) * 
##                                                         im$Height))
##              }







            
            
##             if ((idtype == "None") | (organism == "None")) 
##                 onclf <- function(num, iii) paste("fixedtooltip(", num,
##                                                   ", '<a class=\'tip\' href=\'javascript:hidetip(1)\'>X</a>",
##                                                   "<a class=\'tip2\'", geneNames[iii],
##                                                   ", this, event, '90px'); return false",
##                                                   sep = "")
##             else
##                 onclf <- function(num, iii)
##                     paste("fixedtooltip(", num,
##                           ", '<a class=\'tip\' href=\'javascript:hidetip(1)\'>X</a>",
##                           "<a class=\'tip2\' href=\'http://idclight.bioinfo.cnio.es/idclight.prog?",
##                           geneNames[iii], "&idtype=",idtype, "&org=", organism,
##                           "\'>", geneNames[iii], "</a>', this, event, '90px'); return false",
##                           sep = "")
                
##             oncl <- mapply(onclf, pnum, idc)
##             f2 <- function(i1, click)
##                 imCircle(location[i1], logratio[i1],
##                          r = rr, title = geneNames[i1],
##                          onClick = click)
##             circlist <- mapply(f2, idc, oncl, SIMPLIFY=FALSE)

##             faddreg <- function(im, value) {
##                 im$HTML[[length(im$HTML)+1]] <- toHTML(value,im)
##                 im$Objects[[length(im$Objects)+1]] <- value
##                 return(im)
##             }
##             mapply(faddreg, im2, circlist)
            
##             for(lll in 1:length(circlist))
##                 addRegion(im2) <- circlist[[lll]]

            ## divs
##             f1 <- function(x)
##                    paste("document.write('<div id=\"",x,
##                     "\" class=\"fixedtipdiv\" style=\"visibility:hidden;width:'",
##                     "+tipwidth+';background-color:'+tipbgcolor+'\" ></div>')", sep ="")
##             write(mapply(f1, pnum), file = "divs")
            
##             ## put together hmtl, divs, and javascript initial stuff
##             system(paste("cat js.adacgh.squeleton1.html divs js.adacgh.squeleton2.html ",
##                    paste("Chr", chrom.nums[cnum], "@", nameIm, ".html", sep =""),
##                    " final.adacgh.squeleton.html > tmphtml; mv tmphtml ",
##                    paste("Chr", chrom.nums[cnum], "@", nameIm, ".html", sep =""), sep =""))




## Old code            
##             browser()
##             cat("\n length(indexchr) ", length(indexchr))
##             for(pnum in 1:length(indexchr)) {
##                 cat(" \n .... doing pnum", pnum)
##                 if ((idtype == "None") | (organism == "None")) {
##                     oncl <- paste("fixedtooltip(", pnum,
##                                   ", '<a class=\'tip\' href=\'javascript:hidetip(1)\'>X</a>",
##                                   "<a class=\'tip2\'", geneNames[indexchr[pnum]],
##                                   ", this, event, '90px'); return false",
##                                   sep = "")
##                     addRegion(im2) <- imCircle(location[indexchr[pnum]], logratio[indexchr[pnum]],
##                                                r = pixels.point/2, title = geneNames[indexchr[pnum]],
##                                                alt = geneNames[indexchr[pnum]], 
##                                                onClick = oncl)
##                 } else {
##                     oncl <- paste("fixedtooltip(", pnum,
##                                   ", '<a class=\'tip\' href=\'javascript:hidetip(1)\'>X</a>",
##                                   "<a class=\'tip2\' href=\'http://idclight.bioinfo.cnio.es/idclight.prog?",
##                                   geneNames[indexchr[pnum]], "&idtype=",idtype, "&org=", organism,
##                                   "\'>", geneNames[indexchr[pnum]], "</a>', this, event, '90px'); return false",
##                                   sep = "")
##                     addRegion(im2) <- imCircle(location[indexchr[pnum]], logratio[indexchr[pnum]],
##                                                r = pixels.point/2, title = geneNames[indexchr[pnum]],
##                                                alt = geneNames[indexchr[pnum]], 
##                                                onClick = oncl)
##                 }
##             } ## /the javascript calls per point
##             cat("\n yyyyy por aquí\n")
##             createIM(im2, file = paste("Chr", chrom.nums[cnum], "@", nameIm, ".html", sep =""))
##             imClose(im2)

          
sw.plot4 <- function (logratio, location = seq(length(logratio)),
                      threshold.func = function(x) median(x) + 
                      0.2 * mad(x), sign = -1, highest = TRUE, expected = NULL, 
                      rob = NULL, legend = TRUE, xlab = "Chromosomal location", 
                      ylab = "log ratio",
                      swt.run,
                      swt.perm,
                      p.crit,
                      chrom, main,
                      html = TRUE,
                      nameIm, geneNames,
                      idtype = idtype, organism = organism) {

    ## this puts the gene names
    
    ## geneNames often = positions.merge1$name

    if(html) {
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
        for(i in 1:length(sign.segments)) {
            ## could index by i; but this prevents mistakes
            ## if code changes and signif. no longer in order.
            x0 <- swt.run$start[sign.segments[i]] - 0.5
            x1 <- x0 + swt.run$length[sign.segments[i]]
            y0 <- red.pos
            y1 <- red.pos
            segments(x0, y0, x1, y1, col = "red", lwd = 2)
        }
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
        for(nd in 1:length(logratio)) {
            addRegion(im1) <- imCircle(location[nd], logratio[nd],
                                       r = 5, title = geneNames[nd],
                                       alt = geneNames[nd], 
                                       onClick = paste("_onclick_(", geneNames[nd],
                                       "*", idtype,
                                       "*", organism, ")", sep = ""))
        }
        createIM2(im1, file = paste(nameIm, ".html", sep = ""))
        imClose(im1)
    }

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



mpiACE <- function() {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(ADaCGH2))
}



#data("file.aux.RData")


ace.analysis.C <- function(x, coefs=file.aux, Sdev, array.names="x")
{
  #Coefficients of file
##  coefs <- read.table(file=file.aux, header=TRUE, sep="\t")

  if(is.null(array.names)) array.names <- "x"
  
  #Obtain the null distribution of the (L,H)-pairs
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
                
  ##Clean up the variables and add 1 to first and last
  ##for the difference in array enumeration in C and R
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
    
                                        #Step 1		
                                        #Segmentation
                                        #Compute running mean
    
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
    
                                        #Find change points
    
                                        #nknots <- 0
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
    
    
                                        #Step 2
                                        #Feature extraction
                                        #Compute size, height, first, last tuples
    
    Nclusters <- length(knots) - 1
    first <- knots[1:Nclusters] + 1
    last <- knots[2:(Nclusters+1)] 
    size <- last - first + 1
    height <- mapply(function(first, last) mean(x[first:last]), first=first, last=last)
    z1 <- size
    z2 <- abs(height)
    
    
	
                                        #Step 3
                                        #Obtain the null distribution of the (L,H)-pairs
    
                                        #Step 6
                                        #Report genes
                                        #Adjust start/end positions of clusters
                                        #Could be optimized
    for (j in 1:Nclusters) {
        fi <- first[j]
        la <- last[j]
        err.opt <- 999999999
        p.opt <- fi
        q.opt <- la
        r <- floor(max(0, min(16, (la-fi)/2)-1))
        for (p in fi:(fi+r)) {
            for(q in la:(la-r)) {
                                        #Hay que vigilar que se cumplan las condiciones de los for
                e1 <- ifelse(fi<p,sum(x[fi:(p-1)]^2),0)
                                        #Hay que vigilar que se cumplan las condiciones de los for
                e2 <- ifelse(p<=q,sum((x[p:q]-height[j])^2),0)
                                        #Hay que vigilar que se cumplan las condiciones de los for
                e3 <- ifelse(q+1<=la,sum(x[(q+1):la]^2),0)
                err <- (e1 + e2 + e3) /(la-fi)
                ##In case that (la-fi)==0
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
                                        #Step 4
                                        #Find significant genes
                                        #Determine which clusters are inside/outside
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
                                        #number of rejections
    called <- !(v1<=0 & v2>=0)
    
                                        #Step 5
                                        #Estimate the positive false discovery rate
    
    
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
		#Estimation of the variance parameter: with all samples and chromosomes, 
		#so we save the number of genes to make weighted mean
		#Indexes have changed from 21 to 41 in this version
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
 	invisible(res)
 }

## ACE <- function(x, Chrom, coefs = file.aux, Sdev=0.2) {
## 	####### x is log2.ratio
## 	####### Chrom ---MUST BE NUMERIC-- 
## 	###### only for 1 array at a time
## 	if (!is.numeric(Chrom)) {
## 		stop("Chromosome variable must be numeric")
## 	}

## 	array.names <- colnames(x)
		
## 	nchrom <- length(unique(Chrom))
## 	if(is.null(dim(x))) {
## 		genes <- split(x, Chrom)
## 		first.estimate <- lapply(genes, ace.analysis.C, Sdev=Sdev, array.names=array.names)
## 		Sdevs.estimate <- lapply(first.estimate, sd.ACE.analysis)
## 		Sdevs <- unlist(lapply(Sdevs.estimate, "[", 1))
## 		Sdevs <- mean(Sdevs[Sdevs>0])
## 		res <- lapply(genes, ace.analysis.C, Sdev=Sdevs, array.names=array.names)
## 		class(res) <- "ACE"
## 		}
## 	else {
## 		Sdevs <- matrix(NA, ncol(x),nchrom)
## 		res <- list()
## 		for(i in 1:ncol(x)) {
			
## 			genes <- split(x[,i], Chrom)
## 			first.estimate <- lapply(genes, ace.analysis.C, Sdev=Sdev, array.names=array.names[i])
## 			Sdevs.estimate <- lapply(first.estimate, sd.ACE.analysis)
## 			Sdevs[i,] <- unlist(lapply(Sdevs.estimate, "[", 1))
## 		}
## 		Sdevs <- mean(Sdevs[Sdevs>0])
## 		for (i in 1:ncol(x)) {
## 			genes <- split(x[,i], Chrom)
## 			res[[i]] <- lapply(genes, ace.analysis.C, Sdev=Sdevs, array.names=array.names[i])	
## 			class(res[[i]]) <- "ACE"
## 		}
## 		class(res) <- "ACE.array"
## 	}
## 	invisible(res)
## }

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

summary.ACE <- function(object, fdr=NULL, html = TRUE,
                        outhtml ="ace.fdrtable.html") {

	nchrom <- length(object)
	FDR.table <- get.FDR(object, nchrom)

	if (is.null(fdr)) fdr <- 0.15

        index <- which.min(abs(FDR.table[, 3] - fdr))
	print(FDR.table)
        print(paste("Selected index", index))
        cat(FDR.table[index, 3], file ="aceFDR")

        aceFDR.for.output <- FDR.table[index, 3]
        
        if(html) {
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
        ## sacar esto como html y añadir, en output,
        ## la form de cgi, etc, para seleccionar el fdr.
        ## ese cgi genera los nuevos dibujos, pero no toca lo de antes,
        ## y pone todo en una sola página html (genera de nuevo el results. html)


        ## we need to call, in the miniACE.R, this (or the array equiv) and
        ## the plotting and general table results.


        
## 	####Select a FDR level
## 	if (is.null(index)) {
## 		print(FDR.table)
## 		index <- readline(paste("Please select an index (1 -", nrow(FDR.table), ") for FDR level:"))
## 		if (!index %in% 1:nrow(FDR.table)) {
## 			stop("Not a valid index")
## 			}
## 		index <- as.numeric(index)
## 		}
		
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
	class(res) <- c("summary.ACE", "CGH.ACE.summary")
	res <- do.call("rbind", res)
	rownames(res) <- 1:nrow(res)
        ncr <- ncol(res) - 1
        attr(res, "aceFDR.for.output") <- aceFDR.for.output
	res
}




print.summary.ACE <- function(obj) {
	class(obj) <- "data.frame"
	rownmames(obj) <- 1:nrow(obj)
	print(obj)
	invisible(obj)
}


summary.ACE.array <- function(object, fdr=NULL, html = TRUE,
                              outhtml ="ace.fdrtable.html") {

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
        
        if(html) {
            ## Place this at bottom of page
            ## what gets change is this itself and the figures
            cat("<TABLE border>\n", file = outhtml, append=FALSE)
            cat("<tr><td>Index</td><td>Number of Genes with gains/losses</td><td>FDR</td></tr>\n",
                file = outhtml, append = TRUE)
            for(ni in 1:dim(FDR.table)[1]) {
                if(ni == index) 
                    cat(paste("<tr><td><b>", FDR.table[ni, 1], "</b></td><td><b>",
                              FDR.table[ni, 2], "</b></td><td><b>", FDR.table[ni, 3],
                              "</b></td></tr></b>\n"), file = outhtml, append = TRUE)
                else
                    cat(paste("<tr><td>", FDR.table[ni, 1], "</td><td>",
                              FDR.table[ni, 2], "</td><td>", FDR.table[ni, 3],
                              "</td></tr>\n"), file = outhtml, append = TRUE)
            }
            cat("</TABLE>", file = outhtml, append= TRUE)
        }


        
	####Select a FDR level
## 	if (is.null(index)) {
## 		print(FDR.table)
## 		index <- readline(paste("Please select an index (1 -", length(FDR), ") for FDR level:"))
## 		if (!index %in% 1:length(FDR)) {
## 			stop("Not a valid index")
## 			}
## 		index <- as.numeric(index)
## 		}
		
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
		names(res[[i]])[2] <- array.names[i]
                names(res[[i]])[3] <- paste("Gain.Loss.", array.names[i], sep = "")
		}
		class(res) <- c("summary.ACE.array", "CGH.ACE.summary")
        attr(res, "aceFDR.for.output") <- aceFDR.for.output
		res
}

print.summary.ACE.array <- function(obj) {
	for (i in 1:length(obj)) {
		class(obj[[i]]) <- "data.frame"
		}
	obj <- do.call("cbind", obj)
	#### Delete redundant info: chromosomes from 2nd array on
	obj <- obj[,-seq(from=4, to=ncol(obj), by=3)]
	rownames(obj) <- 1:nrow(obj)
	print(obj)
	invisible(obj)
}

plot.ace <- function(res, chrom,
                          arraynum, main = NULL,
                          colors = c("orange", "red", "green", "blue"),
                          pch = 20, ylim = NULL) {
                                        #res is the results
                                        # color code for region status

    logr <- res[[arraynum]][, 2]
    res.dat <- res[[arraynum]][, 3]

    col <- rep("orange",length(res.dat))
    col[which(res.dat == -1)] <- "green"
    col[which(res.dat == 1)] <- "red"

    simplepos <- 1:length(res.dat)
    plot(logr ~ simplepos, col= col, ylab = "log ratio",
         xlab ="", axes = FALSE, cex = 0.7, main = main,
         pch = pch, ylim = ylim)
    box()
    axis(2)
    abline(h = 0, lty = 2, col = "blue")
    
    ## Limit between chromosomes
    LimitChr <- tapply(simplepos,
                       chrom, max)
    abline(v=LimitChr, col="grey", lty=2)
    chrom.nums <- unique(chrom)
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)

    lines(I(res.dat/2) ~ simplepos, col="black",
           lwd = 2)
}







plot.ace2 <- function(res, chrom,
                      arraynum, main = NULL,
                      colors = c("orange", "red", "green", "blue"),
                      pch = 20, ylim = NULL,
                      geneNames = positions.merge1$name,
                      html = TRUE,
                      idtype = idtype, organism = organism) {
                                        #res is the results
                                        # color code for region status

    logr <- res[[arraynum]][, 2]
    res.dat <- res[[arraynum]][, 3]
    col <- rep("orange",length(res.dat))
    col[which(res.dat == -1)] <- "green"
    col[which(res.dat == 1)] <- "red"
    simplepos <- 1:length(res.dat)

    nameIm <- main
    cat(" nameIm is ", nameIm)
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
                         width = imwidth, ps = 12)
    }
##    browser()
    plot(logr ~ simplepos, col= col, ylab = "log ratio",
         xlab ="Chromosome location", axes = FALSE, cex = 0.7, main = main,
         pch = pch, ylim = ylim)
    box()
    rug(simplepos, ticksize = 0.01)
    axis(2)
    abline(h = 0, lty = 2, col = "blue")
    
    ## Limit between chromosomes
    LimitChr <- tapply(simplepos,
                       chrom, max)
    abline(v=LimitChr, col="grey", lty=2)
    chrom.nums <- unique(chrom)
    d1 <- diff(LimitChr)
    pos.labels <- c(round(LimitChr[1]/2),
                    LimitChr[-length(LimitChr)] + round(d1/2))
    axis(1, at = pos.labels, labels = chrom.nums)

    lines(I(res.dat/2) ~ simplepos, col="black",
           lwd = 2)

    if(html) {
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
        createIM2(im1, file = paste(nameIm, ".html", sep = ""))
        imClose(im1)
    }

    if(html) { ## here is chromosome specific code
        pixels.point <- 3
        chrheight <- 500
        for(cnum in 1:length(chrom.nums)) {
            cat(" .... doing chromosome ", cnum, "\n")
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
            plot(logr[indexchr] ~ simplepos[indexchr], col=col[indexchr], cex = 1,
                 xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE,
                 main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                 pch = pch, ylim = ylim)
            box()
            axis(2)
            abline(h = 0, lty = 2, col = "blue")
            rug(simplepos[indexchr], ticksize = 0.01)
            ## segments
            lines(I(res.dat[indexchr]/2) ~ simplepos[indexchr],
                  col = "black", lwd = 2, type = "l")

            ## The within chromosome map for gene names
            ## in superimposed cases only in first map
            usr2pngCircle <- function(x, y, rr = 2) {
                xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
                r <- abs(xyrc[2, 1] - xyrc[3, 1])
                return(c(xyrc[1, 1], xyrc[1, 2], r))
            }
            
            ccircle <- mapply(usr2pngCircle, simplepos[indexchr],
                              logr[indexchr])
            write(ccircle, file = "pngCoordChr",
                  sep ="\t", ncolumns = 3)
            write(as.character(geneNames[indexchr]), file = "geneNamesChr")
            imClose(im2)
            system(paste(.python.toMap.py,
                         paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         idtype, organism, sep = " "))
        }        ## looping over chromosomes
    } ## if html
}

plot.ace3 <- function(res, chrom, arraynums = 1:numarrays, main = NULL,
                         colors = c("orange", "red", "green", "blue"),
                         pch = "", ylim =c(ymin, ymax), html = TRUE,
                         geneNames = positions.merge1$name,
                      idtype = idtype, organism = organism) {
    ## For superimposed: only all genome plot with map to chromosomes

    nameIm <- main
    if(html) {
        imheight <- 500
        imwidth <- 1600
        im1 <- imagemap3(nameIm, height = imheight,
                         width = imwidth, ps = 12)
    }
    nfig <- 1
    for (arraynum in arraynums) {
        logr <- res[[arraynum]][, 2]
        res.dat <- res[[arraynum]][, 3]
        segmented <- res.dat/2
        col <- rep("orange",length(res.dat))
        col[which(res.dat == -1)] <- "green"
        col[which(res.dat == 1)] <- "red"
        simplepos <- 1:length(logr)
        
        if(nfig == 1) {
            plot(logr ~ simplepos, col=col, ylab = "log ratio", 
                 xlab ="Chromosome location", axes = FALSE, main = main,
                 pch = "", ylim = ylim)
            box()
            axis(2)
            abline(h = 0, lty = 2, col = "blue")
            rug(simplepos, ticksize = 0.01)
        }
        lines(segmented ~ simplepos,
              col = "black", lwd = 2, type = "l")
    

        if(nfig == 1) {
            ## Limit between chromosomes
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
        }
        nfig <- nfig + 1
        par(new = TRUE)
    }
    createIM2(im1, file = paste(nameIm, ".html", sep = ""))
    imClose(im1)
}


plot.ace4 <- function(res, chrom, arraynums = 1:numarrays,
                           main = "All_arrays",
                           colors = c("orange", "red", "green", "blue"),
                           pch = 20, ylim =NULL, html = TRUE,
                           geneNames = positions.merge1$name,
                      idtype = idtype, organism = organism) {
    ## For superimposed: one plot per chr
    
    nameIm <- main
    
    pixels.point <- 3
    chrheight <- 500
    chrom.nums <- unique(chrom)
    for(cnum in 1:length(chrom.nums)) {
        cat(" .... doing chromosome ", cnum, "\n")
        indexchr <- which(chrom == chrom.nums[cnum])
        chrwidth <- round(pixels.point * (length(indexchr) + .10 * length(indexchr)))
        chrwidth <- max(chrwidth, 800)
        im2 <- imagemap3(paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                         height = chrheight, width = chrwidth,
                         ps = 12)
        ## The following seems needed (also inside sw.plot2) for the coords.
        ## of points to work OK
        nfig <- 1
        for(arraynum in arraynums) { ## first, plot the points
            cat(" ........ for points doing arraynum ", arraynum, "\n")
            logr <- res[[arraynum]][, 2]
            res.dat <- res[[arraynum]][, 3]
            col <- rep("orange",length(res.dat))
            col[which(res.dat == -1)] <- "green"
            col[which(res.dat == 1)] <- "red"
            simplepos <- 1:length(logr)
            if(nfig == 1) {
                par(xaxs = "i")
                par(mar = c(5, 5, 5, 5))
                par(oma = c(0, 0, 0, 0))
                plot(logr[indexchr] ~ simplepos[indexchr], col=col[indexchr], 
                     xlab ="Chromosomal location", ylab = "log ratio", axes = FALSE, cex = 1,
                     main = paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     pch = "", ylim = ylim)
                box()
                axis(2)
                abline(h = 0, lty = 2, col = "blue")
                rug(simplepos[indexchr], ticksize = 0.01)
            }
            points(logr[indexchr] ~ simplepos[indexchr], col=col[indexchr],
                   cex = 1, pch = 20)
            nfig <- nfig + 1
        }
        for(arraynum in arraynums) { ## now, do the segments
            cat(" ........ for segments doing arraynum ", arraynum, "\n")
            res.dat <- res[[arraynum]][indexchr, 3]
            segmented <- res.dat/2
            lines(segmented ~ simplepos[indexchr], col = "black", lwd = 2, type = "l")
        }

        usr2pngCircle <- function(x, y, rr = 2) {
            xyrc <- usr2png(cbind(c(x, rr, 0), c(y, 0, 0)), im2)
            r <- abs(xyrc[2, 1] - xyrc[3, 1])
            return(c(xyrc[1, 1], xyrc[1, 2], r))
        }
        ccircle <- mapply(usr2pngCircle, simplepos[indexchr],
                          0)
        write(ccircle, file = "pngCoordChr",
              sep ="\t", ncolumns = 3)
        write(as.character(geneNames[indexchr]), file = "geneNamesChr")
        imClose(im2)
        system(paste(.python.toMap.py,
                     paste("Chr", chrom.nums[cnum], "@", nameIm, sep =""),
                     idtype, organism, sep = " "))
    }
}


print.ACE.results <- function(res,
                              output = paste("ACE.results.FDR=",
                              attr(res, "aceFDR.for.output"), ".txt", sep ="")){
    ## This function "stretches out" the output and creates a table
    ## that matches the original names, etc.

    out <- data.frame(ID = positions.merge1$name,
                      Chromosome = positions.merge1$chromosome,
                      Start = positions.merge1$start,
                      End = positions.merge1$end,
                      MidPoint = positions.merge1$MidPoint)

    if(!is.null(dim(xcenter))) {
        for(i in 1: length(res)) {
            outtmp <- cbind(res[[i]][, 2],
                            res[[i]][, 3])
            colnames(outtmp) <- paste(names(res[[i]])[2],
                                      c(".Original", ".State"),
                                      sep = "")
            out <- cbind(out, outtmp)
        }
    } else {
        outtmp <- cbind(res[, 2],
                        res[, 3])
        colnames(outtmp) <- paste(names(res)[2],
                                  c(".Original", ".State"),
                                  sep = "")
        out <- outtmp
    }
        
    write.table(out, file = output,
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
}



##### cghMCR examples

## library(cghMCR)
## data("sampleData")
## segments <- getSegments(sampleData)

## names(segments)
## segments2 <- segments
## segments2$data <- ""

## mc1 <- cghMCR(segments)
## mc2 <- cghMCR(segments2)



## dim(segments[[1]])
## ##
## > segments[[1]][1:10,1:4]
##             probe chrom  maploc     array1
## 960  A_14_P136500     1  799681      .23


## segments[[2]][1:10,]
##                                                           ID chrom loc.start
## 1  X.home.john.lib64.R.library.cghMCR.sampledata.sample1.agi     1    799681
## 2  X.home.john.lib64.R.library.cghMCR.sampledata.sample1.agi    10   1036908
##      loc.end num.mark seg.mean
## 1  244806850      424   0.1789
## 2 





###############################3  Parallelization notes

### DNAcopy: look into segment, the loop over nsample.

### PSW: we already loop over samples... but my.sw2 loops over chromos.
### Use both?
### each sample is assigned as many CPUs as chromosomes, and then we do simultaneously as
### many samples as possible.


### ACE:
##  we loop over samples (inside ACE); trivial to use a papply instead of lapply
##     inside ace.analysis there are several loops, but none is trivial to parallelize.


#### Serial code

#### wave.aCGH <- function(dat, chrom, minDiff) {
#### ## level to use for wavelet decomposition and thresholding
#### ## The 'recommended' level is floor(log2(log(N)+1)), which
#### ## equals 3 for:  21 <= N <= 1096
####     thrLvl <- 3

####     ncloneschrom <- tapply(dat[, 1], chrom, function(x) length(x))
####     if((max(ncloneschrom) > 1096) | (min(ncloneschrom) < 21))
####         warningsForUsers <-
####             c(warningsForUsers,
####               paste("The number of clones/genes is either",
####                     "larger than 1096 or smaller than 21",
####                     "in at least one chromosome. The wavelet",
####                     "thresholding of 3 might not be appropriate."))
        
    
####     Nsamps  <- ncol(dat)
####     Nclones <- nrow(dat)
    
####     state <- matrix(NA, nrow=Nclones, ncol=Nsamps)
####     pred  <- matrix(NA, nrow=Nclones, ncol=Nsamps)

####     uniq.chrom <- unique(chrom)
####     for (i in 1:Nsamps)
####     {
####         ratio.i <- dat[,i]
####         noNA  <- !is.na(ratio.i)
        
####         for (j in uniq.chrom)
####         {
####             chr.j <- (chrom == j)
####             use.ij <- which(noNA & chr.j)
            
####             ratio.ij <- ratio.i[use.ij]
            
####             ## do wavelet deconstruction, threshold, then reconstruct
####             wc   <- modwt(ratio.ij, "haar", n.levels=thrLvl)
            
####             ## These are the three different thresholding functions used
####             ##thH  <- our.sure(wc, max.level=thrLvl, hard=FALSE)
####             thH  <- our.hybrid(wc, max.level=thrLvl, hard=FALSE)
####             ##thH  <- nominal.thresh(wc, max.level=thrLvl, hard=FALSE, sig=.05)
            
####             ## reconstruct the thresheld ('denoised') data
####             recH <- imodwt(thH)
            
####             ## Categorize the denoised data then combine ("merge") levels that
####             ## have predicted values with an absolute difference < 'minDiff' 
####             pred.ij <- segmentW(ratio.ij, recH, minDiff=minDiff)
            
####             pred[use.ij, i]  <- pred.ij
            
####             labs <- as.character(1:length(unique(pred.ij)))
####             state[use.ij, i] <- as.integer(factor(pred.ij, labels=labs))
####         }
####     }
####     out <- list(Predicted =pred, State = state)
#### ##    rownames(out) <- rownames(dat)
####     return(out)
#### }




##### sw.plot2 <- function (logratio, location = seq(length(logratio)),
#####                       threshold.func = function(x) median(x) + 
#####                       0.2 * mad(x), sign = -1, highest = TRUE, expected = NULL, 
#####                       rob = NULL, legend = TRUE, xlab = "Chromosomal location", 
#####                       ylab = "log ratio", ...) 
##### {
#####     my.line <- function(x, y, ...) {
#####         len <- length(x)
#####         run <- rle(y)[[1]]
#####         run.len <- length(run)
#####         j <- 1
#####         m <- 2 * x[1] - x[2]
#####         if (run.len == 1) 
#####             lines(x = c(3/2 * x[1] - 1/2 * x[2], 3/2 * x[len] - 
#####                 1/2 * x[len - 1]), y = c(y[1], y[len]), ...)
#####         else {
#####             for (i in 1:(run.len - 1)) {
#####                 k <- run[i]
#####                 lines(x = c((x[j] + m)/2, (x[j + k - 1] + x[j + 
#####                   k])/2), y = c(y[j], y[j]), ...)
#####                 lines(x = rep((x[j + k - 1] + x[j + k])/2, 2), 
#####                   y = c(y[j], y[j + k]), ...)
#####                 m <- x[j + k - 1]
#####                 j <- j + k
#####             }
#####             lines(x = c((m + x[j])/2, 3/2 * x[len] - 1/2 * x[len - 
#####                 1]), y = c(y[j], y[j]), ...)
#####         }
#####     }
#####     island.line <- function(x, y, start = 1, len = length(x), 
#####         edge = 0, ...) {
#####         if (is.null(start) || is.null(len) || length(start) == 
#####             0 || length(len) == 0 || len <= 0) 
#####             return
#####         lenx <- length(x)
#####         x1 <- c(2 * x[1] - x[2], x, 2 * x[lenx] - x[lenx - 1])
#####         x2 <- x1[1:(lenx + 1)] + diff(x1)/2
#####         lines(x = c(rep(x2[start], 2), rep(x2[start + len], 2)), 
#####             y = c(y - edge, y, y, y - edge), ...)
#####     }
#####     log2 <- log(2)
#####     len <- length(logratio)
#####     par <- par()
#####     par(mar = c(5, 5, 5, 5), yaxs = "r")
#####     par(oma = c(0, 0, 0, 0))
#####     threshold <- threshold.func(sign * logratio)
#####     plot(y = logratio, x = location, type = "n", xlab = "", ylab = "", 
#####         ..., axes = FALSE)
    
#####     rug(location, ticksize = 0.01)
#####     maxlr <- max(logratio)
#####     minlr <- min(logratio)
#####     axis(side = 4, at = seq(minlr, maxlr, length = 5), labels = c("0", 
#####         ".25", " .50", ".75", "1"))
#####     box()
#####     axis(2)
#####     mtext(xlab, side = 1, line = 3, cex = 0.8)
#####     mtext(ylab, side = 2, line = 3, cex = 0.8)
##### ##    abline(h = maxlr, lty = 2)
##### ##    abline(h = minlr, lty = 2)
#####     if (!is.null(rob)) {
#####         mtext("Robustness", side = 4, line = 3, cex = 0.8)
#####         my.line(y = (maxlr - minlr) * rob + minlr, x = location, 
#####             col = "#99ffff", lwd = 2)
#####     }
#####     if (highest) {
#####         x <- sign * logratio - threshold
#####         swx <- sw(x)
#####         if (length(swx$score)) {
#####             island.line(x = location, y = maxlr + (maxlr - minlr) * 
#####                 0.02, start = swx$start[1], len = swx$length[1], 
#####                 edge = (maxlr - minlr) * 0.01, col = "green")
#####         }
#####     }
##### ##    my.line(y = rep(sign * threshold, len), x = location, col = "#00ff00")
#####     if (!is.null(expected)) {
#####         my.line(y = log2(expected) - 1, x = location, col = "#ff0000")
#####     }
##### ##    points(y = logratio, x = location, pch = 20, col = "orange", cex = 0.5)
#####     points(y = logratio, x = location, pch = 20, col = "orange", cex = 1)
#####     if (legend) { ## I make changes here; essentially, I break a few things!
##### ##         legend.str <- c("highest-scoring island", "robustness", "signif.")
##### ##         legend.col <- c("green", "#99ffff", "red")
##### ##         legend(location[1], minlr + (maxlr - minlr) * 0.25, legend.str, 
##### ##             lty = rep(1, 3), col = legend.col, cex = 0.8)
#####         legend.str <- c("robustness", "signif.")
#####         legend.col <- c("#99ffff", "red")
#####         legend(location[1], minlr + (maxlr - minlr) * 0.25, legend.str, 
#####             lty = rep(1, 2), col = legend.col, cex = 0.8)

#####     }
#####     abline(h = 0, lty = 2, col = "blue")
#####     par(mar = par$mar, yaxs = par$yaxs)
##### }
