####  Copyright (C) 2005, 2006, 2007, Ramon Diaz-Uriarte <rdiaz02@gmail.com>

#### This program is free software; you can redistribute it and/or
#### modify it under the terms of the GNU General Public License
#### as published by the Free Software Foundation; either version 2
#### of the License, or (at your option) any later version.

#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### GNU General Public License for more details.

#### You should have received a copy of the GNU General Public License
#### along with this program; if not, write to the Free Software
#### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
#### USA.

######################################################################
######################################################################
######################################################################


### FIXME: qué pasa con un solo array o un solo chrom???





assign(".__ADaCGH_SERVER_APPL", TRUE)


library(ADaCGH21, verbose = FALSE)
cat("\nADaCGH2 Version :\n")
packageDescription("ADaCGH21")$Version
cat("\n\n")





##############################################
##############################################
######                              ##########
######         Error checking       ##########
######          utilities           ##########
######            and               ##########
######          other functions     ##########
######                              ##########
##############################################
##############################################


doCheckpoint <- function(num, to.save, delete.rest = TRUE) {
##    checkpoint.num.new <- num
  if(!is.null(to.save)) {
    save(file = ".RData", list = to.save, envir = .GlobalEnv,
         compress = FALSE)
    if(delete.rest) {
      to.delete <- setdiff(ls(envir = .GlobalEnv), to.save)
      rm(list = to.delete, envir = .GlobalEnv)
    }
  }
##    checkpoint.num <<- num
    sink("checkpoint.num")
    cat(num)
    sink()
    return(num)
}

### We will use hidden stuff from ADaCGH
#doCheckpoint <- ADaCGH:::doCheckpoint


caughtUserError.Web <- ADaCGH21:::caughtUserError.Web
caughtOurError.Web <- ADaCGH21:::caughtOurError.Web

NormalTermination <- function(){
    ADaCGH21:::snowfall.clean.quit.Web()
    status <- file("R_Status.txt", "w")
    cat("Normal termination\n", file = status)
    flush(status)
    close(status)
    cat("\n\n Normal termination\n")
    quit(save = "no", status = 0, runLast = FALSE)
}

readOptions <- function(x) {
    ### return a list, because many more indexing options
    tmp <- read.table(x, sep = "\t", stringsAsFactors = FALSE)
    l1 <- sapply(tmp[, 2], as.list)
    names(l1) <- tmp[, 1]
    return(l1)
}

checkAssign <- function(value, rangeOK, lista) {
    if(is.null(lista[[value]]))
       caughtUserError.Web(paste(value, "has no value"))
    if(lista[[value]] %in% rangeOK)
        return(lista[[value]])
    else
        caughtUserError.Web(paste(value, "has a value we do not accept"))
}

checkConvertMethodOptions <- function(method.options, options) {
  methodaCGH <- options[["method"]]
  indexopts <- which(names(method.options) == methodaCGH)
  if(length(indexopts) > 0) {
    all.method.options <- method.options[[indexopts]]
    for(mo in all.method.options) {
      if(!(mo %in% names(options)))
        caughtUserError.Web(paste(mo, " not among the options."))
      try1 <- try(options[[mo]] <- as.numeric(options[[mo]]))
      if(inherits(try1, "try-error"))
        caughtUserError.Web(paste("User failure for option ", mo))
    }
  }
  return(options)
}


acceptedIDTypes <- c('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl',
                     'entrez', 'ug', 'rsrna', 'rspeptide', 'hugo')
acceptedOrganisms <- c('None', 'Hs', 'Mm', 'Rn')
acceptedMethodaCGH <- c ('Wavelets', 'DNAcopy', 'GLAD', 'HMM', 'BioHMM',
                      'CGHseg', 'HaarSeg')
methodOptions <- list('Wavelets' = c('Wave.minDiff'),
                      'CGHseg'   = c('CGHseg.s'),
                      'HaarSeg'  = c('HaarSeg.m')
                      )
acceptedColors <- colors()





put.part.rdata.together <- function(i, pos.start, pos.end) {

  print(system.time( {
    load("probeNames.RData")
    probeNames <- probeNames[pos.start[i]:pos.end[i]]
  }))


  print(system.time( {
    
    load("chromData.RData")
    open(chromData, readonly = TRUE)
    chrs.part <- chromData[pos.start[i]:pos.end[i]]
    close(chromData)
    rm(chromData)
  }))
        
  print(system.time( {
    
    load("posData.RData")
    open(posData, readonly = TRUE)
    poss.part <- posData[pos.start[i]:pos.end[i]]
    close(posData)
    rm(posData)
  }))



  print(system.time( {
    
    load("segmres.RData")
    open(segmres[[2]], readonly = TRUE)
    oname <- paste("calls.out.", i, sep = "")
    assign(oname,
           data.frame(
                      ProbeName = probeNames,
                      Chr = chrs.part,
                      Position = poss.part,
                      segmres[[2]][pos.start[i]:pos.end[i], , drop = FALSE])
           )
    save(file = paste(oname, ".RData", sep = ""), list = c(oname))
    write.table(file = paste(oname, ".txt", sep = ""), get(oname),
                row.names = FALSE, sep = "\t", quote = FALSE)
    rm(list = c(oname))
    close(segmres[[2]])
  }))

  print(system.time( {
    
    open(segmres[[1]], readonly = TRUE)
    oname <- paste("segmented.out.", i, sep = "")
    assign(oname,
           data.frame(
                      ProbeName = probeNames,
                      Chr = chrs.part,
                      Position = poss.part,
                      segmres[[1]][pos.start[i]:pos.end[i], , drop = FALSE])
           )
    save(file = paste(oname, ".RData", sep = ""), list = c(oname))
    write.table(file = paste(oname, ".txt", sep = ""), get(oname),
                row.names = FALSE, sep = "\t", quote = FALSE)
    rm(list = c(oname))
    close(segmres[[1]])
    rm(segmres)
  }))
  print(system.time( {
    
    load("cghData.RData")
    open(cghData, readonly = TRUE)
    oname <- paste("original.", i, sep = "")
    assign(oname,
           data.frame(
                      ProbeName = probeNames,
                      Chr = chrs.part,
                      Position = poss.part,
                      cghData[pos.start[i]:pos.end[i], , drop = FALSE])
           )
    save(file = paste(oname, ".RData", sep = ""), list = c(oname))
    rm(list = c(oname))
    close(cghData)
    rm(cghData)
  }))
 
}




new.custom2 <- function(segmresRDataName = "segmres.RData",
                       cghRDataName = "cghData.RData",
                       chromRDataName = "chromData.RData",
                       posRDataName = "posData.RData",
                       probeNamesRDataName = "probeNames.RData",
                       nround = 6) {

  load(chromRDataName)
  rle.chr <- intrle(as.integer(chromData[]))
  chr.end <- cumsum(rle.chr$lengths)
  chr.start <- c(1, chr.end[-length(chr.end)] + 1)
  seqc <- seq.int(1, length(chr.start))

  null <- sfClusterApplyLB(seqc, put.part.rdata.together,
                           chr.start, chr.end)

  ## using OS cat all calls and segmented
  os.call.1 <- paste("cat ",
                     paste("calls.out.", seqc, ".txt", sep = "", collapse = " "),
                     " > calls.out.txt")
  system(os.call.1)

  system("head -1 calls.out.1.txt > tmphead")
  system("grep -P -f tmphead -v calls.out.txt > tmp.calls.out.txt")
  system("cat tmphead tmp.calls.out.txt > calls.out.txt")
  system("chmod 777 calls.out.txt")
  system("rm tmp.calls.out.txt")
  
  os.call.2 <- paste("cat ",
                     paste("segmented.out.", seqc, ".txt", sep = "", collapse = " "),
                     " > segmented.out.txt")
  system(os.call.2)
  system("head -1 segmented.out.1.txt > tmphead")
  system("grep -P -f tmphead -v segmented.out.txt > tmp.segmented.out.txt")
  system("cat tmphead tmp.segmented.out.txt > segmented.out.txt")
  system("chmod 777 segmented.out.txt")
  system("rm tmp.segmented.out.txt")
  system("rm tmphead")
}
  



## new.custom <- function(segmresRDataName = "segmres.RData",
##                        cghRDataName = "cghData.RData",
##                        chromRDataName = "chromData.RData",
##                        posRDataName = "posData.RData",
##                        probeNamesRDataName = "probeNames.RData",
##                        nround = 6) {
##   load(posRDataName)
##   load(chromRDataName)
##   load(cghRDataName)
##   load(segmresRDataName)
##   load(probeNamesRDataName)
  
##   narrays <- ncol(segmres[[1]])
##   seqi <- seq.int(1, narrays)

##   if(!is.factor(probeNames))
##     probeNames <- factor(probeNames)
  
##   list.1 <- list(ProbeName = as.ff(probeNames, vmode = NULL),
##                  Chr = chromData,
##                  Position = posData)
##   ## close(list.1[[1]]) ## Don't close; must be open for
##   ## later writing
##   rm(probeNames)
##   gc()
##   l.smoothed <- lapply(seqi,
##                        function(i) segmres[[1]][[i]])
##   l.calls <- lapply(seqi,
##                     function(i) segmres[[2]][[i]])
##   l.original <- lapply(seqi,
##                        function(i) cghData[[i]])

##   names(l.smoothed) <- names(l.calls) <- names(l.original) <-
##     names(segmres[[1]])

##   l.smoothed <- c(list.1, l.smoothed)
##   l.calls <- c(list.1, l.calls)
##   l.original <- c(list.1, l.original)

##   segmentedffdf <- do.call("ffdf", l.smoothed)
##   callsffdf <- do.call("ffdf", l.calls)
##   originalffdf <- do.call("ffdf", l.original)

##   callsffdf.no.names <- do.call("ffdf", l.calls[-1])
  
##   open(segmentedffdf)
##   open(callsffdf)
##   open(originalffdf)
  
##   ## write.table.ffdf(segmentedffdf, file = "segmented.out.txt",
##   ##                  sep = "\t", quote = FALSE)
##   ## write.table.ffdf(callsffdf, file = "calls.out.txt",
##   ##                  sep = "\t", quote = FALSE)


##   ############  Change this to something faster/better/before
##   ############  How does it work with many columns?
##   ############  Could use same logic for chromosome and position
##   write("ProbeName", file = "pn2.txt")
##   write(probeNames, file = "pn2.txt",
##         append = TRUE)
##   write.table.ffdf(callsffdf.no.names, file = "call.no.name.txt",
##                    sep = "\t", quote = FALSE)
##   system("paste pn2.txt call.no.name,txt > calls.out.txt")
##   #################################

  
##   rle.chr <- intrle(as.integer(chromData[]))
##   chr.end <- cumsum(rle.chr$lengths)
##   chr.start <- c(1, chr.end[-length(chr.end)] + 1)
##   seqc <- seq.int(1, length(chr.start))

##   f1 <- function(i,  objectin, nameout) {
##     oname <- paste(nameout, i, sep = "")
##     assign(oname,
##            objectin[ri(chr.start[i], chr.end[i]), ])
##     save(file = paste(oname, ".RData", sep = ""),
##          list = c(oname), compress = FALSE)
##   }

##   ## the following three lines take their time!!
##   ## parallelize!! yes, via clusterApplyLB??
##   null <- sapply(seqc, f1, segmentedffdf, "segmented.out.")
##   null <- sapply(seqc, f1, callsffdf, "calls.out.")
##   null <- sapply(seqc, f1, originalffdf, "original.")
  
##   null <- close(segmentedffdf)
##   null <- close(callsffdf)
##   null <- close(originalffdf)
##   return(0)
## }
  




##############################################
##############################################
######                              ##########
######         Start execution      ##########
######                              ##########
##############################################
##############################################

if (.__ADaCGH_SERVER_APPL) {
  version
  system("date")
  system("hostname")
  cat("\nRunning\n", file = "R_Status.txt")
  hostn <- system("hostname", intern = TRUE)
  pid <- Sys.getpid()
  cat("\nPID is ", pid, "\n")
  
  startExecTime <- format(Sys.time())
  write.table(file = "pid.txt", pid,
              row.names = FALSE,
              col.names = FALSE)
  
  sink(file = "current_R_proc_info")
  cat(hostn)
  cat("  ")
  cat(pid)
  cat("\n")
  sink()

## attach pid to name in R.running.procs
  new.name1 <- unlist(strsplit(getwd(), "/"))
  new.name1 <- paste(new.name1[length(new.name1)], "@", hostn, sep = "")
  new.name <- paste("R.", new.name1, "%", pid, sep = "")
  new.name1 <- paste("R.", new.name1, sep = "")
  system(paste("mv ../../R.running.procs/", new.name1,
               " ../../R.running.procs/", new.name,
               sep = ""))
  
  checkpoint.num <- scan("checkpoint.num", what = double(0), n = 1)
} else {
  checkpoint.num <- 0
}


#######################################################
#######################################################
#######################################################
###
###       Getting and checking options
###
#######################################################
#######################################################
#######################################################


  
trythis <- try({WaviOptions <- readOptions("options.txt")})
if(inherits(trythis, "try-error"))
  caughtUserError.Web(trythis)

## checkAssign("idtype", acceptedIDTypes, WaviOptions)
## checkAssign("organism", acceptedOrganisms, WaviOptions)
checkAssign("method", acceptedMethodaCGH, WaviOptions)
checkAssign("colorNoChange", acceptedColors, WaviOptions)
checkAssign("colorGain", acceptedColors, WaviOptions)
checkAssign("colorLoss", acceptedColors, WaviOptions)
checkAssign("colorSmooth", acceptedColors, WaviOptions)


WaviOptions$colorsWavi <- c(WaviOptions$colorNoChange, WaviOptions$colorGain,
                        WaviOptions$colorLoss, WaviOptions$colorSmooth,
                        "black")

## this is not a good idea, but I leave so as not to have to change
## the options files. But we do more processing of options below
WaviOptions <- checkConvertMethodOptions(methodOptions, WaviOptions)    

## merging options as mergeRes, etc, are ignored.

## we do not want to mess around with options anymore
try2 <- try({
  
  if(WaviOptions$method == "HaarSeg") {
    if(!is.null(WaviOptions$HaarSeg.m)) {
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS mad.threshold in HaarSeg\n\n\n")
      WaviOptions$mad.threshold <- WaviOptions$HaarSeg.m
    }
  }
  
  if(!is.null(WaviOptions$mad.threshold)) {
    WaviOptions$mad.threshold <- as.numeric(WaviOptions$mad.threshold)
  } else {
    cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS mad.threshold\n\n\n")
    WaviOptions$mad.threshold <- 2
  }
  
  
  
  
  if(is.null(WaviOptions[["merge"]])) {
    if(WaviOptions$method == "DNAcopy")
      WaviOptions$merge <- "mergeLevels"
    if(WaviOptions$method == "CGHseg") 
      WaviOptions$merge <- "MAD"
    if(WaviOptions$method == "Wavelets")
      WaviOptions$merge <- "MAD"
  }

  if(!is.null(WaviOptions$DNAcopy.min.width)) {
    WaviOptions$DNAcopy.min.width <- as.numeric(WaviOptions$DNAcopy.min.width)
    if((WaviOptions$DNAcopy.min.width > 5) |
       (WaviOptions$DNAcopy.min.width < 2))
      caughtUserError.Web("min.width for DNAcopy must be between 2 and 5")
  } else {
    if(WaviOptions$method == "DNAcopy")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS DNAcopy.min.width \n\n\n")
    WaviOptions$DNAcopy.min.width <- 2
  }
  if(!is.null(WaviOptions$DNAcopy.alpha)) {
    WaviOptions$DNAcopy.alpha <- as.numeric(WaviOptions$DNAcopy.alpha)
    if((WaviOptions$DNAcopy.alpha >= 1) |
       (WaviOptions$DNAcopy.alpha <= 0))
      caughtUserError.Web("alpha for DNAcopy must be between 0 and 1")
  } else {
    if(WaviOptions$method == "DNAcopy")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS DNAcopy.alpha\n\n\n")
    WaviOptions$DNAcopy.alpha <- 0.01
  }
  if(!is.null(WaviOptions$DNAcopy.nperm)) {
    WaviOptions$DNAcopy.nperm <- as.numeric(WaviOptions$DNAcopy.nperm)
    if(WaviOptions$DNAcopy.nperm < 1)
      caughtUserError.Web("nperm for DNAcopy must > 1")
  } else {
    if(WaviOptions$method == "DNAcopy")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS DNAcopy.nperm\n\n\n")
    WaviOptions$DNAcopy.nperm <- 10000
  }

  if(is.null(WaviOptions$GLAD.deltaN)) {
    if(WaviOptions$method == "GLAD")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS GLAD.deltaN\n\n\n")
    WaviOptions$GLAD.deltaN <- 0.10
  } else {
    WaviOptions$GLAD.deltaN <- as.numeric(WaviOptions$GLAD.deltaN)
  }
  if(is.null(WaviOptions$GLAD.forceGL1)) {
    if(WaviOptions$method == "GLAD")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS GLAD.forceGL1\n\n\n")
    WaviOptions$GLAD.forceGL1 <- -0.15
  } else {
    WaviOptions$GLAD.forceGL1 <- as.numeric(WaviOptions$GLAD.forceGL1)
  }
  if(is.null(WaviOptions$GLAD.forceGL2)) {
    if(WaviOptions$method == "GLAD")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS GLAD.forceGL2\n\n\n")
    WaviOptions$GLAD.forceGL2 <- 0.15
  } else {
    WaviOptions$GLAD.forceGL2 <- as.numeric(WaviOptions$GLAD.forceGL2)
  }
  WaviOptions$GLAD.forceGL <- c(WaviOptions$GLAD.forceGL1,
                                WaviOptions$GLAD.forceGL2)
  
  if(is.null(WaviOptions$GLAD.deletion)) {
    if(WaviOptions$method == "GLAD")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS GLAD.deletion\n\n\n")
    WaviOptions$GLAD.deletion <- -5
  } else {
    WaviOptions$GLAD.deletion <- as.numeric(WaviOptions$GLAD.deletion)
  }
  if(is.null(WaviOptions$GLAD.amplicon)) {
    if(WaviOptions$method == "GLAD")
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS GLAD.amplicon\n\n\n")
    WaviOptions$GLAD.amplicon <- 1
  } else {
    WaviOptions$GLAD.amplicon <- as.numeric(WaviOptions$GLAD.amplicon)
  }
  if(is.null(WaviOptions$HMM.AIC.BIC)) {
    if(WaviOptions$method %in% c("HMM", "BioHMM"))
      cat("\n\n\n WARNING!!!! OLD FORMAT FOR OPTIONS HMM.AIC.BIC\n\n\n")
    WaviOptions$HMM.AIC.BIC <- "AIC"
  }
  
})


if(inherits(try2, "try-error"))
  caughtUserError.Web("Options file incorrect")


## if(is.null(WaviOptions$mergeRes)) {
##   if(WaviOptions$method == "DNAcopy")
##     WaviOptions$merge <- "mergeLevels"
##   if(WaviOptions$method == "CGHseg")
##     WaviOptions$merge <- "MAD"
##   if(WaviOptions$method == "Wavelets")
##     WaviOptions$merge <- "MAD"
## }



cat("\n\n All WaviOptions are: \n")

WaviOptions


#######################################################
#######################################################
#######################################################
###
###       Read data and initial stuff
###
#######################################################
#######################################################
#######################################################

if(checkpoint.num < 1) {
  
  ## With ff and new functions. Reading data, etc, is done in a child process
  ## that does only that. That way, main process does not use a lot of RAM
  ## Verify we are doing OK killing the child, or whatever. I think we are.

  ## Recall multicore leaves zombies. But we cannot use "fork" in
  ## package fork, since fork does not permit
  ## getting back the utput of a function.
  library(multicore)
  parallel(inputDataToADaCGHData(), silent = FALSE)
  tableChromArray <- collect()[[1]]
  if(inherits(tableChromArray, "try-error")) {
    caughtOurError.Web("ERROR in input data conversion")
  }

  numarrays <- max(tableChromArray$ArrayNum)
  chromnum <- max(tableChromArray$ChromNum)

  cat("\n gc right before checkpoint 1 \n")
  print(gc())

  to.save <- c("numarrays", "chromnum", "new.custom2",
               "put.part.rdata.together",
               "tableChromArray",
               "WaviOptions", ".__ADaCGH_SERVER_APPL",
               "doCheckpoint", "NormalTermination",
               "caughtUserError.Web", "caughtOurError.Web")
  checkpoint.num <- doCheckpoint(1, to.save)
  
  cat("\n gc right after checkpoint 1 \n")
  gc()
}

################################################################
## MPI, LAM, etc.
## enter info into lam suffix log table


if (.__ADaCGH_SERVER_APPL) {
  trylam <- try(
                lamSESSION <- scan("lamSuffix", sep = "\t",
                                   what = "",
                                   strip.white = TRUE)
                )
  tmpDir <- getwd()
  hostn <- system("hostname", intern = TRUE)
  pid <- Sys.getpid()

  sed.command <- paste("sed -i 's/RprocessPid\\t",
                       lamSESSION, "\\t", hostn, "/",
                       pid, "\\t",
                       lamSESSION, "\\t", hostn, "/' ",
                       "/http/adacgh-server/runs-tmp/logs/LAM_SUFFIX_Log",
                       ##                     "/http/mpi.log/LAM_SUFFIX_Log",
                       sep = "")
  
  system(sed.command)
  cat("\n\n Did sed.command ")
  cat(sed.command)
}


#######################################################
#######################################################
#######################################################
###
###       MPI and snowfall: launch
###
#######################################################
#######################################################
#######################################################

options(warn = -1)

### Launch Rmpi as late as possible with only the minimum possible slaves

library(Rmpi)

print(system("lamnodes"))
print(paste("Universe size is ", mpi.universe.size()))

usize <- min(numarrays * chromnum, mpi.universe.size())
## make sure at least two, o.w. rsprng won't work, and
## we do not want to hack my mpiInit.
print(paste("usize is", usize))

if (.__ADaCGH_SERVER_APPL) {
  usize <- max(2, usize)
  snowfallInit(universeSize = usize, exit_on_fail = TRUE)
  print("after mpiInit")
  cat("\n\nAbout to print mpiOK file\n")
  sink(file = "mpiOK")
  cat("MPI started OK\n")
  sink()
} else {
  snowfallInit(universeSize = mpi.universe.size(), exit_on_fail = FALSE)
}


#######################################################
#######################################################
#######################################################
###
###       Segmentation
###
#######################################################
#######################################################
#######################################################


if(checkpoint.num < 3) {

  trythis <- try({
    fseg <- get(paste("pSegment", WaviOptions$method, sep = ""))
    segmres <- fseg(cghRDataName = "cghData.RData",
                    chromRDataName = "chromData.RData",
                    posRDataName = "posData.RData",
                    merging = WaviOptions$merge,
                    minDiff = WaviOptions$Wave.minDiff,
                    CGHseg.thres = WaviOptions$CGHseg.s,
                    mad.threshold = WaviOptions$mad.threshold,
                    min.width = WaviOptions$DNAcopy.min.width,
                    alpha = WaviOptions$DNAcopy.alpha,
                    nperm = WaviOptions$DNAcopy.nperm,
                    deltaN = WaviOptions$GLAD.deltaN,
                    forceGL = WaviOptions$GLAD.forceGL,
                    deletion = WaviOptions$GLAD.deletion,
                    amplicon = WaviOptions$GLAD.amplicon,
                    aic.or.bic = WaviOptions$HMM.AIC.BIC)
  })
  
  if(inherits(trythis, "try-error"))
    caughtOurError.Web(trythis)
  cat("\n\n Segmentation done \n\n")
  save(segmres, file = "segmres.RData", compress = FALSE)
 
  cat("\n gc right before checkpoint 3 \n")
  print(gc())

  to.save <- c("numarrays", "chromnum",
               "tableChromArray", "new.custom2",
               "put.part.rdata.together",
               "WaviOptions", ".__ADaCGH_SERVER_APPL",
               "doCheckpoint", "NormalTermination",
               "caughtUserError.Web", "caughtOurError.Web")
  
  checkpoint.num <- doCheckpoint(3, to.save)
   cat("\n gc right after checkpoint 3 \n")
  print(gc())



  
}





if(checkpoint.num < 5) {

  trythis <- try(
                 pChromPlot(outRDataName = "segmres.RData",
                            cghRDataName = "cghData.RData",
                            chromRDataName = "chromData.RData",
                            posRDataName = "posData.RData",
                            probenamesRDataName = "probeNames.RData",
                            colors = WaviOptions$colorsWavi,
                            imgheight = 350)
                 )
  if(inherits(trythis, "try-error"))
    caughtOurError.Web(trythis)
  cat("\n\n Color plotting done \n\n")
  cat("\n gc right after plotting \n")
  print(gc())

  
  ## Plots without colours
  ## As plots would get overwritten I create a new dir, etc.
  dir1 <- getwd()
  dir.create("BW")
  setwd("BW")
  print(getwd())
  sfClusterEval(setwd("BW"))
##  mpi.remote.exec(setwd("BW"))
  trythis <- try(
                 pChromPlot(outRDataName = "../segmres.RData",
                            cghRDataName = "../cghData.RData",
                            chromRDataName = "../chromData.RData",
                            posRDataName = "../posData.RData",
                            probenamesRDataName = "../probeNames.RData",
                            colors = c(rep("black", 3), "blue"),
                            imgheight = 350)
                 )
  if(inherits(trythis, "try-error"))
    caughtOurError.Web(trythis)
  files.in.BW <- dir()
  for (ffbw in files.in.BW) file.rename(ffbw, paste("BW_", ffbw, sep = ""))
  null <- file.copy(from = dir(), to = dir1)
  ##         system('mmv "*" "BW_#1"')
  ##         system('cp * ../.')
  setwd(dir1)
  sfExport("dir1")
  sfClusterEval(setwd(dir1))
##  mpi.remote.exec(setwd(dir1))
  
  cat("\n\n BW plotting done \n\n")
  cat("\n gc right after plotting \n")
  print(gc())

 to.save <- c("numarrays", "chromnum",
               "tableChromArray", "new.custom2",
               "put.part.rdata.together",
               "WaviOptions", ".__ADaCGH_SERVER_APPL",
               "doCheckpoint", "NormalTermination",
               "caughtUserError.Web", "caughtOurError.Web")
  
  checkpoint.num <- doCheckpoint(5, to.save)
   cat("\n gc right after checkpoint 5 \n")
  print(gc())
}


if(checkpoint.num < 7) {
  cat("\n Before writing output files\n")
  print(date())
  new.custom2()
  cat("\n After writing output files\n")
  print(date())
  print(gc())
  system("tar -czf results.txt.tar.gz calls.out*txt segmented.out*txt")
  system("tar -czf results.RData.tar.gz calls.out.*.RData segmented.out.*.RData original.*.RData")
  NormalTermination()
}

