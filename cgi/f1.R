####  Copyright (C) 2005, 2006, Ramon Diaz-Uriarte <rdiaz02@gmail.com>

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



rm(list = ls()) ## Just in case.

.Last <- function(){
    ## the next four lines are a way to ensure that the OS writes
    ## a file that could be immediately read by Python to see
    ## if we are done
    RterminatedOK <- file("RterminatedOK", "w")
    cat("\nNormal termination\n", file = RterminatedOK)
    flush(RterminatedOK)
    close(RterminatedOK)
    ##save.image()
    if (is.loaded("mpi_initialize")){ 
        if (mpi.comm.size(1) > 0){ 
        try(print("Please use mpi.close.Rslaves() to close slaves."), silent = TRUE)
        try(mpi.close.Rslaves() , silent = TRUE)
        } 
        try(print("Please use mpi.quit() to quit R"), silent = TRUE)
        cat("\n\n Normal termination\n")
        try(mpi.quit(save = "yes"), silent = TRUE)
    }
    cat("\n\n Normal termination\n")
    ## In case the CGI is not called (user kills browser)
    ## have a way to stop lam
    try(system(paste("/http/mpi.log/killLAM.py", lamSESSION, "&")))
    try(mpi.quit(save = "yes"), silent = TRUE)
}


##startExecTime <- format(Sys.time())

pid <- Sys.getpid()
write.table(file = "pid.txt", pid,
            row.names = FALSE,
            col.names = FALSE)


## attach pid to name in R.running.procs
hostn <- system("hostname", intern = TRUE)
new.name1 <- unlist(strsplit(getwd(), "\/"))
new.name1 <- paste(new.name1[length(new.name1)], "@", hostn, sep = "")
new.name <- paste("R.", new.name1, "%", pid, sep = "")
new.name1 <- paste("R.", new.name1, sep = "")
system(paste("mv ../../R.running.procs/", new.name1,
             " ../../R.running.procs/", new.name,
             sep = ""))




########################################################

########   Start MPI here to check if everything OK

#########################################################

library(Rmpi)
mpi.spawn.Rslaves(nslaves= mpi.universe.size())
library(papply)
sink(file = "mpiOK")
cat("MPI started OK\n")
sink()

assign(".__ADaCGH_WEB_APPL", TRUE)
print("testing existence of indicator")
print(exists(".__ADaCGH_WEB_APPL"))
library(Hmisc)
library(cgh)
library(aCGH)
library("waveslim") ## we will have to load ADaCGH soon,
## but we must mask certain defs. in waveslim. So load
## waveslim here
library("cluster") 
library(cghMCR)
library(DNAcopy)
library(GDD)
library(imagemap)
library(ADaCGH)



startExecTime <- format(Sys.time())

pid <- Sys.getpid()
write.table(file = "pid.txt", pid,
            row.names = FALSE,
            col.names = FALSE)
trylam <- try(
              lamSESSION <- scan("lamSuffix", sep = "\t",
                                 strip.white = TRUE)
              )


png.width = 400
png.height = 400
png.pointsize = 10
# png.family = "Helvetica"

## defaults for DNA copy
DNA.undo.splits = "prune" ## don't touch this
DNA.undo.sd = 3  ## not needed, really


##############################################
##############################################
######                              ##########
######         Error checking       ##########
######          utilities           ##########
######                              ##########
##############################################
##############################################


caughtUserError <- function(message) {
    GDD("ErrorFigure.png", width = png.width,
           height = png.height, 
           ps = png.pointsize)
##           family = png.family)
    plot(x = c(0, 1), y = c(0, 1),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.7, "There was a PROBLEM with your data.")
    text(0.5, 0.5,
    "Please read carefully the error messages,")
    
    text(0.5, 0.3, "fix the problem, and try again.")
    dev.off()
    sink(file = "results.txt")
    cat(message)
    sink()
    sink(file = "exitStatus")
    cat("Error\n\n")
    cat(message)
    sink()
    quit(save = "no", status = 11, runLast = TRUE)
}

## FIXME: this is here and in the package code. Can we eliminate
##        from here? But it gives a clear message structure: User and
##        ours.
caughtOurError <- function(message) {
    GDD("ErrorFigure.png", width = png.width,
           height = png.height, 
           ps = png.pointsize)
##           family = png.family)
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
}

warningsForUsers <- vector()

#######################################################
#######################################################
#######################################################
###
###       Read data and initial stuff
###
#######################################################
#######################################################
#######################################################



## The xdata can have either a single copy of each clone, or
## multiple. It does not matter. Averaging is done w.r.t Name in
## "positionInfo". Thus, if a single clone, all clones of a
## multi-clone gene have same weight. This can ve achieved by first
## using the preprocessor. Otherwise, each clone can have a different
## weight if diff. clones have different number of copies in the
## array.


## zz: future: allow merging by clone.then map clones, then merge by
## name and pos.


idtype <- try(scan("idtype", what = "", n = 1))
organism <- try(scan("organism", what = "", n = 1))

MCR.gapAllowed <- try(scan("MCR.gapAllowed", what = double(0), n = 1))
MCR.alteredLow <- try(scan("MCR.alteredLow", what = double(0), n = 1))
MCR.alteredHigh <- try(scan("MCR.alteredHigh", what = double(0), n = 1))
MCR.recurrence <- try(scan("MCR.recurrence", what = double(0), n = 1))



methodaCGH <- scan("methodaCGH", what = "", n = 1)
if (methodaCGH == "CBS") {
    DNA.smooth.region    <- scan("DNA.smooth.region", what = double(0), n = 1)   
    DNA.outlier.SD.scale <- scan("DNA.outlier.SD.scale", what = double(0), n = 1)
    DNA.smooth.SD.scale  <- scan("DNA.smooth.SD.scale", what = double(0), n = 1) 
    DNA.trim             <- scan("DNA.trim", what = double(0), n = 1)
    DNA.nmin             <- scan("DNA.nmin", what = double(0), n = 1)
    DNA.kmax             <- scan("DNA.kmax", what = double(0), n = 1)   
    DNA.copy.alpha       <- scan("DNA.copy.alpha", what = double(0), n = 1)      
    DNA.nperm            <- scan("DNA.nperm", what = double(0), n = 1)           
    DNA.overlap          <- scan("DNA.overlap", what = double(0), n = 1)         
    DNA.undo.prune       <- scan("DNA.undo.prune", what = double(0), n = 1)
    DNA.merge            <- scan("DNA.merge", what = "", n = 1)          
} else if (methodaCGH == "WS") {
    Wave.minDiff <-  scan("Wave.minDiff", what = double(0), n = 1)
} else if (methodaCGH == "PSW") {
    PSW.nIter <- scan("PSW.nIter", what = double(0), n = 1)
    PSW.prec <- scan("PSW.prec", what = double(0), n = 1)
    PSW.p.crit <- scan("PSW.p.crit", what = double(0), n = 1)
} else if (methodaCGH == "ACE") {
    ACE.fdr <- scan("ACE.fdr", what = double(0), n = 1)
} else { ## nothing else for now
    caughtUserError("This method is not yet implemented.")
}


twoFiles <- try(scan("twofiles", what = "", n = 1))
centering <- try(scan("centering", what = "", n = 1))

## positionInfo         has  name, chromosome, start, end
trypositionInfo <-
    try(
        positionInfo <- read.table("positionInfo", header = FALSE,
                                   sep = "\t",
                                   strip.white = TRUE,
                                   comment.char = "#",
                                   quote = ""))
if(class(trypositionInfo) == "try-error")
    caughtUserError(paste("The position file is not of the appropriate format\n",
                    "In case it helps this is the error we get\n",
                    trypositionInfo, sep =""))

if(ncol(positionInfo) > 4) positionInfo <- positionInfo[, 1:4]

arrayNames <- scan("arrayNames", sep = "\t", what = "char", quote = "")
if(length(arrayNames) > 0) {
    ##arrayNames <- arrayNames[-1]
    if(length(unique(arrayNames)) < length(arrayNames)) {
        dupnames <- which(duplicated(arrayNames))
        message <- paste("Array names are not unique.\n",
                         "Please change them so that they are unique.\n",
                         "The duplicated names are ", dupnames, "\n")
        caughtUserError(message)
    }
}

tryxdata <-
    try(
        xdata <- scan("covarR", what = double(0), sep = "\t")
        )
if(class(tryxdata) == "try-error")
    caughtUserError(paste("The acgh data file is not of the appropriate format\n",
                          "In case it helps this is the error we get\n",
                          tryxdata, sep =""))

xdata <- matrix(xdata, ncol = length(arrayNames), byrow = TRUE)

if(ncol(xdata) < 1)
    caughtUserError(paste("The acgh data file does not contain any data\n",
                          "(recall that the first column is only identifiers)\n"))

if(nrow(xdata) != nrow(positionInfo))
    caughtUserError(paste("Different number of genes/clones in your\n",
                          "data and position (coordinate) files.\n",
                          nrow(xdata), "genes/clones in you data file\n",
                          nrow(positionInfo), "genes/clones in you positions file.\n"))
colnames(xdata) <- arrayNames

colnames(positionInfo) <- c("name", "chromosome", "start", "end")

if(!length(arrayNames))
    cat(paste(colnames(xdata), collapse = "\t"), file = "arrayNames")

if(!(is.numeric(as.matrix(xdata)))) {
    caughtUserError("Your aCGH file contains non-numeric data. \n That is not allowed.\n")
}
if(any(is.na(xdata))) {
    caughtUserError("Your aCGH file contains missing values. \n That is not allowed.\n")
}
if(any(is.na(positionInfo))) {
    caughtUserError("The position (coordinate) information contains missing values.\n That is not allowed.\n")
}
if(any(!is.numeric(as.matrix(positionInfo[, c(3, 4)])))) {
    caughtUserError(paste("Your position information contains non-numeric values \n",
                          "for the start and/or end positions."))
}

## Get rid of possible chr, Chr, etc. and possible " in the chr
positionInfo$chromosome <- as.character(positionInfo$chromosome)
positionInfo$chromosome <- sub("\"", "", positionInfo$chromosome)
positionInfo$chromosome <- sub("\"", "", positionInfo$chromosome)
positionInfo$chromosome <- sub("chr", "", positionInfo$chromosome)
positionInfo$chromosome <- sub("chr ", "", positionInfo$chromosome)
positionInfo$chromosome <- sub("Chr", "", positionInfo$chromosome)
positionInfo$chromosome <- sub("Chr ", "", positionInfo$chromosome)

weird.chromos <- which(!positionInfo$chromosome %in%
                       c("X", "Y", 1:100))

if(length(weird.chromos)) {
    warningsForUsers <-
        c(warningsForUsers,
          paste("There were", length(weird.chromos),
                "clones/genes with a chromosome which was neither",
                "1 to 100 or X or Y; these have been excluded ",
                "from further analyses."))
    if (length(weird.chromos) > (dim(positionInfo)[1]/2)) {
        caughtUserError("More than half of your data have chromosomes with values that are neither an integer 1:100 or X, Y")
    }
    positionInfo <- positionInfo[-weird.chromos, ]
    xdata <- xdata[-weird.chromos, , drop = FALSE]
}

positionInfo$chromosome <- factor(positionInfo$chromosome)
positionInfo$MidPoint <- positionInfo$start +
    0.5 * (positionInfo$end - positionInfo$start)

chrom.numeric <- as.numeric(as.character(positionInfo$chromosome))
chrom.numeric[positionInfo$chromosome == "X"] <- 23
chrom.numeric[positionInfo$chromosome == "Y"] <- 24
positionInfo$chrom.numeric <- chrom.numeric
ncrom <- length(unique(chrom.numeric))
rm(chrom.numeric)
reorder <- order(positionInfo$chrom.numeric,
                 positionInfo$MidPoint,
                 positionInfo$start,
                 positionInfo$end,
                 positionInfo$name)

positionInfo <- positionInfo[reorder, ]
xdata <- xdata[reorder, , drop = FALSE]

#### Make chromosome a numeric variable
### make more sophisticated later: return X and Y in plots. zz



#######################################################
#######################################################
#######################################################
###
###       Averaging by clones
###
#######################################################
#######################################################
#######################################################

### ordering variable

ov <- paste(positionInfo$chrom.numeric,
            positionInfo$MidPoint,
            positionInfo$name,
            sep = "*")
positionInfo$ov <- factor(ov, levels = unique(ov),
                          labels = unique(ov))
                          
### By name
xdata.merge1 <- apply(xdata, 2,
                      function(x) {
                          unlist(tapply(x,
                                        positionInfo$ov,
                                        function(z) {mean(z)}))})
positions.merge1 <- positionInfo[!duplicated(positionInfo$ov), ]

## Do we have any identical MidPos in the same chromosome??  Just to solve
## it quickly and without nasty downstream consequences, we add a runif to
## midPos.

tmp <- paste(positions.merge1$chromosome, positions.merge1$MidPoint, sep = ".")
tmp <- factor(tmp, levels = unique(tmp), labels = unique(tmp))
if (sum(duplicated(tmp))) {
    ## add a random variate, to break ties:
    positions.merge1$MidPoint[duplicated(tmp)] <-
        positions.merge1$MidPoint[duplicated(tmp)] +
            runif(sum(duplicated(tmp)))

    ## Reorder, just in case
    reorder <- order(positions.merge1$chrom.numeric,
                     positions.merge1$MidPoint,
                     positions.merge1$start,
                     positions.merge1$end,
                     positions.merge1$name)
    
    positions.merge1 <- positions.merge1[reorder, ]
    xdata <- xdata[reorder, , drop = FALSE]
}

tmp <- paste(positions.merge1$chromosome, positions.merge1$MidPoint, sep = ".")
if (sum(duplicated(tmp)))
    stopOurError("still duplicated MidPoints; shouldn't happen")


## below, we fix the distance between first obs. of each chromos
positions.merge1$DistanceClones <- c(NA, diff(positions.merge1$MidPoint))

## get ranks
obs.per.chrom <- as.vector(table(positions.merge1$chrom.numeric))
positions.merge1$rank <- unlist(sapply(obs.per.chrom,
                                    function(x) {seq(from = 1, to = x)}))

## fixing the distance of the first clone in each chromosome
positions.merge1$rank[positions.merge1$rank == 1] <- NA


#######################################################
#######################################################
#######################################################
###
###       Rescaling data
###
#######################################################
#######################################################
#######################################################

### Provide simple statistics: by array, and by array by chromosome.


means <- apply(xdata.merge1, 2, mean)
medians <- apply(xdata.merge1, 2, median)
mads <- apply(xdata.merge1, 2, mad)


sink(file = "results.for.html")
cat("<h3>Settings</h3>\n")

cat("<h4>Method:              ", methodaCGH,"</h4>\n")

cat("<h4>Method parameters </h4>\n")
if (methodaCGH == "CBS") {
    cat("<p>        DNA.smooth.region\t\t:     ",    DNA.smooth.region,"</p>\n")   
    cat("<p>        DNA.outlier.SD.scale\t\t:  ", DNA.outlier.SD.scale,"</p>\n")
    cat("<p>        DNA.smooth.SD.scale\t\t:   ",  DNA.smooth.SD.scale,"</p>\n") 
    cat("<p>        DNA.trim\t\t:              ",             DNA.trim,"</p>\n")            
    cat("<p>        DNA.nmin\t\t:              ",             DNA.nmin,"</p>\n")            
    cat("<p>        DNA.kmax\t\t:              ",             DNA.kmax,"</p>\n")            
    cat("<p>        DNA.copy.alpha\t\t:        ",       DNA.copy.alpha,"</p>\n")      
    cat("<p>        DNA.nperm\t\t:             ",            DNA.nperm,"</p>\n")           
    cat("<p>        DNA.overlap\t\t:           ",          DNA.overlap,"</p>\n")         
    cat("<p>        DNA.undo.prune\t\t:        ",       DNA.undo.prune,"</p>\n")
    cat("<p>        DNA.merge\t\t:             ",       DNA.merge,"</p>\n")
} else if (methodaCGH == "WS") {
    cat("<p>        Wave.minDiff\t\t:          ",         Wave.minDiff,"</p>\n")
} else if (methodaCGH == "PSW") {
    cat("<p>        PSW.nIter\t\t:             ",             PSW.nIter,"</p>\n")
    cat("<p>        PSW.prec\t\t:              ",             PSW.prec,"</p>\n")
    cat("<p>        PSW.p.crit\t\t:            ",             PSW.p.crit,"</p>\n")
}


cat("<h4>Minimal common regions </h4>\n")
cat("<p>        MCR.gapAllowed\t\t:             ",       MCR.gapAllowed,"</p>\n")
cat("<p>        MCR.alteredLow\t\t:             ",       MCR.alteredLow,"</p>\n")            
cat("<p>        MCR.alteredHigh\t\t:             ",       MCR.alteredHigh,"</p>\n")
cat("<p>        MCR.recurrence\t\t:             ",       MCR.recurrence,"</p>\n")            




cat("<h4>Centering:              </h4>", centering,"\n")
sink()






sink(file = "results.txt")

cat("\n\n\n*********************************************************************\n")
cat("*********************************************************************\n")
cat("*********                                        ********************\n")
cat("*********              Settings                  ********************\n")
cat("*********                                        ********************\n")
cat("*********************************************************************\n")
cat("*********************************************************************\n\n")

cat("\n\n Method:              ", methodaCGH)

cat("\n\n Parameters:")
if (methodaCGH == "CBS") {
    cat("\n\nDNA.smooth.region\t\t:",    DNA.smooth.region)   
    cat("\n\nDNA.outlier.SD.scale\t\t:  ", DNA.outlier.SD.scale)
    cat("\n\nDNA.smooth.SD.scale\t\t:   ",  DNA.smooth.SD.scale) 
    cat("\n\nDNA.trim\t\t:              ",             DNA.trim)            
    cat("\n\nDNA.nmin\t\t:              ",             DNA.nmin)            
    cat("\n\nDNA.kmax\t\t:              ",             DNA.kmax)            
    cat("\n\nDNA.copy.alpha\t\t:        ",       DNA.copy.alpha)      
    cat("\n\nDNA.nperm\t\t:             ",            DNA.nperm)           
    cat("\n\nDNA.overlap\t\t:           ",          DNA.overlap)         
    cat("\n\nDNA.undo.prune\t\t:        ",       DNA.undo.prune)
    cat("\n\nDNA.merge\t\t:             ",       DNA.merge)        
} else if (methodaCGH == "WS") {
    cat("\n\nWave.minDiff\t\t:          ",         Wave.minDiff)
} else if (methodaCGH == "PSW") {
    cat("\n\nPSW.nIter\t\t:             ",             PSW.nIter)
    cat("\n\nPSW.prec\t\t:              ",             PSW.prec)
    cat("\n\nPSW.p.crit\t\t:            ",             PSW.p.crit)
}

cat("\n\n\n\nMinimal common regions\n")
cat("\n\nMCR.gapAllowed\t\t:             ",       MCR.gapAllowed)
cat("\n\nMCR.alteredLow\t\t:             ",       MCR.alteredLow)            
cat("\n\nMCR.alteredHigh\t\t:             ",       MCR.alteredHigh)
cat("\n\nMCR.recurrence\t\t:             ",       MCR.recurrence)            

cat("\n\n\n\n Centering:              ", centering)




cat("\n\n\n*********************************************************************\n")
cat("*********************************************************************\n")
cat("*********                                        ********************\n")
cat("*********        Clones/genes per chromosome     ********************\n")
cat("*********                                        ********************\n")
cat("*********************************************************************\n")
cat("*********************************************************************\n\n")

tmp <- table(positions.merge1$chrom.numeric)
tmp2 <- data.frame(as.vector(tmp))
rownames(tmp2) <- names(tmp)
colnames(tmp2) <- "Number of genes/clones per chromosome"
tmp2
tmphtml <- html(tmp2, file = "clones.per.chrom.html", dec = 4, first.col = "Chromosome")
                
rm(tmp, tmp2, tmphtml)
    


cat("\n\n\n*********************************************************************\n")
cat("*********************************************************************\n")
cat("*********                                        ********************\n")
cat("*********       Basic statistics by array        ********************\n")
cat("*********        and array by chromosome         ********************\n")
cat("*********           BEFORE  CENTERING            ********************\n")
cat("*********                                        ********************\n")
cat("*********************************************************************\n")
cat("*********************************************************************\n\n")


tmpdf <- data.frame(means, medians, mads)
colnames(tmpdf) <- c("Mean", "Median", "MAD")
tmphtml <- html(tmpdf, file = "stats.before.centering.html", dec = 4, first.col = "Array name")


cat("\n\n Means, medians, MAD of log ratios  per subject/array\n")
round(data.frame(tmpdf), 5)

## Same, by chromosome:
## means
cat("\n\n Means per subject and chromosome\n")
a1 <- apply(xdata.merge1, 2,
            function(x) {
                unlist(tapply(x, positions.merge1$chrom.numeric,
                              function(z) {mean(z)}))
            }
            )
round(a1, 3)
cat("\n\n Medians per subject and chromosome\n")
### medians
a2 <- apply(xdata.merge1, 2,
            function(x) {
                unlist(tapply(x, positions.merge1$chrom.numeric,
                              function(z) {median(z)}))
            }
            )
round(a2, 3)
cat("\n\n MADs per subject and chromosome\n")
### madss
a3 <- apply(xdata.merge1, 2,
            function(x) {
                unlist(tapply(x, positions.merge1$chrom.numeric,
                              function(z) {mad(z)}))
            }
            )

round(a3, 3)

tmphtml <- html(a1, file = "stats.subj.by.chrom.mean.BEFORE.html", dec = 4,
                first.col = "Chromosome")
tmphtml <- html(a2, file = "stats.subj.by.chrom.median.BEFORE.html", dec = 4,
                first.col = "Chromosome")
tmphtml <- html(a3, file = "stats.subj.by.chrom.mad.BEFORE.html", dec = 4,
                first.col = "Chromosome")



#### Centering data
if(centering == "Mean") {
    xcenter <- xdata.merge1
    for(i in 1:ncol(xcenter))
        xcenter[, i] <- xdata.merge1[, i] - means[i]
} else if(centering == "Median") {
    xcenter <- xdata.merge1
    for(i in 1:ncol(xcenter))
        xcenter[, i] <- xdata.merge1[, i] - medians[i]
} else {
    xcenter <- xdata.merge1
}
## in case a single-array object
if(is.null(dim(xcenter))) xcenter <- matrix(xcenter, ncol = 1)
    
means <- apply(xcenter, 2, mean)
medians <- apply(xcenter, 2, median)
mads <- apply(xcenter, 2, mad)


cat("\n\n\n*********************************************************************\n")
cat("*********************************************************************\n")
cat("*********                                        ********************\n")
cat("*********       Basic statistics by array        ********************\n")
cat("*********        and array by chromosome         ********************\n")
cat("*********           AFTER  CENTERING            ********************\n")
cat("*********                                        ********************\n")
cat("*********************************************************************\n")
cat("*********************************************************************\n\n")

tmpdf <- data.frame(means, medians, mads)
colnames(tmpdf) <- c("Mean", "Median", "MAD")
tmphtml <- html(tmpdf, file = "stats.after.centering.html", dec = 4,
                first.col = "Array name")


cat("\n\n Means, medians, MAD of log ratios  per subject/array\n")
round(data.frame(tmpdf), 5)


## Same, by chromosome:

### means

cat("\n\n Means per subject and chromosome\n")
a1 <- apply(xcenter, 2,
            function(x) {
                unlist(tapply(x, positions.merge1$chrom.numeric,
                              function(z) {mean(z)}))
            })
           
round(a1, 3)

cat("\n\n Medians per subject and chromosome\n")
### medians
a2 <- apply(xcenter, 2,
            function(x) {
                unlist(tapply(x, positions.merge1$chrom.numeric,
                              function(z) {median(z)}))
            }
            )
round(a2, 3)

cat("\n\n MADs per subject and chromosome\n")
### madss
a3 <- apply(xcenter, 2,
            function(x) {
                unlist(tapply(x, positions.merge1$chrom.numeric,
                              function(z) {mad(z)}))
            }
            )
round(a3, 3)


sink()

tmphtml <- html(a1, file = "stats.subj.by.chrom.mean.AFTER.html", dec = 4, first.col = "Chromosome")
tmphtml <- html(a2, file = "stats.subj.by.chrom.median.AFTER.html", dec = 4, first.col = "Chromosome")
tmphtml <- html(a3, file = "stats.subj.by.chrom.mad.AFTER.html", dec = 4, first.col = "Chromosome")


### Checking not weird data
if(min(xcenter) < -1000)
    caughtUserError("The minimum value is < -1000; \n there is likely an error with the data (this would mean \n that the log2 ratio is < -1000 !!!)")

if(max(xcenter) > 1000)
    caughtUserError("The maximum value is > 1000; \n there is likely an error with the data (this would mean \n that the log2 ratio is > 1000 !!!)")

if(quantile(xcenter, 0.25) < -300)
    caughtUserError("At least 25 % of your vales are below < -300; \n there is likely an error with the data (this would mean \n that you have many values with log2 ratio is < -300 !!!)")

if(quantile(xcenter, 0.75) > 300)
    caughtUserError("At least 25 % of your vales are larger > 300; \n there is likely an error with the data (this would mean \n that you have many values with log2 ratio >  300 !!!)")


if(length(warningsForUsers)) {
    sink(file = "results.txt", append = TRUE)
cat("\n\n\n*********************************************************************\n")
cat("*********************************************************************\n")
cat("*********                                        ********************\n")
cat("*********    WARNINGS IN THE INTIAL PROCESSING   ********************\n")
cat("*********     Please make sure you are OK        ********************\n")
cat("*********           with all this               ********************\n")
cat("*********                                        ********************\n")
cat("*********************************************************************\n")
cat("*********************************************************************\n\n")
    
    print(warningsForUsers)
    sink()


    sink(file = "results.for.html", append = TRUE)
    cat("<h2> WARNINGS </h2>\n")
    cat("<pre>\n")
    print(warningsForUsers)
    cat("</pre>\n")
    sink()
}

#####################################################################
#####################################################################


ymax <- max(as.matrix(xcenter))
ymin <- min(as.matrix(xcenter))

numarrays <- ncol(xcenter)

options(warn = -1)

if(methodaCGH == "CBS") {
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
    
    ADaCGH:::mpiCBS()

    trythis <- try(
                   CNA.object <- CNA(as.matrix(xcenter), chrom = positions.merge1$chrom.numeric,
                                     ##                  maploc = positions.merge1$MidPoint,
                                     maploc = 1:length(positions.merge1$MidPoint),
                                     ## plotting is simpler
                                     data.type = "logratio",
                                     sampleid = colnames(xcenter))
                   )
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function CNA bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))
        
    
    trythis <-
        try(
            smoothed.CNA.object <- smooth.CNA(CNA.object,
                                              smooth.region = DNA.smooth.region,
                                              outlier.SD.scale = DNA.outlier.SD.scale,
                                              smooth.SD.scale = DNA.smooth.SD.scale,
                                              trim = DNA.trim)
            )
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function smooth.CNA bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))



    trythis <- try({
## FIXME: eventually, use the DNAcopyDiagnosticPlots function
                                        #    pdf("CBS.diagnostic.plots.%03d.pdf", height = 12, width = 18)
        pdf("CBS.diagnostic.plots.pdf", height = 12, width = 18)
        par(mfrow = c(4, 6))  ## zz: change this to a smaller number when few arrays.
        par(pty = "s")
        ## diagnostic plots for DNA copy; think about output zz
        for(i in 1:numarrays)
            plot(CNA.object[, (i + 2)], smoothed.CNA.object[, (i + 2)],
                 main = colnames(CNA.object[i + 2]),
                 xlab = "Original data", ylab = "Smoothed data")
        dev.off()
        
        ## so we will probably run into problems of huge pdf files. will deal with
        ## it later. Can use
        ## first, creat multiple pdf files.
        ##    tmp.npl <- ceiling(numarrays/24)
        ##   convert -quality 70 Gains.V2.pdf g2.png
        ## and insert those into an html
        
        pdf("CBS.diagnostic.plots.per.array.and.chromosome.pdf",
            height = 12, width = 18)
        par(pty = "s")
        for(i in 1:numarrays) {
            par(mfrow = c(4, 6))
            ncr <- unique(positions.merge1$chrom.numeric)
            for(j in ncr) {
                x <- CNA.object[positions.merge1$chrom.numeric == j, (i + 2)]
                y <- smoothed.CNA.object[positions.merge1$chrom.numeric == j, (i + 2)]
                plot(x, y,
                     main = paste(colnames(CNA.object[i+2]), "; Chr ", j, sep = ""),
                     xlab = "Original data", ylab = "Smoothed data")
            }
        }
        dev.off()
    })
#    if((!is.null(class(trythis))) & (class(trythis == "try-error")))
    if(class(trythis) == "try-error")
        caughtOurError(paste("Error in diagnostic plots with error",
                             trythis, ". \n Please let us know so we can fix the code."))


    
    ## p.method= we use hybrid
    trythis <-
        try({
            if(numarrays == 1) {
                segment.smoothed.CNA.object <- segment(smoothed.CNA.object,
                                                       alpha = DNA.copy.alpha,
                                                       kmax = DNA.kmax,
                                                       nmin = DNA.nmin,
                                                       nperm = DNA.nperm,
                                                       overlap = DNA.overlap,
                                                       trim = DNA.trim,
                                                       undo.splits = DNA.undo.splits,
                                                       undo.prune = DNA.undo.prune,
                                                       undo.SD = DNA.undo.sd,
                                                       p.method = "hybrid",
                                                       verbose = 2)
            } else {
                segment.smoothed.CNA.object <-
                    pSegmentDNAcopy(smoothed.CNA.object,
                                    alpha = DNA.copy.alpha,
                                    kmax = DNA.kmax,
                                    nmin = DNA.nmin,
                                    nperm = DNA.nperm,
                                    overlap = DNA.overlap,
                                    trim = DNA.trim,
                                    undo.splits = DNA.undo.splits,
                                    undo.prune = DNA.undo.prune,
                                    undo.SD = DNA.undo.sd,
                                    p.method = "hybrid")
##                                    verbose = 2)
            }
            })
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function segment (in CBS) bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))
save.image()                             
    ### Minimal common regions
    if(numarrays > 1) {
        ## cghMCR seems needs different classes.
        res <-  segment.smoothed.CNA.object
        class(res[[1]]) <- "data.frame"
        cghmcr <- cghMCR(res,
                         gapAllowed = MCR.gapAllowed,
                         alteredLow = MCR.alteredLow,
                         alteredHigh = MCR.alteredHigh)
        mcrs <- MCR(cghmcr)
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


        sink(file = "mcr.results.html")
        ##cat("<h3>Minimal common regions</h3>\n")
##         cat("<p>(Yes, this output might be ugly. We want your comments on how to",
##             "make the output more useful to you.)</p>\n")
        if (nrow(mcrsc) == 0) {
          cat("\n<p> No common minimal regions found.</p>\n")
        } else {
          html.data.frame(mcrsc, first.col = "Case",
                          file = "mcr.results.html", append = TRUE)
        }
        sink()
        
        sink(file = "results.txt")
        cat("\n\n\nMinimal common regions\n")
##         cat("(Yes, this output is ugly. We want your comments on how to",
##             "make the output more useful for you.)\n")
        if (nrow(mcrsc) == 0)
          cat("\n No common minimal regions found.\n")
        else 
          print(mcrsc)
        sink()
      }


    if(DNA.merge == "No") {
        trythis <- try({
            ## The segmented plots, one per array
            segmentPlot(segment.smoothed.CNA.object,
                        arraynames = colnames(xcenter),
                        chrom.numeric = positions.merge1$chrom.numeric,
                        idtype = idtype,
                        organism = organism,
                        geneNames = positions.merge1$name,
                        yminmax = c(ymin, ymax),
                        superimposed = FALSE)
            ## Supperimposed
            segmentPlot(segment.smoothed.CNA.object,
                        arraynames = colnames(xcenter),
                        chrom.numeric = positions.merge1$chrom.numeric,
                        idtype = idtype,
                        organism = organism,
                        geneNames = positions.merge1$name,
                        yminmax = c(ymin, ymax),
                        superimposed = TRUE)
            
        })
        if(class(trythis) == "try-error")
            caughtOurError(paste("Error in segment plots  with error",
                                 trythis, ". \n Please let us know so we can fix the code."))

############ Plateau plots
        
        trythis <- try({ ## fix this later, using GDD
            pdf("CBS.plateau.plots.pdf", height = 6, width = 9)
            plateauPlot(segment.smoothed.CNA.object)
            dev.off()
        })
    } else { ## If we used mergeLevels, then we can pretend we are using ACE for plotting
        ## there is some data duplication here.
        merged_segments <- mergeDNAcopy(segment.smoothed.CNA.object)
        trythis <- try({
          save.image()
          segmentPlot(merged_segments,
                      arraynames = colnames(xcenter),
                      geneNames = positions.merge1$name,
                      idtype = idtype,
                      organism = organism,
                      yminmax = c(ymin, ymax),
                      superimposed = FALSE)
          segmentPlot(merged_segments,
                      arraynames = colnames(xcenter),
                      geneNames = positions.merge1$name,
                      idtype = idtype,
                      organism = organism,
                      yminmax = c(ymin, ymax),
                      superimposed = TRUE)
      })
      }
    
    if(class(trythis) == "try-error")
        caughtOurError(paste("Error in plateau plots with error",
                             trythis, ". \n Please let us know so we can fix the code."))

    ## return final results
    merged_segments <- if(exists("merged_segments")) merged_segments else NULL
    writeResults(segment.smoothed.CNA.object, xcenter,
                 commondata = positions.merge1,
                 merged = merged_segments,
                 file = "CBS.results.txt")
    quit()
    
    
} else if(methodaCGH == "WS") {
    
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


    ADaCGH:::mpiWave()
    
    ## lets try the diagnostic plots code
    pdf("Autocorrelation.plots.pdf", width = 17.6, height = 12.5)
    par(mfrow = c(6,4))
    par(oma = c(2, 2, 2, 2))
    ADaCGH:::WaveletsDiagnosticPlots(xcenter,
                                     chrom.numeric = positions.merge1$chrom.numeric)
    dev.off()
    
    trythis <- try(
                   out.wave <-
                   pSegmentWavelets(as.matrix(xcenter),
                                    chrom.numeric = positions.merge1$chrom.numeric,
                                    minDiff = Wave.minDiff)
                   )
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function wave.aCGH bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))

    trythis <- try({
        segmentPlot(out.wave,
                    chrom.numeric = positions.merge1$chrom.numeric,
                    cghdata = xcenter,
                    arraynames = colnames(xcenter),
                    geneNames = positions.merge1$name,
                    idtype = idtype,
                    organism = organism,
                    yminmax = c(ymin, ymax),
                    superimposed = FALSE)
        segmentPlot(out.wave,
                    chrom.numeric = positions.merge1$chrom.numeric,
                    cghdata = xcenter,
                    arraynames = colnames(xcenter),
                    geneNames = positions.merge1$name,
                    idtype = idtype,
                    organism = organism,
                    yminmax = c(ymin, ymax),
                    superimposed = TRUE)
    })
    if(class(trythis) == "try-error")
        caughtOurError(paste("Error in segment plots  with error",
                             trythis, ". \n Please let us know so we can fix the code."))
    
    trythis <- try(
                   writeResults(out.wave, xcenter, commondata = positions.merge1,
                                file = "Wavelets.results.txt")
                   )
    if(class(trythis) == "try-error")
            caughtOurError(paste("Function print.wavelets.results bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))
    
    
#### The plateau plots:

    trythis <- try({ ## we should change this to use GDD, etc. Not for now.
        pdf("WS.plateau.plots.pdf", height = 6, width = 9)
        plateauPlot(out.wave, xcenter)
        dev.off()
    })
    if(class(trythis) == "try-error")
            caughtOurError(paste("An error printing the plateau plots with error",
                             trythis, ". \n Please let us know so we can fix the code."))
    ##save(file = "WS", list = ls())
    quit()
    
} else if(methodaCGH == "PSW") {
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

    ADaCGH:::mpiPSW()

print("testing existence of indicator")
print(exists(".__ADaCGH_WEB_APPL"))


    
## zz: PSW.nIter <- 1000; PSW.prec <- 100; PSW.p.crit <- 0.1
## p.crit is the largest p-value for
    ## which we want a region to be shown, in red,
    ## in the plot.

    out.gains <- data.frame(ID = positions.merge1$name,
                            Chromosome = positions.merge1$chromosome,
                            Start = positions.merge1$start,
                            End = positions.merge1$end,
                            MidPoint = positions.merge1$MidPoint)
    
    out.losses <- data.frame(ID = positions.merge1$name,
                             Chromosome = positions.merge1$chromosome,
                             Start = positions.merge1$start,
                             End = positions.merge1$end,
                             MidPoint = positions.merge1$MidPoint)


    print("testing existence of indicator before gains")
    print(exists(".__ADaCGH_WEB_APPL"))

### Gains
    trythis <- try({
        out.gains <-
            pSegmentPSW(out.gains,
                        xcenter,
                        chrom.numeric =  positions.merge1$chrom.numeric,
                        sign = +1, p.crit = PSW.p.crit,
                        nIter = PSW.nIter,
                        prec = PSW.prec,
                        name = "Gains.")
        segmentPlot(out.gains, geneNames = positions.merge1$name,
                    cghdata = xcenter,
                    idtype = idtype, organism = organism)
    })
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function pSegmentPSW (positive) bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))
    writeResults(out.gains, file = "Gains.Price.Smith.Waterman.results.txt")

    print("testing existence of indicator before losses")
    print(exists(".__ADaCGH_WEB_APPL"))

### Losses
    trythis <- try({
        out.losses <-
            pSegmentPSW(out.losses,
                        xcenter,
                        chrom.numeric =  positions.merge1$chrom.numeric,
                        sign = -1, p.crit = PSW.p.crit,
                        nIter = PSW.nIter,
                        prec = PSW.prec,
                        name = "Losses.")
        segmentPlot(out.losses, geneNames = positions.merge1$name,
                    cghdata = xcenter,
                    idtype = idtype, organism = organism)
    })
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function pSegmentPSW (negative) bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))
    writeResults(out.losses, file = "Losses.Price.Smith.Waterman.results.txt")
    
    save(file = "PSW.RData", list = ls(all.names = TRUE))
    ADaCGH:::PSWtoPaLS()
    
    quit()
    
    
} else if(methodaCGH == "ACE") {
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




  ADaCGH:::mpiACE()
    
    ## zz: ugly hack: it it is a 1 dimension array, make it a vector
    ## so that the correct methods are used.

    if(dim(xcenter)[2] == 1) {
        one.name <- colnames(xcenter)[1]
        xcenter <- xcenter[, 1]
    }
    
    trythis <- try(
                   ACE.object <- pSegmentACE(as.matrix(xcenter),
                                     chrom.numeric = positions.merge1$chrom.numeric)
                   )
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function pSegmentACE bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))
        

    trythis <- try(
                   ACE.summ <- summary(ACE.object, fdr = ACE.fdr)
                   )
    if(class(trythis) == "try-error")
        caughtOurError(paste("Function summary.ACE bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))

    save(file = "ace.RData", list = ls())

  
    trythis <- try(
                   writeResults(ACE.summ, commondata = positions.merge1,
                                file = NULL)
                   )
    if(class(trythis) == "try-error")
            caughtOurError(paste("Function writeResults.CGH.ACE.summary bombed unexpectedly with error",
                             trythis, ". \n Please let us know so we can fix the code."))


    ## re-hack. or re-do the kuldge:
    if(is.null(dim(xcenter))) {
        xcenter <- matrix(xcenter, ncol = 1)
        colnames(xcenter) <- one.name
    }

    save(file = "ace.RData", list = ls())
  
    trythis <- try({
        ## The segmented plots, one per array
        segmentPlot(ACE.summ,
                    chrom.numeric = positions.merge1$chrom.numeric,
                    geneNames = positions.merge1$name,
                    cghdata = xcenter,
                    idtype = idtype, organism = organism)
        ## Supperimposed
        segmentPlot(ACE.summ,
                    chrom.numeric = positions.merge1$chrom.numeric,
                    geneNames = positions.merge1$name,
                    cghdata = xcenter,
                    idtype = idtype, organism = organism,
                    superimposed = TRUE)
    })
    if(class(trythis) == "try-error")
        caughtOurError(paste("Error in ACE plots  with error",
                             trythis, ". \n Please let us know so we can fix the code."))
    save(file = "ace.RData", list = ls())
    quit()
}
