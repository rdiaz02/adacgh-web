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


#### For easier debugging: we go saving results as things go along.


## rm(list = ls()) ## Just in case.


checkpoint.num <- scan("checkpoint.num", what = double(0), n = 1)



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
##    try(system(paste("/http/mpi.log/killLAM.py", lamSESSION, "&")))
    try(mpi.quit(save = "yes"), silent = TRUE)
}

doCheckpoint <- function(num) {
    save.image()
    sink("checkpoint.num")
    cat(num)
    sink()
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

sink(file = "hostname")
cat(hostn)
sink()


### The above is no longer really needed. What follows is what buryThem2.py
##  uses
sink(file = "current_R_proc_info")
cat(hostn)
cat("   ")
cat(pid)
sink()



########################################################

########   Start MPI here to check if everything OK

#########################################################

methodaCGH <- scan("methodaCGH", what = "", n = 1)


assign(".__ADaCGH_WEB_APPL", TRUE)
print("testing existence of indicator")
print(exists(".__ADaCGH_WEB_APPL"))
library(Hmisc)
library("waveslim") ## we will have to load ADaCGH soon,
## but we must mask certain defs. in waveslim. So load
## waveslim here
library(ADaCGH)



## I am not sure this is really needed, since we do check this things elsewhere too,
## when we start the MPI universe from Python.

if (! ((methodaCGH == "PSW") & (checkpoint.num >= 4))) {
## we don't use MPI with PSW at the end
    library(Rmpi)
    
    MPI_MIN_UNIVERSE_SIZE <- 15
    
    if (mpi.universe.size () < MPI_MIN_UNIVERSE_SIZE) {
        cat("\n\n mpi.universe.size () < MPI_MIN_UNIVERSE_SIZE \n\n")
        quit(save = "no", status = 11, runLast = TRUE)
    }
    
    try({
        mpiInit()
        cat("\n\nAbout to print mpiOK file\n")
        sink(file = "mpiOK")
        cat("MPI started OK\n")
        sink()
    }
        )
} else { ## just because we need this for the controlling code
    sink(file = "mpiOK")
    cat("MPI started OK\n")
    sink()
}



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

if(checkpoint.num < 1) {

idtype <- try(scan("idtype", what = "", n = 1))
organism <- try(scan("organism", what = "", n = 1))

MCR.gapAllowed <- try(scan("MCR.gapAllowed", what = double(0), n = 1))
MCR.alteredLow <- try(scan("MCR.alteredLow", what = double(0), n = 1))
MCR.alteredHigh <- try(scan("MCR.alteredHigh", what = double(0), n = 1))
MCR.recurrence <- try(scan("MCR.recurrence", what = double(0), n = 1))



if (methodaCGH == "WS") {
    Wave.minDiff <-  scan("Wave.minDiff", what = double(0), n = 1)
    mergeRes <- Wave.merge            <- scan("Wave.merge", what = "", n = 1)          
} else if (methodaCGH == "PSW") {
    PSW.nIter <- scan("PSW.nIter", what = double(0), n = 1)
    PSW.p.crit <- scan("PSW.p.crit", what = double(0), n = 1)
} else if (methodaCGH == "ACE") {
    ACE.fdr <- scan("ACE.fdr", what = double(0), n = 1)
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

cat("<h4>Methods parameters </h4>\n")

if (methodaCGH == "WS") {
    cat("<p>        Wave.minDiff\t\t:          ",         Wave.minDiff,"</p>\n")
    cat("<p>        Wave.merge\t\t:            ",         Wave.merge,"</p>\n")
} else if (methodaCGH == "PSW") {
    cat("<p>        PSW.nIter\t\t:             ",             PSW.nIter,"</p>\n")
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


if (methodaCGH == "WS") {
    cat("\n\nWave.minDiff\t\t:          ",         Wave.minDiff)
    cat("\n\nWave.merge\t\t:             ",       Wave.merge)        
} else if (methodaCGH == "PSW") {
    cat("\n\nPSW.nIter\t\t:             ",             PSW.nIter)
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
tmp2 <- data.frame(as.vector(round(tmp, 0)))
rownames(tmp2) <- names(tmp)
colnames(tmp2) <- "Number of genes/clones per chromosome"
tmp2[, 1] <- as.character(round(tmp2[, 1]))
tmphtml <- ADaCGH:::my.html.data.frame(tmp2, file = "clones.per.chrom.html", dec = 0, first.col = "Chromosome")
                
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
tmphtml <- ADaCGH:::my.html.data.frame(tmpdf, file = "stats.before.centering.html", digits = 4, first.col = "Array name")


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

tmphtml <- ADaCGH:::my.html.data.frame(a1, file = "stats.subj.by.chrom.mean.BEFORE.html", digits = 4,
                first.col = "Chromosome")
tmphtml <- ADaCGH:::my.html.data.frame(a2, file = "stats.subj.by.chrom.median.BEFORE.html", digits = 4,
                first.col = "Chromosome")
tmphtml <- ADaCGH:::my.html.data.frame(a3, file = "stats.subj.by.chrom.mad.BEFORE.html", digits = 4,
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
tmphtml <- ADaCGH:::my.html.data.frame(tmpdf, file = "stats.after.centering.html", digits = 4,
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

tmphtml <- ADaCGH:::my.html.data.frame(a1, file = "stats.subj.by.chrom.mean.AFTER.html", digits = 4, first.col = "Chromosome")
tmphtml <- ADaCGH:::my.html.data.frame(a2, file = "stats.subj.by.chrom.median.AFTER.html", digits = 4, first.col = "Chromosome")
tmphtml <- ADaCGH:::my.html.data.frame(a3, file = "stats.subj.by.chrom.mad.AFTER.html", digits = 4, first.col = "Chromosome")


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

doCheckpoint(1)

}
#####################################################################
#####################################################################
options(warn = -1)

### This all there is to execution itself

  

if(! (methodaCGH %in% c("PSW", "ACE"))) {




    if(checkpoint.num < 2) {
        if(!(exists("mergeRes"))) mergeRes <- TRUE
        
        common.data <- data.frame(ID = positions.merge1$name,
                                  Chromosome = positions.merge1$chromosome,
                                  Start = positions.merge1$start,
                                  End = positions.merge1$end,
                                  MidPoint = positions.merge1$MidPoint)
        ##save.image()
        doCheckpoint(2)

    }
    if(checkpoint.num < 3) {

        SegmentPlotWrite(as.matrix(xcenter),
                         chrom = positions.merge1$chrom.numeric,
                         mergeSegs = mergeRes,
                         Pos = positions.merge1$MidPoint,
                         idtype = idtype,
                         organism = organism,
                         method = methodaCGH,
                         geneNames = positions.merge1$name,
                         commondata = common.data, ## zz?
                         MCR.gapAllowed = MCR.gapAllowed,
                         MCR.alteredLow = MCR.alteredLow,
                         MCR.alteredHigh = MCR.alteredHigh,
                         MCR.recurrence = MCR.recurrence)
        
##        save.image()
        doCheckpoint(3)
        quit()
    }
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

    if(checkpoint.num < 2) {
        
        common.data <- data.frame(ID = positions.merge1$name,
                                Chromosome = positions.merge1$chromosome,
                                Start = positions.merge1$start,
                                End = positions.merge1$end,
                                MidPoint = positions.merge1$MidPoint)
        print("testing existence of indicator before gains")
        print(exists(".__ADaCGH_WEB_APPL"))
        
        ## save.image()
        doCheckpoint(2)
    }

    if(checkpoint.num < 3) {

## change order: first segmentation of both. then plotting after a "mpi.exit()

        
        ## Gains
        trythis <- try({
            out.gains <-
                pSegmentPSW(xcenter,
                            chrom.numeric =  positions.merge1$chrom.numeric,
                            common.data = common.data,
                            sign = +1, p.crit = PSW.p.crit,
                            nIter = PSW.nIter,
                            prec = 100,
                            name = "Gains.")
            cat("\n ************ done segmentation positive \n")
            save.image()
            
        })
        if(class(trythis) == "try-error")
            caughtOurError(paste("Function pSegmentPSW (positive) bombed unexpectedly with error",
                                 trythis, ". \n Please let us know so we can fix the code."))
        writeResults(out.gains, acghdata = xcenter, commondata = common.data,
                     file = "Gains.Price.Smith.Waterman.results.txt")

        doCheckpoint(3)
    }
    if(checkpoint.num < 4) {
        
        print("testing existence of indicator before losses")
        print(exists(".__ADaCGH_WEB_APPL"))
        
        ## Losses
        trythis <- try({
            out.losses <-
                pSegmentPSW(xcenter,
                            chrom.numeric =  positions.merge1$chrom.numeric,
                            common.data = common.data,
                            sign = -1, p.crit = PSW.p.crit,
                            nIter = PSW.nIter,
                            prec = 100,
                            name = "Losses.")
            save.image()
            cat("\n ************ done segmentation negative \n")
            
            save.image()
            
        })
        if(class(trythis) == "try-error")
            caughtOurError(paste("Function pSegmentPSW (negative) bombed unexpectedly with error",
                                 trythis, ". \n Please let us know so we can fix the code."))
        writeResults(out.losses, acghdata = xcenter, commondata = common.data,
                     file = "Losses.Price.Smith.Waterman.results.txt")
        
        save(file = "PSW.RData", list = ls(all.names = TRUE))
        PSWtoPaLS()
        doCheckpoint(4)
    }
    if(checkpoint.num < 5) {
        try(mpi.exit()) ## now papply0 is called inside segmentPlot
        segmentPlot(out.gains, geneNames = positions.merge1$name,
                    cghdata = xcenter,
                    idtype = idtype, organism = organism)
        
        segmentPlot(out.losses, geneNames = positions.merge1$name,
                    cghdata = xcenter,
                    idtype = idtype, organism = organism)
        
        doCheckpoint(5)
        quit()
    }
    
} else if(methodaCGH == "ACE") {

    if(checkpoint.num < 2) {

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
        
        ##save.image()
        doCheckpoint(2)
    }

    if(checkpoint.num < 3) {
        
        trythis <- try(
                       ACE.summ <- summary(ACE.object, fdr = ACE.fdr)
                       )
        if(class(trythis) == "try-error")
            caughtOurError(paste("Function summary.ACE bombed unexpectedly with error",
                                 trythis, ". \n Please let us know so we can fix the code."))
        
        save(file = "ace.RData", list = ls())

        common.data <- data.frame(ID = positions.merge1$name,
                                  Chromosome = positions.merge1$chromosome,
                                  Start = positions.merge1$start,
                                  End = positions.merge1$end,
                                  MidPoint = positions.merge1$MidPoint)
       
        trythis <- try(
                       writeResults(ACE.summ,
                                    acghdata = as.matrix(xcenter),
                                    commondata = common.data,
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
        
        trythis <- try(doMCR(ACE.summ$segm,
                             chrom.numeric = positions.merge1$chrom.numeric,
                             data = xcenter,
                             MCR.gapAllowed = MCR.gapAllowed,
                             MCR.alteredLow = MCR.alteredLow,
                             MCR.alteredHigh = MCR.alteredHigh,
                             MCR.recurrence = MCR.recurrence,
                             Pos = positions.merge1$MidPoint)
                       )
        if(class(trythis) == "try-error")
            caughtOurError(trythis)

        save.image()
        
        trythis <- try({
            ## The segmented plots, one per array
            segmentPlot(ACE.summ,
                        chrom.numeric = positions.merge1$chrom.numeric,
                        geneNames = positions.merge1$name,
                        cghdata = xcenter,
                        idtype = idtype, organism = organism)
        })
        if(class(trythis) == "try-error")
            caughtOurError(paste("Error in ACE plots  with error",
                                 trythis, ". \n Please let us know so we can fix the code."))
        save(file = "ace.RData", list = ls())
        doCheckpoint(3)
        quit()
    }
}
