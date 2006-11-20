load("ace.RData")
#load(".RData")
.Last <- function(){
    save.image()
    cat("\n\n Normal ACE termination\n")
    ## In case the CGI is not called (user kills browser)
    ## have a way to stop lam
}
library(CGIwithR)
library(imagemap)
library(GDD)


graphDir <- paste(getwd(), "/", sep = "")
png.width = 7
png.height = 7
png.pointsize = 14
png.family = "Helvetica"


pid <- Sys.getpid()
write.table(file = "ACEpid.txt", pid,
            row.names = FALSE,
            col.names = FALSE)

library(ADaCGH2)

caughtUserError <- function(message) {
    webPNG("ErrorFigure.png", width = png.width,
           height = png.height, 
           pointsize = png.pointsize,
           family = png.family)
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


caughtOurError <- function(message) {
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
}

ACE.fdr <- scan("fdrace")


trythis <- try(
               ACE.summ <- summary(ACE.object, fdr = ACE.fdr)
               )
if(class(trythis) == "try-error")
    caughtOurError(paste("Function ACE.summary bombed unexpectedly with error",
                         trythis, ". \n Please let us know so we can fix the code."))


trythis <- try(
               print.ACE.results(ACE.summ)
               )
if(class(trythis) == "try-error")
    caughtOurError(paste("Function print.ACE.results bombed unexpectedly with error",
                         trythis, ". \n Please let us know so we can fix the code."))



trythis <- try({
    ## The segmented plots, one per array
    for(i in 1:numarrays) {
        plot.ace2(ACE.summ, positions.merge1$chrom.numeric, arraynum = i,
                  main = colnames(xcenter)[i])
    }
    
    ## Supperimposed
    plot.ace3(ACE.summ, positions.merge1$chrom.numeric, 
              main = "All_arrays",
              ylim = c(ymin, ymax),
              pch = "")
    plot.ace4(ACE.summ, positions.merge1$chrom.numeric, 
              main = "All_arrays",
              ylim = c(ymin, ymax))
    
})
if(class(trythis) == "try-error")
    caughtOurError(paste("Error in ACE plots  with error",
                         trythis, ". \n Please let us know so we can fix the code."))

