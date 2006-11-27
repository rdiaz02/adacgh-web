\name{summary.ACE}
\alias{summary.ACE, summary.ACE.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.ACE(object, fdr = NULL, html = TRUE, outhtml = "ace.fdrtable.html")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{fdr}{ ~~Describe \code{fdr} here~~ }
  \item{html}{ ~~Describe \code{html} here~~ }
  \item{outhtml}{ ~~Describe \code{outhtml} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(object, fdr=NULL, html = TRUE,
                        outhtml ="ace.fdrtable.html") {

 nchrom <- length(object)
 FDR.table <- get.FDR(object, nchrom)

 if (is.null(fdr)) fdr <- 0.15

        index <- which.min(abs(FDR.table[, 3] - fdr))
 print(FDR.table)
        print(paste("Selected index", index))
        cat(FDR.table[index, 3], file ="aceFDR")

        .aceFDR.for.output <<- FDR.table[index, 3] ## xx: oh yes, ugly
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


        
##  ####Select a FDR level
##  if (is.null(index)) {
##   print(FDR.table)
##   index <- readline(paste("Please select an index (1 -", nrow(FDR.table), ") for FDR level:"))
##   if (!index \%in\% 1:nrow(FDR.table)) {
##    stop("Not a valid index")
##    }
##   index <- as.numeric(index)
##   }
  
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
 res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
