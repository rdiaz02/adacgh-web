library(lattice)
freadall7 <- function(name) {
   fns <- paste("web.bnchmk7.", name, "_",
                c("1.1", "1.2", "1.3", "1.4", "1.5",
                  "5.1", "10.1", "20.1"), ".txt",
                sep = "")
   tms <- unlist(mapply(function(x, y) scan(x, n = y, what = double(0)),
                 fns, c(rep(1, 5), 5, 10, 20)))
   names(tms) <- NULL
   return(tms)
}
   
nusers <- factor(c(rep(1, 5), rep(5, 5), rep(10, 10), rep(20, 20)))


cbs.s <- freadall7("CBS_small")
hmm.s <- freadall7("HMM_small")
glad.s <- freadall7("GLAD_small")
wavelets.s <- freadall7("Wavelets_small")
cghseg.s <- freadall7("CGHseg_small")

small.timings7 <- data.frame(timings = c(cbs.s, hmm.s, glad.s, cghseg.s),
                            nusers = rep(nusers, 4),
                            method = c(rep("CBS", 40),
                            rep("HMM", 40),
                            rep("GLAD", 40),
                            rep("CGHseg", 40)))


cbs.m <- freadall7("CBS_medium")
hmm.m <- freadall7("HMM_medium")
glad.m <- freadall7("GLAD_medium")
wavelets.m <- freadall7("Wavelets_medium")
cghseg.m <- freadall7("CGHseg_medium")

medium.timings7 <- data.frame(timings = c(cbs.m, hmm.m, glad.m, cghseg.m),
                            nusers = rep(nusers, 4),
                            method = c(rep("CBS", 40),
                            rep("HMM", 40),
                            rep("GLAD", 40),
                            rep("CGHseg", 40)))

cbs.m <- freadall7("CBS_large")
hmm.m <- freadall7("HMM_large")
glad.m <- freadall7("GLAD_large")
wavelets.m <- freadall7("Wavelets_large")
cghseg.m <- freadall7("CGHseg_large")

large.timings7 <- data.frame(timings = c(cbs.m, hmm.m, glad.m, cghseg.m),
                             nusers = rep(nusers, 4),
                             method = c(rep("CBS", 40),
                             rep("HMM", 40),
                             rep("GLAD", 40),
                             rep("CGHseg", 40)))






medium.timings <- medium.timings7
small.timings <- small.timings7
large.timings <- large.timings7



## freadall6 <- function(name) {
##    fns <- paste("web.bnchmk6.", name, "_",
##                 c("1.1", "1.2", "1.3", "1.4", "1.5",
##                   "5.1", "10.1", "20.1"), ".txt",
##                 sep = "")
##    tms <- unlist(mapply(function(x, y) scan(x, n = y, what = double(0)),
##                  fns, c(rep(1, 5), 5, 10, 20)))
##    names(tms) <- NULL
##    return(tms)
## }

## nusers <- factor(c(rep(1, 5), rep(5, 5), rep(10, 10), rep(20, 20)))


## cbs.s <- freadall6("CBS_small")
## hmm.s <- freadall6("HMM_small")
## glad.s <- freadall6("GLAD_small")
## wavelets.s <- freadall6("Wavelets_small")
## cghseg.s <- freadall6("CGHseg_small")

## small.timings6 <- data.frame(timings = c(cbs.s, hmm.s, glad.s, cghseg.s),
##                             nusers = rep(nusers, 4),
##                             method = c(rep("CBS", 40),
##                             rep("HMM", 40),
##                             rep("GLAD", 40),
##                             rep("CGHseg", 40)))


## cbs.m <- freadall6("CBS_medium")
## hmm.m <- freadall6("HMM_medium")
## glad.m <- freadall6("GLAD_medium")
## wavelets.m <- freadall6("Wavelets_medium")
## cghseg.m <- freadall6("CGHseg_medium")

## medium.timings6 <- data.frame(timings = c(cbs.m, hmm.m, glad.m, cghseg.m),
##                             nusers = rep(nusers, 4),
##                             method = c(rep("CBS", 40),
##                             rep("HMM", 40),
##                             rep("GLAD", 40),
##                             rep("CGHseg", 40)))




## small.timings <- rbind(small.timings6, small.timings7)
## medium.timings <- rbind(medium.timings6, medium.timings7)



library(lattice)
pdf(file = "bench_web_small.pdf", height = 9, width = 12,
           onefile = FALSE, paper = "special")
 par(new = TRUE)
 plot.new()
##  mtext("Small data set (2271 genes, 10 arrays)", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
##  mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 2, cex = 1.5)
##  par(las = 0)
##  mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.8, cex = 1.5)
##  par(new = TRUE)
par(oma = c(3, 3, 4, 2))
bwplot(timings ~ nusers | method, groups = nusers, data = small.timings,
       xlab = "", ylab = "",
#       ylim = c(50, 450),
       scales = list(cex = 1.2,
       y = list(log = TRUE,
       at = c(75, 100, 200, 400))),
       cex.lab = 1.5,
       panel = function(...) {
           panel.bwplot(fill = "gray", cex.lab = 1.5, ...)
       }
       )
par(new = TRUE)
mtext("2271 genes, 10 arrays", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 1.5, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 1, cex = 1.5)
dev.off()



pdf(file = "bench_web_medium.pdf", height = 9, width = 12,
           onefile = FALSE, paper = "special")
 par(new = TRUE)
 plot.new()
##  mtext("Small data set (2271 genes, 10 arrays)", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
##  mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 2, cex = 1.5)
##  par(las = 0)
##  mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.8, cex = 1.5)
##  par(new = TRUE)
par(oma = c(3, 3, 4, 2))
bwplot(timings ~ nusers | method, groups = nusers, data = medium.timings,
       xlab = "", ylab = "",
       scales = list(cex = 1.2,
       y = list(log = TRUE,
       at = c(250, 500, 1000, 2000, 4000))),
       cex.lab = 1.5,
       panel = function(...) {
           panel.bwplot(fill = "gray", cex.lab = 1.5, ...)
       }
       )
par(new = TRUE)
mtext("15000 genes, 40 arrays", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 1.5, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 1.5, cex = 1.5)
dev.off()


pdf(file = "bench_web_large.pdf", height = 9, width = 12,
           onefile = FALSE, paper = "special")
 par(new = TRUE)
 plot.new()
##  mtext("Small data set (2271 genes, 10 arrays)", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
##  mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 2, cex = 1.5)
##  par(las = 0)
##  mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.8, cex = 1.5)
##  par(new = TRUE)
par(oma = c(3, 3, 4, 2))
bwplot(timings ~ nusers | method, groups = nusers, data = large.timings,
       xlab = "", ylab = "",
       scales = list(cex = 1.2,
       y = list(log = TRUE,
       at = c(500, 1000, 2000, 5000, 10000, 20000))),
       cex.lab = 1.5,
       panel = function(...) {
           panel.bwplot(fill = "gray", cex.lab = 1.5, ...)
       }
       )
par(new = TRUE)
mtext("42325 genes, 40 arrays", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 1.5, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 1.5, cex = 1.5)
dev.off()
