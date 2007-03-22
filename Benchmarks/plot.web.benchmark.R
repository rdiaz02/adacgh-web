## Small data sets



freadall <- function(name) {
   fns <- paste("web.bnchmk2.", name, "_",
                c("1.1", "1.2", "1.3", "1.4", "1.5",
                  "2.1", "2.2", "2.3", "5.1"), ".txt",
                sep = "")
   tms <- unlist(mapply(function(x, y) scan(x, n = y, what = double(0)),
                 fns, c(rep(1, 5), rep(2, 3), 5)))
   names(tms) <- NULL
   return(tms)
}
   

cbs <- freadall("CBS_small")
hmm <- freadall("HMM_small")
glad <- freadall("GLAD_small")
wavelets <- freadall("Wavelets_small")

nusers <- factor(c(rep(1, 5), rep(2, 6), rep(5, 5)))
                 

small.timings <- data.frame(timings = c(cbs, hmm, glad, wavelets),
                            nusers = rep(nusers, 4),
                            method = c(rep("CBS", 16),
                            rep("HMM", 16),
                            rep("GLAD", 16),
                            rep("Wavelets", 16)))

pdf(file = "bench_web_small.pdf", height = 9, width = 12,
           onefile = FALSE, paper = "special")
 par(new = TRUE)
 plot.new()
##  mtext("Small data set (2271 genes, 10 arrays)", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
##  mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 2, cex = 1.5)
##  par(las = 0)
##  mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.8, cex = 1.5)
##  par(new = TRUE)
par(oma = c(4, 4, 4, 4))
bwplot(timings ~ nusers | method, groups = nusers, data = small.timings,
       xlab = "", ylab = "",
       scales = list(cex = 1.2),
       cex.lab = 1.5,
       panel = function(...) {
           panel.bwplot(fill = "gray", cex.lab = 1.5, ...)
       }
       )
par(new = TRUE)
mtext("Small data set (2271 genes, 10 arrays)", side = 3, outer = TRUE, line = 2.5, cex = 1.5)
mtext("Number of simultaneous users", side = 1, outer = TRUE, line = 2, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.8, cex = 1.5)
dev.off()

