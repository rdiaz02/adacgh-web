load("/home/ramon/bzr-local-repos/adacgh2/Benchmarks/smallTiming10.RData")
load("/home/ramon/bzr-local-repos/adacgh2/Benchmarks/smallTiming30.RData")
load("/home/ramon/bzr-local-repos/adacgh2/Benchmarks/smallTiming60.RData")
load("/home/ramon/bzr-local-repos/adacgh2/Benchmarks/smallTiming120.RData")
load("/home/ramon/bzr-local-repos/adacgh2/Benchmarks/smallTimingNone.RData")

smallTiming <- rbind(smallTimingNone, smallTiming10, smallTiming30,
                     smallTiming60, smallTiming120)## yes, we get a warning from factors.

summary(smallTiming) ##1.27 to 91.4

px <- function(x) {
    tmp <- tapply(x$out, x$NumberArrays, mean)
    return(cbind(as.integer(names(tmp)),
                 tmp))
}  



##eh, que esto es por método!!

postscript(file = "bench.r-code.eps", height = 8, width = 12,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")
par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.2)

plot(px(smallTimingNone),
     ylim = c(1, 70), lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "Number of arrays (samples)",
     ylab = "User wall time (seconds)",
     main = "Small data set (2271 genes)",
     xlim = c(1, 65))
axis(1, at = c(5, 10, 20, 50),
     labels = TRUE)
points(px(smallTiming120), type = "b", col = "blue",  lwd = 2)
points(px(smallTiming60), type = "b", col = "green",  lwd = 2)
points(px(smallTiming30), type = "b", col = "orange",  lwd = 2)
points(px(smallTiming10), type = "b", col = "red",  lwd = 2)


text(20, 50, "Original sequential code")
text(20, 5, "Parallelized code")
text(cbind(2, 0) + px(smallTiming120)[4, ], "120 slaves", col = "blue", adj = 0)
text(cbind(2, 0) + px(smallTiming60)[4, ], "60 slaves", col = "green", adj = 0)
text(cbind(2, 0) + px(smallTiming30)[4, ], "30 slaves", col = "orange", adj = 0)
text(cbind(2, 0) + px(smallTiming10)[4, ], "10 slaves", col = "red", adj = 0)


dev.off()
