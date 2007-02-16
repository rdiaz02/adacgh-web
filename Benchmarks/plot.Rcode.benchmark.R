load("smallTiming10_all.RData")
load("smallTiming30_all.RData")
load("smallTiming60_all.RData")
load("smallTiming120_all.RData")
load("smallTimingNone_all.RData")


px <- function(x, method) {
    x <- x[x$Method == method,]
    tmp <- tapply(x$out, x$NumberArrays, mean)
    return(cbind(as.integer(names(tmp)),
                 tmp))
}  



pdf(file = "bench_r-code_small.pdf", height = 9, width = 12,
    onefile = FALSE, paper = "special")
ylim <- c(1, 290)
xlim <- c(2, 165)
par(mfrow = c(1, 3))
par(oma = c(0.5, 4, 2 , 3))
par(las = 1)
par(cex = 1.2)
par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "HMM"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "User wall time (seconds)",
     main = "HMM",
     xlim = xlim, axes = FALSE)
axis(1, at = c(5, 10, 20, 50, 100, 150),
     labels = TRUE)
axis(2, at = c(1, 2, 5, 10, 20, 50, 100, 150, 200, 300))
box()
points(px(smallTiming120, "HMM"), type = "b", col = "blue",  lwd = 2)
points(px(smallTiming60, "HMM"), type = "b", col = "green",  lwd = 2)
points(px(smallTiming30, "HMM"), type = "b", col = "orange",  lwd = 2)
points(px(smallTiming10, "HMM"), type = "b", col = "red",  lwd = 2)
text(50, 250, "Sequential code")
text(80, 30, "Parallelized code")
text(cbind(4, 0) + px(smallTiming120, "HMM")[6, ], "120", col = "blue", adj = 0)
text(cbind(4, 0) + px(smallTiming60, "HMM")[6, ], "60", col = "green", adj = 0)
text(cbind(4, 0) + px(smallTiming30, "HMM")[6, ], "30", col = "orange", adj = 0)
text(cbind(4, 0) + px(smallTiming10, "HMM")[6, ], "10", col = "red", adj = 0)


par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "GLAD"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(5, 10, 20, 50, 100, 150), 
     labels = TRUE)
points(px(smallTiming120, "GLAD"), type = "b", col = "blue",  lwd = 2)
points(px(smallTiming60, "GLAD"), type = "b", col = "green",  lwd = 2)
points(px(smallTiming30, "GLAD"), type = "b", col = "orange",  lwd = 2)
points(px(smallTiming10, "GLAD"), type = "b", col = "red",  lwd = 2)

par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "CBS"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "CBS",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(5, 10, 20, 50, 100, 150),
     labels = TRUE)
axis(4, at = c(1, 2, 5, 10, 20, 50, 100, 150, 200, 300))
points(px(smallTiming120, "CBS"), type = "b", col = "blue",  lwd = 2)
points(px(smallTiming60, "CBS"), type = "b", col = "green",  lwd = 2)
points(px(smallTiming30, "CBS"), type = "b", col = "orange",  lwd = 2)
points(px(smallTiming10, "CBS"), type = "b", col = "red",  lwd = 2)

mtext("Small data set (2271 genes)", side = 3, outer = TRUE, line = 0.2, cex = 1.5)
mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.5, cex = 1.5)
dev.off()






load("mediumTiming10_all.RData")
load("mediumTiming30_all.RData")
load("mediumTiming60_all.RData")
load("mediumTiming120_all.RData")
load("mediumTimingNone_all.RData")


pdf(file = "bench_r-code_medium.pdf", height = 9, width = 12,
    onefile = FALSE, paper = "special")
ylim <- c(5, 900)
xlim <- c(2, 165)
par(mfrow = c(1, 3))
par(oma = c(0.5, 4, 2 , 3))
par(las = 1)
par(cex = 1.2)
par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "HMM"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "User wall time (seconds)",
     main = "HMM",
     xlim = xlim, axes = FALSE)
axis(1, at = c(5, 10, 20, 50, 100, 150),
     labels = TRUE)
axis(2, at = c(5, 10, 20, 50, 100, 150, 200, 300, 500, 800))
box()
points(px(mediumTiming120, "HMM"), type = "b", col = "blue",  lwd = 2)
points(px(mediumTiming60, "HMM"), type = "b", col = "green",  lwd = 2)
points(px(mediumTiming30, "HMM"), type = "b", col = "orange",  lwd = 2)
points(px(mediumTiming10, "HMM"), type = "b", col = "red",  lwd = 2)
text(50, 500, "Sequential code")
text(80, 80, "Parallelized code")
text(cbind(4, 0) + px(mediumTiming120, "HMM")[6, ], "120", col = "blue", adj = 0)
text(cbind(4, 0) + px(mediumTiming60, "HMM")[6, ], "60", col = "green", adj = 0)
text(cbind(4, 0) + px(mediumTiming30, "HMM")[6, ], "30", col = "orange", adj = 0)
text(cbind(4, 0) + px(mediumTiming10, "HMM")[6, ], "10", col = "red", adj = 0)


par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "GLAD"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(5, 10, 20, 50, 100, 150), 
     labels = TRUE)
points(px(mediumTiming120, "GLAD"), type = "b", col = "blue",  lwd = 2)
points(px(mediumTiming60, "GLAD"), type = "b", col = "green",  lwd = 2)
points(px(mediumTiming30, "GLAD"), type = "b", col = "orange",  lwd = 2)
points(px(mediumTiming10, "GLAD"), type = "b", col = "red",  lwd = 2)

par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "CBS"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "CBS",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(5, 10, 20, 50, 100, 150),
     labels = TRUE)
axis(4, at = c(5, 10, 20, 50, 100, 150, 200, 300, 500, 800))
points(px(mediumTiming120, "CBS"), type = "b", col = "blue",  lwd = 2)
points(px(mediumTiming60, "CBS"), type = "b", col = "green",  lwd = 2)
points(px(mediumTiming30, "CBS"), type = "b", col = "orange",  lwd = 2)
points(px(mediumTiming10, "CBS"), type = "b", col = "red",  lwd = 2)

mtext("Medium data set (10,000 genes)", side = 3, outer = TRUE, line = 0.2, cex = 1.5)
mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.5, cex = 1.5)
dev.off()
