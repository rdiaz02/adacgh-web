load("smallTiming10.RData")
load("smallTiming30.RData")
load("smallTiming60.RData")
load("smallTiming120.RData")
load("smallTimingNone.RData")


px <- function(x, method) {
    x <- x[x$Method == method,]
    tmp <- tapply(x$out, x$NumberArrays, mean)
    return(cbind(as.integer(names(tmp)),
                 tmp))
}  



## ##eh, que esto es por método!!

## postscript(file = "bench.r-code.eps", height = 8, width = 12,
##            horizontal = FALSE,
##            onefile = FALSE, paper = "special")
## par(mfrow = c(1, 3))
## par(las = 1)
## par(cex = 1.2)

## plot(px(smallTimingNone, "HMM"),
##      ylim = c(2, 80), lwd = 2, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
## ##     xlab = "Number of arrays (samples)",
##      ylab = "User wall time (seconds)",
##      main = "HMM",
##      xlim = c(2, 58))
## axis(1, at = c(5, 10, 20, 50),
##      labels = TRUE)
## points(px(smallTiming120, "HMM"), type = "b", col = "blue",  lwd = 2)
## points(px(smallTiming60, "HMM"), type = "b", col = "green",  lwd = 2)
## points(px(smallTiming30, "HMM"), type = "b", col = "orange",  lwd = 2)
## points(px(smallTiming10, "HMM"), type = "b", col = "red",  lwd = 2)


## text(30, 70, "Sequential code")
## text(40, 10, "Parallelized code")
## text(cbind(2, 0) + px(smallTiming120, "HMM")[4, ], "120", col = "blue", adj = 0)
## text(cbind(2, 0) + px(smallTiming60, "HMM")[4, ], "60", col = "green", adj = 0)
## text(cbind(2, 0) + px(smallTiming30, "HMM")[4, ], "30", col = "orange", adj = 0)
## text(cbind(2, 0) + px(smallTiming10, "HMM")[4, ], "10", col = "red", adj = 0)





## plot(px(smallTimingNone, "GLAD"),
##      ylim = c(1, 50), lwd = 2, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "",
##      main = "GLAD",
##      xlim = c(2, 58))
## axis(1, at = c(5, 10, 20, 50),
##      labels = TRUE)
## points(px(smallTiming120, "GLAD"), type = "b", col = "blue",  lwd = 2)
## points(px(smallTiming60, "GLAD"), type = "b", col = "green",  lwd = 2)
## points(px(smallTiming30, "GLAD"), type = "b", col = "orange",  lwd = 2)
## points(px(smallTiming10, "GLAD"), type = "b", col = "red",  lwd = 2)


## plot(px(smallTimingNone, "CBS"),
##      ylim = c(3, 50), lwd = 2, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "",
##      main = "CBS",
##      xlim = c(2, 58))
## axis(1, at = c(5, 10, 20, 50),
##      labels = TRUE)
## points(px(smallTiming120, "CBS"), type = "b", col = "blue",  lwd = 2)
## points(px(smallTiming60, "CBS"), type = "b", col = "green",  lwd = 2)
## points(px(smallTiming30, "CBS"), type = "b", col = "orange",  lwd = 2)
## points(px(smallTiming10, "CBS"), type = "b", col = "red",  lwd = 2)

## dev.off()




## postscript(file = "bench.r-code.eps", height = 9, width = 12,
##            horizontal = FALSE,
##            onefile = FALSE, paper = "special")

pdf(file = "bench_r-code_small.pdf", height = 9, width = 12,
    onefile = FALSE, paper = "special")
ylim <- c(1, 80)
xlim <- c(2, 58)
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
axis(1, at = c(5, 10, 20, 50),
     labels = TRUE)
axis(2, at = c(1, 2, 5, 10, 20, 50))
box()
points(px(smallTiming120, "HMM"), type = "b", col = "blue",  lwd = 2)
points(px(smallTiming60, "HMM"), type = "b", col = "green",  lwd = 2)
points(px(smallTiming30, "HMM"), type = "b", col = "orange",  lwd = 2)
points(px(smallTiming10, "HMM"), type = "b", col = "red",  lwd = 2)
text(20, 70, "Sequential code")
text(30, 10, "Parallelized code")
text(cbind(2, 0) + px(smallTiming120, "HMM")[4, ], "120", col = "blue", adj = 0)
text(cbind(2, 0) + px(smallTiming60, "HMM")[4, ], "60", col = "green", adj = 0)
text(cbind(2, 0) + px(smallTiming30, "HMM")[4, ], "30", col = "orange", adj = 0)
text(cbind(2, 0) + px(smallTiming10, "HMM")[4, ], "10", col = "red", adj = 0)


par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "GLAD"),
     ylim = ylim, lwd = 2, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(5, 10, 20, 50),
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
axis(1, at = c(5, 10, 20, 50),
     labels = TRUE)
axis(4, at = c(1, 2, 5, 10, 20, 50))
points(px(smallTiming120, "CBS"), type = "b", col = "blue",  lwd = 2)
points(px(smallTiming60, "CBS"), type = "b", col = "green",  lwd = 2)
points(px(smallTiming30, "CBS"), type = "b", col = "orange",  lwd = 2)
points(px(smallTiming10, "CBS"), type = "b", col = "red",  lwd = 2)

mtext("Small data set (2271 genes)", side = 3, outer = TRUE, line = 0.2, cex = 1.5)
mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 1.5)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 2.5, cex = 1.5)
dev.off()




## smallTiming <- rbind(smallTimingNone, smallTiming10, smallTiming30,
##                      smallTiming60, smallTiming120)## yes, we get a warning from factors.
## smallTiming$Method <- as.factor(c(as.character(smallTimingNone$Method),
##                                   as.character(smallTiming10$Method),
##                                   as.character(smallTiming30$Method),
##                                   as.character(smallTiming60$Method),
##                                   as.character(smallTiming120$Method)))

## smallTiming$MPI <- factor(c(as.character(smallTimingNone$MPI),
##                     as.character(smallTiming10$MPI),
##                     as.character(smallTiming30$MPI),
##                     as.character(smallTiming60$MPI),
##                     as.character(smallTiming120$MPI)),
##                   levels = c("None", "10", "30", "60", "120"))

## uu <- tapply(smallTiming$out, list(smallTiming$MPI, smallTiming$Method), mean)
## vv <- uu
## dim(vv) <- NULL
## vv <- data.frame(out = vv)
## vv$MPI <- rep(rownames(uu), ncol(uu))
## vv$Metho

## xyplot(out ~ NumberArrays | Method, groups = MPI, data = smallTiming,
##        type = "b",
##              scales= list(y = list(at = c(1, 2, 5, 10, 20, 50), log = 10)),
##        auto.key = TRUE)


## xyplot(out ~ NumberArrays | Method, groups = MPI, data = smallTiming,
##        panel = "panel.superpose",
##        panel.groups = "panel.linejoin",
##              scales= list(y = list(at = c(1, 2, 5, 10, 20, 50), log = 10)),
##        auto.key = TRUE)

