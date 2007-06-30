load("smallTiming10_all.RData")
load("smallTiming30_all.RData")
load("smallTiming60_all.RData")
load("smallTiming120_all.RData")
load("smallTimingNone_all.RData")


load("timings.RData")


smallTimingNone <- smallTimings[smallTimings$MPI == "None", ]
smallTiming10 <- smallTimings[smallTimings$MPI == "10", ]
smallTiming30 <- smallTimings[smallTimings$MPI == "30", ]
smallTiming60 <- smallTimings[smallTimings$MPI == "60", ]


mediumTimingNone <- mediumTimings[mediumTimings$MPI == "None", ]
mediumTiming10 <- mediumTimings[mediumTimings$MPI == "10", ]
mediumTiming30 <- mediumTimings[mediumTimings$MPI == "30", ]
mediumTiming60 <- mediumTimings[mediumTimings$MPI == "60", ]

largeTimingNone <- largeTimings[largeTimings$MPI == "None", ]
largeTiming10 <- largeTimings[largeTimings$MPI == "10", ]
largeTiming30 <- largeTimings[largeTimings$MPI == "30", ]
largeTiming60 <- largeTimings[largeTimings$MPI == "60", ]




px <- function(x, method) {
    x <- x[x$Method == method,]
    tmp <- tapply(x$out, x$NumberArrays, mean)
    return(cbind(as.integer(names(tmp)),
                 tmp))
}  



postscript(file = "bench_r-code_small.eps", height = 9, width = 17,
    onefile = FALSE, paper = "special")
ylim <- c(4, 6500)
xlim <- c(10 , 165)
par(mfrow = c(1, 4))
par(oma = c(0, 7, 2 , 4))
par(las = 1)
par(cex = 1)
par(mar = c(5, 0, 4, 0))
par(cex.lab = 2)
par(cex.main = 2.2)
par(cex.axis = 1.6)
plot(px(smallTimingNone, "HMM"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "User wall time (seconds)",
     main = "HMM",
     xlim = xlim, axes = FALSE)
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
axis(1, at = c(20),
     labels = TRUE, cex.axis = 0.9)
axis(2, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000))
box()
points(px(smallTiming60, "HMM"), type = "b", col = "blue",  lwd = 1.5)
points(px(smallTiming30, "HMM"), type = "b", col = "orange",  lwd = 1.5)
points(px(smallTiming10, "HMM"), type = "b", col = "red",  lwd = 1.5)
text(60, 850, "Sequential code", cex = 1.8)
text(80, 90, "Parallelized code", cex = 1.8)
text(cbind(4, 0) + px(smallTiming60, "HMM")[5, ], "60", col = "blue", adj = 0, cex = 1.9)
text(cbind(4, 0) + px(smallTiming30, "HMM")[5, ], "30", col = "orange", adj = 0, cex = 1.9)
text(cbind(4, 0) + px(smallTiming10, "HMM")[5, ], "10", col = "red", adj = 0, cex = 1.9)


par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "GLAD"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
points(px(smallTiming60, "GLAD"), type = "b", col = "blue",  lwd = 1.5)
points(px(smallTiming30, "GLAD"), type = "b", col = "orange",  lwd = 1.5)
points(px(smallTiming10, "GLAD"), type = "b", col = "red",  lwd = 1.5)



par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "CBS"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "CBS",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
points(px(smallTiming60, "CBS"), type = "b", col = "blue",  lwd = 1.5)
points(px(smallTiming30, "CBS"), type = "b", col = "orange",  lwd = 1.5)
points(px(smallTiming10, "CBS"), type = "b", col = "red",  lwd = 1.5)


par(mar = c(5, 0, 4, 0))
plot(px(smallTimingNone, "BioHMM"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "BioHMM",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
axis(4, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000))
points(px(smallTiming60, "BioHMM"), type = "b", col = "blue",  lwd = 1.5)
points(px(smallTiming30, "BioHMM"), type = "b", col = "orange",  lwd = 1.5)
points(px(smallTiming10, "BioHMM"), type = "b", col = "red",  lwd = 1.5)

mtext("a) 10,000 genes", side = 3, outer = TRUE, line = 0.2, cex = 2.5, adj = 0)
## mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 2)
## par(las = 0)
## mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 3, cex = 2)
dev.off()






postscript(file = "bench_r-code_medium.eps", height = 9, width = 17,
    onefile = FALSE, paper = "special")
ylim <- c(15, 12000)
xlim <- c(10 , 165)
par(mfrow = c(1, 4))
par(oma = c(0.5, 7, 2 , 4))
par(las = 1)
par(cex = 1)
par(mar = c(5, 0, 4, 0))

par(cex.lab = 2)
par(cex.main = 2.2)
par(cex.axis = 1.6)

plot(px(mediumTimingNone, "HMM"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "User wall time (seconds)",
     main = "HMM",
     xlim = xlim, axes = FALSE)
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
axis(1, at = c(20),
     labels = TRUE, cex.axis = 0.7)
axis(2, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000, 10000))
box()
points(px(mediumTiming60, "HMM"), type = "b", col = "blue",  lwd = 1.5)
points(px(mediumTiming30, "HMM"), type = "b", col = "orange",  lwd = 1.5)
points(px(mediumTiming10, "HMM"), type = "b", col = "red",  lwd = 1.5)
## text(60, 1400, "Sequential code")
## text(80, 150, "Parallelized code")
## text(cbind(4, 0) + px(mediumTiming60, "HMM")[5, ], "60", col = "blue", adj = 0)
## text(cbind(4, 0) + px(mediumTiming30, "HMM")[5, ], "30", col = "orange", adj = 0)
## text(cbind(4, 0) + px(mediumTiming10, "HMM")[5, ], "10", col = "red", adj = 0)


par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "GLAD"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
points(px(mediumTiming60, "GLAD"), type = "b", col = "blue",  lwd = 1.5)
points(px(mediumTiming30, "GLAD"), type = "b", col = "orange",  lwd = 1.5)
points(px(mediumTiming10, "GLAD"), type = "b", col = "red",  lwd = 1.5)


par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "CBS"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "CBS",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
points(px(mediumTiming60, "CBS"), type = "b", col = "blue",  lwd = 1.5)
points(px(mediumTiming30, "CBS"), type = "b", col = "orange",  lwd = 1.5)
points(px(mediumTiming10, "CBS"), type = "b", col = "red",  lwd = 1.5)


par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "BioHMM"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "BioHMM",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
axis(4, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000, 10000))
points(px(mediumTiming60, "BioHMM"), type = "b", col = "blue",  lwd = 1.5)
points(px(mediumTiming30, "BioHMM"), type = "b", col = "orange",  lwd = 1.5)
points(px(mediumTiming10, "BioHMM"), type = "b", col = "red",  lwd = 1.5)

mtext("b) 20,000 genes", side = 3, outer = TRUE, line = 0.2, cex = 2.5,
      adj = 0)
## mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 2)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 5, cex = 2.5)
dev.off()






postscript(file = "bench_r-code_large.eps", height = 9, width = 17,
    onefile = FALSE, paper = "special")
ylim <- c(50, 25000)
xlim <- c(10 , 165)
par(mfrow = c(1, 4))
par(oma = c(1, 7, 2 , 4))
par(las = 1)
par(cex = 1)
par(mar = c(5, 0, 4, 0))

par(cex.lab = 2)
par(cex.main = 2.2)
par(cex.axis = 1.6)



plot(px(largeTimingNone, "HMM"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "User wall time (seconds)",
     main = "HMM",
     xlim = xlim, axes = FALSE)
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
axis(1, at = c(20),
      labels = TRUE, cex.axis = 0.7)
axis(2, at = c(5, 20, 50, 100, 200, 500, 1000, 2500, 5000, 10000, 25000))
box()
points(px(largeTiming60, "HMM"), type = "b", col = "blue",  lwd = 1.5)
points(px(largeTiming30, "HMM"), type = "b", col = "orange",  lwd = 1.5)
points(px(largeTiming10, "HMM"), type = "b", col = "red",  lwd = 1.5)
## text(60, 3500, "Sequential code")
## text(80, 400, "Parallelized code")
## text(cbind(4, 0) + px(largeTiming60, "HMM")[5, ], "60", col = "blue", adj = 0)
## text(cbind(4, 0) + px(largeTiming30, "HMM")[5, ], "30", col = "orange", adj = 0)
## text(cbind(4, 0) + px(largeTiming10, "HMM")[5, ], "10", col = "red", adj = 0)


par(mar = c(5, 0, 4, 0))
plot(px(largeTimingNone, "GLAD"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
points(px(largeTiming60, "GLAD"), type = "b", col = "blue",  lwd = 1.5)
points(px(largeTiming30, "GLAD"), type = "b", col = "orange",  lwd = 1.5)
points(px(largeTiming10, "GLAD"), type = "b", col = "red",  lwd = 1.5)



par(mar = c(5, 0, 4, 0))
plot(px(largeTimingNone, "CBS"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "CBS",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
points(px(largeTiming60, "CBS"), type = "b", col = "blue",  lwd = 1.5)
points(px(largeTiming30, "CBS"), type = "b", col = "orange",  lwd = 1.5)
points(px(largeTiming10, "CBS"), type = "b", col = "red",  lwd = 1.5)


par(mar = c(5, 0, 4, 0))
plot(px(largeTimingNone, "BioHMM"),
     ylim = ylim, lwd = 1.5, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "BioHMM",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
axis(4, at = c(5, 20, 50, 100, 200, 500, 1000, 2500, 5000, 10000, 25000))
points(px(largeTiming60, "BioHMM"), type = "b", col = "blue",  lwd = 1.5)
points(px(largeTiming30, "BioHMM"), type = "b", col = "orange",  lwd = 1.5)
points(px(largeTiming10, "BioHMM"), type = "b", col = "red",  lwd = 1.5)

mtext("c) 42,325 genes", side = 3, outer = TRUE, line = 0.2, cex = 2.5,
      adj = 0)
mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -1,
      cex = 2.5)
## par(las = 0)
## mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 3.5,
##       cex = 2)
dev.off()


system("epstopdf bench_r-code_small.eps")
system("epstopdf bench_r-code_medium.eps")
system("epstopdf bench_r-code_large.eps")
system("pdftk bench_r-code_small.pdf bench_r-code_medium.pdf bench_r-code_large.pdf cat output bench_preall.pdf")
system("pdfnup --orient landscape --nup 3x1 bench_preall.pdf --outfile bench_r_all.pdf")
system("pdf2ps bench_r_all.pdf")
system("cp bench_r_all.pdf ~/Proyectos/ADaCGH-paper/GenomeBiology/.")
system("cp bench_r_all.ps ~/Proyectos/ADaCGH-paper/GenomeBiology/bench_r_all.eps")




#### Plot for the ERCIM (Geneva) talk

pdf(file = "talk_bench_r-code_medium.pdf", height = 9, width = 17,
    onefile = FALSE, paper = "special")
ylim <- c(10, 12000)
xlim <- c(10 , 165)
par(mfrow = c(1, 4))
par(oma = c(0.3, 6.1, 2 , 4.8))
par(las = 1)
par(cex = 1)
par(mar = c(5, 0, 4, 0))

par(cex.lab = 2)
par(cex.main = 2.2)
par(cex.axis = 1.6)

plot(px(mediumTimingNone, "HMM"),
     ylim = ylim, lwd = 3.0, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "User wall time (seconds)",
     main = "HMM",
     xlim = xlim, axes = FALSE)
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
axis(1, at = c(20),
     labels = TRUE, cex.axis = 0.7)
axis(2, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000, 10000))
box()
points(px(mediumTiming60, "HMM"), type = "b", col = "blue",  lwd = 3.0)
points(px(mediumTiming30, "HMM"), type = "b", col = "orange",  lwd = 3.0)
points(px(mediumTiming10, "HMM"), type = "b", col = "red",  lwd = 3.0)
text(60, 1400, "Sequential code", cex = 1.8)
text(80, 150, "Parallelized code", cex = 1.8)
text(cbind(4, 0) + px(mediumTiming60, "HMM")[5, ], "60", col = "blue", adj = 0, cex = 1.9)
text(cbind(4, 0) + px(mediumTiming30, "HMM")[5, ], "30", col = "orange", adj = 0, cex = 1.9)
text(cbind(4, 0) + px(mediumTiming10, "HMM")[5, ], "10", col = "red", adj = 0, cex = 1.9)



par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "GLAD"),
     ylim = ylim, lwd = 3.0, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "GLAD",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
points(px(mediumTiming60, "GLAD"), type = "b", col = "blue",  lwd = 3.0)
points(px(mediumTiming30, "GLAD"), type = "b", col = "orange",  lwd = 3.0)
points(px(mediumTiming10, "GLAD"), type = "b", col = "red",  lwd = 3.0)


par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "CBS"),
     ylim = ylim, lwd = 3.0, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "CBS",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150),
     labels = TRUE)
points(px(mediumTiming60, "CBS"), type = "b", col = "blue",  lwd = 3.0)
points(px(mediumTiming30, "CBS"), type = "b", col = "orange",  lwd = 3.0)
points(px(mediumTiming10, "CBS"), type = "b", col = "red",  lwd = 3.0)


par(mar = c(5, 0, 4, 0))
plot(px(mediumTimingNone, "BioHMM"),
     ylim = ylim, lwd = 3.0, type = "b",
     log = "y", xaxt = "n",
     xlab = "",
     ylab = "",
     main = "BioHMM",
     xlim = xlim, axes = FALSE)
box()
axis(1, at = c(10, 20, 50, 100, 150), 
     labels = TRUE)
axis(4, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000, 10000))
points(px(mediumTiming60, "BioHMM"), type = "b", col = "blue",  lwd = 3.0)
points(px(mediumTiming30, "BioHMM"), type = "b", col = "orange",  lwd = 3.0)
points(px(mediumTiming10, "BioHMM"), type = "b", col = "red",  lwd = 3.0)

mtext("20,000 genes", side = 3, outer = TRUE, line = 0.2, cex = 2.5)
mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 2)
par(las = 0)
mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 4.5, cex = 2)
dev.off()









## load("smallTimingNoSeq30.RData")
## load("smallTimingNoSeq10.RData")
## load("mediumTimingNoSeq30.RData")
## load("mediumTimingNoSeq10.RData")


## postscript(file = "bench_r-codeNS_small.eps", height = 9, width = 17,
##     onefile = FALSE, paper = "special")
## ylim <- c(4, 450)
## xlim <- c(10 , 165)
## par(mfrow = c(1, 3))
## par(oma = c(0, 4.2, 2 , 3))
## par(las = 1)
## par(cex = 1)
## par(mar = c(5, 0, 4, 0))
## plot(px(smallTimingNoSeq10, "Wave"),
##      ylim = ylim, lwd = 1.5, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "User wall time (seconds)",
##      main = "Wavelet-based smoothing",
##      xlim = xlim, axes = FALSE,
##      col = "red")
## axis(1, at = c(10, 20, 50, 100, 150),
##      labels = TRUE, cex.axis = 0.9)
## axis(1, at = c(20),
##      labels = TRUE, cex.axis = 0.7)
## axis(2, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000))
## box()
## points(px(smallTimingNoSeq30, "Wave"), type = "b", col = "orange",  lwd = 1.5)
## text(cbind(4, 0) + px(smallTimingNoSeq30, "Wave")[5, ], "30", col = "orange", adj = 0)
## text(cbind(4, 0) + px(smallTimingNoSeq10, "Wave")[5, ], "10", col = "red", adj = 0)

## plot(px(smallTimingNoSeq10, "ACE"),
##      ylim = ylim, lwd = 1.5, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "User wall time (seconds)",
##      main = "Wavelet-based smoothing",
##      xlim = xlim, axes = FALSE,
##      col = "red")
## axis(1, at = c(10, 20, 50, 100, 150),
##      labels = TRUE, cex.axis = 0.9)
## axis(1, at = c(20),
##      labels = TRUE, cex.axis = 0.7)
## axis(2, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000))
## box()
## points(px(smallTimingNoSeq30, "ACE"), type = "b", col = "orange",  lwd = 1.5)

## plot(px(smallTimingNoSeq10, "CGHseg"),
##      ylim = ylim, lwd = 1.5, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "User wall time (seconds)",
##      main = "CGHseg",
##      xlim = xlim, axes = FALSE,
##      col = "red")
## axis(1, at = c(10, 20, 50, 100, 150),
##      labels = TRUE, cex.axis = 0.9)
## axis(1, at = c(20),
##      labels = TRUE, cex.axis = 0.7)
## axis(2, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000))
## box()
## points(px(smallTimingNoSeq30, "CGHseg"), type = "b", col = "orange",  lwd = 1.5)








## par(mar = c(5, 0, 4, 0))
## plot(px(smallTimingNone, "GLAD"),
##      ylim = ylim, lwd = 1.5, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "",
##      main = "GLAD",
##      xlim = xlim, axes = FALSE)
## box()
## axis(1, at = c(10, 20, 50, 100, 150), 
##      labels = TRUE, cex.axis = 0.9)
## points(px(smallTiming60, "GLAD"), type = "b", col = "blue",  lwd = 1.5)
## points(px(smallTiming30, "GLAD"), type = "b", col = "orange",  lwd = 1.5)
## points(px(smallTiming10, "GLAD"), type = "b", col = "red",  lwd = 1.5)



## par(mar = c(5, 0, 4, 0))
## plot(px(smallTimingNone, "CBS"),
##      ylim = ylim, lwd = 1.5, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "",
##      main = "CBS",
##      xlim = xlim, axes = FALSE)
## box()
## axis(1, at = c(10, 20, 50, 100, 150),
##      labels = TRUE, cex.axis = 0.9)
## points(px(smallTiming60, "CBS"), type = "b", col = "blue",  lwd = 1.5)
## points(px(smallTiming30, "CBS"), type = "b", col = "orange",  lwd = 1.5)
## points(px(smallTiming10, "CBS"), type = "b", col = "red",  lwd = 1.5)


## par(mar = c(5, 0, 4, 0))
## plot(px(smallTimingNone, "BioHMM"),
##      ylim = ylim, lwd = 1.5, type = "b",
##      log = "y", xaxt = "n",
##      xlab = "",
##      ylab = "",
##      main = "BioHMM",
##      xlim = xlim, axes = FALSE)
## box()
## axis(1, at = c(10, 20, 50, 100, 150), 
##      labels = TRUE, cex.axis = 0.9)
## axis(4, at = c(5, 20, 50, 100, 300, 500, 1000, 2500, 5000))
## points(px(smallTiming60, "BioHMM"), type = "b", col = "blue",  lwd = 1.5)
## points(px(smallTiming30, "BioHMM"), type = "b", col = "orange",  lwd = 1.5)
## points(px(smallTiming10, "BioHMM"), type = "b", col = "red",  lwd = 1.5)

## mtext("a) 10,000 genes", side = 3, outer = TRUE, line = 0.2, cex = 2, adj = 0)
## mtext("Number of arrays (samples)", side = 1, outer = TRUE, line = -2, cex = 1.5)
## par(las = 0)
## mtext("User wall time (seconds)", side = 2, outer = TRUE, line = 3, cex = 1.5)
## dev.off()








