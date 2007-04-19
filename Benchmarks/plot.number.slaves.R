### Data from timings to find out whether 60 or 120 slaves
### are faster, and get an idea of memory consumption.
### scripts are "compare.parallel.R"

vars <- expand.grid(replicate = c(1, 2),
                    num.arrays = c(20, 50, 100, 150),
                    version = c("Array", "Array by Chromosome"),
                    method = c("ACE", "HMM", "BioHMM", "CBS"),
                    num.genes = c(20000, 40000),
                    slaves = c("Slaves: 2", "Slaves: 4")
                    )
vars$slaves <- factor(vars$slaves)
#vars$num.genes <- factor(vars$num.genes)

times <- c(29, 29, 75, 75, 152, 151, 232, 231,
           17, 12, 57, 61, 131, 392, 309, 772,
           33, 33, 32, 31, 39, 33, 39, 38,
           27, 32, 55, 56, 162, 162, 243, 286,
           583, 603, 621, 627, 617, 631, 622, 581,
           114, 119, 218, 215, 403, 440, 642, 922,
           123, 87, 81, 124, 127, 90, 87, 99,
           70, 71, 104, 104, 180, 174, 426, 325,
           34, 34, 87, 87, 181, 180, 282, 282,
           28, 26, 78, 96, 313, 375, 568, 598,
           88, 84, 84, 87, 90, 84, 93, 95,
           82, 85, 139, 122, 322, 347, 510, 623,
           1311, 1266, 1220, 1272, 1331, 1296, 1322, 1277,
           257, 248, 444, 402, 986, 973, 1735, 1840,
           191, 191, 193, 192, 204, 200, 205, 205,
           105, 103, 169, 163, 414, 367, 558, 638,
           47, 46, 118, 116, 240, 238, 361, 362,
           18, 19, 52, 60, 132, 154, 971, 339,
           35, 34, 35, 33, 34, 35, 36, 35,
           34, 35, 75, 53, 160, 152, 975, 257,
           601, 592, 592, 582, 604, 604, 608, 609,
           83, 86, 159, 150, 281, 292, 544, 720,
           101, 104, 96, 97, 116, 118, 102, 100,
           71, 36, 114, 101, 199, 181, 310, 266,
           52, 53, 135, 132, 276, 273, 437, 418,
           32, 31, 104, 79, 319, 378, 463, 489,
           93, 87, 87, 84, 87, 80, 78, 83,
           82, 81, 125, 129, 292, 308, 602, 555,
           1215, 1225, 1304, 1270, 1288, 1299, 1270, 1279,
           200, 188, 313, 288, NA, NA, NA, NA,
           190, 174, 180, 176, 184, 185, 190, 190,
           107, 107, 176, 175, 556, 319, 603, 520)

datos <- cbind(vars, times)
datos.20000 <- datos[datos$num.genes == 20000, ]
datos.40000 <- datos[datos$num.genes == 40000, ]

library(lattice)

pdf("120.vs.60.20000.pdf", height = 9, width = 12)
xyplot(times ~ num.arrays|method*slaves,
       groups = version,
       key = list(text = list(c("Array (Chromosome for ACE)",
                  "Array by Chromosome"), cex = 0.8),
       col = c("blue", "magenta"),
       points = list(pch = c(1, 1), cex = 1.3)),
#       space = "top",
       ##        corner = c(0,0),
       ##        x = 0.3, y = 0.7
       scales = list(x = list(at = c(20, 50, 100, 150), cex = 0.6),
       y = list(log = TRUE, at = c(10, 50, 100, 500, 1000, 2000)),
       cex = 0.6),
       xlab = "Number of arrays",
       ylab = "Users' wall time",
       data = datos.20000,
       main = "20,000 genes",
       cex = 1.2,
       ylim = c(8, 2500),
       par.settings = list(fontsize = list(text = 30, points = 10)))
dev.off()


pdf("120.vs.60.40000.pdf", height = 9, width = 12)
xyplot(times ~ num.arrays|method*slaves,
       groups = version,
       auto.key = FALSE,
       scales = list(x = list(at = c(20, 50, 100, 150), cex = 0.6),
       y = list(log = TRUE, at = c(10, 50, 100, 500, 1000, 2000),
       cex = 0.6),
       cex = 1.1),
       xlab = "Number of arrays",
       ylab = "Users' wall time",
       data = datos.40000,
       main = "42,325 genes",
       cex = 1.2,
       ylim = c(8, 2500), 
       par.settings = list(fontsize = list(text = 30, points = 10)))
dev.off()


system("pdftk 120.vs.60.20000.pdf 120.vs.60.40000.pdf cat output 120.vs.60.pdf")
system("pdfnup --nup 1x2 120.vs.60.pdf --outfile 120.vs.60.b.pdf")
system("cp 120.vs.60.b.pdf ~/Proyectos/ADaCGH-paper/GenomeBiology/.")
