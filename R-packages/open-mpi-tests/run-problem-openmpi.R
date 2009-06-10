## I rename .RData to saved.RData

## load("segmres.RData")
## load("saved.RData")


## geneNames <- positions.merge1$name
## chrom.numeric <- positions.merge1$chrom.numeric
## save(file = "min.data.RData", geneNames, chrom.numeric, arrayNames,
##      ymin, ymax, idtype, organism, numarrays, segmres)

library(ADaCGH)
load("min.data.RData")
library(Rmpi)


#dir.to.create <- tempdir2()
#setwd(dir.to.create)
#dir.to.create
# mpiInit()

#setwd('/var/www/adacgh-pruebas')
mpiInit(universeSize = 5)

segmentPlot(segmres, geneNames = geneNames,
            chrom.numeric = chrom.numeric,
            yminmax = c(ymin, ymax),
            idtype = idtype,
            organism = organism)

mpi.close.Rslaves()
mpi.finalize()
