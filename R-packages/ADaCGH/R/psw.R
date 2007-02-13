    funsw0 <- function(zz) {
        browser()
        swt.run <- sw(zz, trace = FALSE)
        swt.perm <- sw.perm.test(zz, max.nIslands = NULL,
                                 nIter = nIter,
                                 seed = NULL,
                                 trace = FALSE)
        swt.rob <- sw.rob(zz, prec = prec)
        list(swt.run = swt.run, swt.perm = swt.perm, swt.rob = swt.rob)
    }

    papout0 <- lapply(swtlist, funsw0)
