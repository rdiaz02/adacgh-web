##########################################################################
##########################################################################

###########   The code below provides reproducible examples in a plain
###########   vanilla R using snapCGH (i.e., no need for ADaCGH).

###########   crashes with snapCGH 1.4.0 and 1.5.0 (details below)

load("dellinger2.RData")

## note: cghdata are just median-centered cgh data.
library(snapCGH)

n <- dim(cghdata)[1]
## this will crash 
out <- fit.model(sample = 1, chrom = 1, dat = cghdata[, 1, drop = FALSE],
                 datainfo = data.frame(Name = 1:n, Chrom = rep(1, n),
                 Position = positions))



## try all of them
out1 <- fit.model(sample = 1, chrom = 1, dat = cghdata,
                 datainfo = data.frame(Name = 1:n, Chrom = rep(1, n),
                 Position = positions))
out2 <- fit.model(sample = 2, chrom = 1, dat = cghdata,
                 datainfo = data.frame(Name = 1:n, Chrom = rep(1, n),
                 Position = positions))
out3 <- fit.model(sample = 3, chrom = 1, dat = cghdata,
                 datainfo = data.frame(Name = 1:n, Chrom = rep(1, n),
                 Position = positions))
out4 <- fit.model(sample = 4, chrom = 1, dat = cghdata,
                 datainfo = data.frame(Name = 1:n, Chrom = rep(1, n),
                 Position = positions))




## > version
##                _                           
## platform       x86_64-pc-linux-gnu         
## arch           x86_64                      
## os             linux-gnu                   
## system         x86_64, linux-gnu           
## status                                     
## major          2                           
## minor          5.1                         
## year           2007                        
## month          06                          
## day            27                          
## svn rev        42083                       
## language       R                           
## version.string R version 2.5.1 (2007-06-27)






## Package:       snapCGH
## Title:         Segmentation, normalisation and processing of aCGH data.
## Version:       1.4.0
## Date:          2007-03-16
## Author:        Mike L. Smith, John C. Marioni, Steven McKinney, Tom J.
##                Hardcastle, Natalie P. Thorne
## Description:   Methods for segmenting, normalising and processing aCGH
##                data; including plotting functions for visualising raw
##                and segmented data for individual and multiple arrays.
## biocViews:     Microarray, DNACopyNumber, TwoChannel, Preprocessing
## Maintainer:    John Marioni <jcm68@cam.ac.uk>
## Depends:       limma, tilingArray, DNAcopy, GLAD, cluster, methods,
##                aCGH
## License:       GPL
## Packaged:      Wed Apr 25 02:30:47 2007; biocbuild
## Built:         R 2.5.1; x86_64-pc-linux-gnu; 2007-09-05 12:20:45; unix





## 		Information on package 'snapCGH'

## Description:

## Package:       snapCGH
## Title:         Segmentation, normalisation and processing of aCGH data.
## Version:       1.5.0
## Date:          2007-03-16
## Author:        Mike L. Smith, John C. Marioni, Steven McKinney, Tom J.
##                Hardcastle, Natalie P. Thorne
## Description:   Methods for segmenting, normalising and processing aCGH
##                data; including plotting functions for visualising raw
##                and segmented data for individual and multiple arrays.
## biocViews:     Microarray, DNACopyNumber, TwoChannel, Preprocessing
## Maintainer:    John Marioni <jcm68@cam.ac.uk>
## Depends:       limma, tilingArray, DNAcopy, GLAD, cluster, methods,
##                aCGH
## License:       GPL
## Packaged:      Tue Apr 24 14:37:13 2007; biocbuild
## Built:         R 2.5.1; x86_64-pc-linux-gnu; 2007-09-05 13:10:17; unix
