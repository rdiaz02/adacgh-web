#!/bin/bash

## since BioC dependencies are such a pain, we first download all the required
## (and maybe some extra stuff, but I no longer now).

wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/Rmpi*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/papply*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/sma_*.gz 
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/GDD*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/gridBase*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/colorspace*.gz 
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/RColorBrewer*.gz 
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/rlecuyer_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/zoo_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/sandwich_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/strucchange_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/cgh_*.gz 
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/Hmisc_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/rsprng_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/combinat*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/snow_*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/R2HTML*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/aws*.gz
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/pixmap_*.gz .
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/waveslim_*.gz .
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/xtable*.gz .
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/RSQLite_* .
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/DBI_* .
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/polynom* .
wget ftp://ftp.stat.math.ethz.ch/Software/CRAN/src/contrib/Matrix_* .


rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/IRanges*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/DynDoc*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/widgetTools*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/convert*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/hexbin*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/Biobase*.gz . 
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/Biostrings*.gz . 
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/affyio_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/affy_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/matchprobes_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/gcrma_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/affydata_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/affyPLM*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/genefilter_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/limma_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/multtest*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/marray_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/vsn_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/annotate_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/geneplotter_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/tkWidgets*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/tilingArray_*.gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/DNAcopy*.gz . 
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/aCGH_*.gz . 
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/arrayQuality_*gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/cghMCR_*.gz . 
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/GLAD_*gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/snapCGH_*gz .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/AnnotationDbi_* .
rsync -rtlv --delete bioconductor.org::2.3/data/experiment/src/contrib/affydata_* .
rsync -rtlv --delete bioconductor.org::2.3/bioc/src/contrib/preprocessCore_* .



## I am tired. I get affydata by hand, copied from another dir





