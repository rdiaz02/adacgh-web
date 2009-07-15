#!/bin/sh
#cd /http/R-custom/bin
## cd /home/sources.programs/R-patched/bin
#cd /home/sources.programs/R-devel/bin


### Instal Rmpi seàrately, specifyin lam directory

alias RR="/Part-ramon/sources.programs/R-tests/R-patched-2009-07-09/bin/R"


RR CMD INSTALL /http/adacgh-required-packs/Rmpi_* --configure-args=--with-mpi=/usr/lib/lam

RR CMD INSTALL /http/adacgh-required-packs/aws*.gz
##RR CMD INSTALL /http/adacgh-required-packs/Rmpi*.gz
RR CMD INSTALL /http/adacgh-required-packs/papply*.gz
RR CMD INSTALL /http/adacgh-required-packs/polynom*.gz 
RR CMD INSTALL /http/adacgh-required-packs/sma_*.gz 
RR CMD INSTALL /http/adacgh-required-packs/GDD*.gz
RR CMD INSTALL /http/adacgh-required-packs/gridBase*.gz
##RR CMD INSTALL /http/adacgh-required-packs/convert*.gz   
RR CMD INSTALL /http/adacgh-required-packs/colorspace*.gz 
RR CMD INSTALL /http/adacgh-required-packs/hexbin*.gz 
RR CMD INSTALL /http/adacgh-required-packs/RColorBrewer*.gz 
RR CMD INSTALL /http/adacgh-required-packs/rlecuyer_*.gz
RR CMD INSTALL /http/adacgh-required-packs/DBI_*.gz
RR CMD INSTALL /http/adacgh-required-packs/RSQLite_*.gz
RR CMD INSTALL /http/adacgh-required-packs/xtable_*.gz
RR CMD INSTALL /http/adacgh-required-packs/waveslim_*.gz
RR CMD INSTALL /http/adacgh-required-packs/zoo_*.gz
RR CMD INSTALL /http/adacgh-required-packs/sandwich_*.gz
RR CMD INSTALL /http/adacgh-required-packs/strucchange_*.gz
RR CMD INSTALL /http/adacgh-required-packs/pixmap_*.gz
RR CMD INSTALL /http/adacgh-required-packs/Hmisc_*.gz
RR CMD INSTALL /http/adacgh-required-packs/combinat*.gz
RR CMD INSTALL /http/adacgh-required-packs/snow_*.gz
RR CMD INSTALL /http/adacgh-required-packs/R2HTML*.gz
RR CMD INSTALL /http/adacgh-required-packs/rsprng_*.gz
RR CMD INSTALL /http/adacgh-required-packs/cgh_*.gz 


RR CMD INSTALL /http/adacgh-required-packs/IRanges*.gz
RR CMD INSTALL /http/adacgh-required-packs/Biobase*.gz
RR CMD INSTALL /http/adacgh-required-packs/DynDoc*.gz 
RR CMD INSTALL /http/adacgh-required-packs/widgetTools*.gz
RR CMD INSTALL /http/adacgh-required-packs/tkWidgets*.gz  
RR CMD INSTALL /http/adacgh-required-packs/Biostrings_*.gz
RR CMD INSTALL /http/adacgh-required-packs/AnnotationDbi_* 
RR CMD INSTALL /http/adacgh-required-packs/preprocessCore*.gz
RR CMD INSTALL /http/adacgh-required-packs/affyio_*.gz
RR CMD INSTALL /http/adacgh-required-packs/affy_*.gz
RR CMD INSTALL /http/adacgh-required-packs/matchprobes_*.gz
RR CMD INSTALL /http/adacgh-required-packs/gcrma_*.gz
RR CMD INSTALL /http/adacgh-required-packs/affydata_*.gz
RR CMD INSTALL /http/adacgh-required-packs/affyPLM*.gz
RR CMD INSTALL /http/adacgh-required-packs/annotate_*
RR CMD INSTALL /http/adacgh-required-packs/genefilter_*.gz
RR CMD INSTALL /http/adacgh-required-packs/limma_*.gz
RR CMD INSTALL /http/adacgh-required-packs/multtest*.gz
RR CMD INSTALL /http/adacgh-required-packs/marray_*.gz
RR CMD INSTALL /http/adacgh-required-packs/vsn_*.gz
RR CMD INSTALL /http/adacgh-required-packs/annotate_*.gz
RR CMD INSTALL /http/adacgh-required-packs/geneplotter_*.gz
RR CMD INSTALL /http/adacgh-required-packs/tilingArray_*.gz
RR CMD INSTALL /http/adacgh-required-packs/DNAcopy*.gz 
RR CMD INSTALL /http/adacgh-required-packs/aCGH_*.gz 
RR CMD INSTALL /http/adacgh-required-packs/convert_*gz
RR CMD INSTALL /http/adacgh-required-packs/arrayQuality_*gz
RR CMD INSTALL /http/adacgh-required-packs/cghMCR_*.gz 
RR CMD INSTALL /http/adacgh-required-packs/GLAD_*gz
RR CMD INSTALL /http/adacgh-required-packs/snapCGH_*gz

RR CMD INSTALL /http/pomelo2/imagemap_current.tar.gz

RR CMD INSTALL /http/adacgh-required-packs/mvtnorm_*
RR CMD INSTALL /http/adacgh-required-packs/modeltools_*
RR CMD INSTALL /http/adacgh-required-packs/vcd_*
RR CMD INSTALL /http/adacgh-required-packs/coin_*
RR CMD INSTALL /http/adacgh-required-packs/zoo_*
RR CMD INSTALL /http/adacgh-required-packs/sandwich_*
RR CMD INSTALL /http/adacgh-required-packs/strucchange_*
RR CMD INSTALL /http/adacgh-required-packs/party_*
RR CMD INSTALL /http/adacgh-required-packs/combinat_*
RR CMD INSTALL /http/adacgh-required-packs/mboost_*


RR CMD INSTALL /http/adacgh-required-packs/ADaCGH_*
RR CMD INSTALL /http/adacgh-required-packs/randomForest_*

RR CMD INSTALL /http/adacgh-required-packs/survival_*

RR CMD INSTALL /http/adacgh-required-packs/gtools_*
RR CMD INSTALL /http/adacgh-required-packs/gdata_*
RR CMD INSTALL /http/adacgh-required-packs/bitops_*
RR CMD INSTALL /http/adacgh-required-packs/caTools_*
RR CMD INSTALL /http/adacgh-required-packs/gplots_*

RR CMD INSTALL ~/bzr-local-repos/signs/R-packages/SignS2
RR CMD INSTALL ~/bzr-local-repos/varSelRF









## just in case a few things did not work

# sudo R CMD INSTALL /http/adacgh-required-packs/DynDoc*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/widgetTools*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/marray*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/genefilter*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/tilingArray*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/arrayQuality*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/cghMCR*.gz
# sudo R CMD INSTALL /http/adacgh-required-packs/snapCGH*.gz




### For BioC 2.4 do


