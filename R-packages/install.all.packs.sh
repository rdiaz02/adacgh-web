#!/bin/bash
#cd /http/R-custom/bin
## cd /home/sources.programs/R-patched/bin
#cd /home/sources.programs/R-devel/bin


### Instal Rmpi seàrately, specifyin lam directory

/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/Rmpi_0.5-7.tar.gz --configure-args=--with-mpi=/usr/lib/lam

/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/aws*.gz
##/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/Rmpi*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/papply*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/polynom*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/sma_*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/GDD*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/gridBase*.gz
##/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/convert*.gz   
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/colorspace*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/hexbin*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/RColorBrewer*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/rlecuyer_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/DBI_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/RSQLite_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/xtable_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/waveslim_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/zoo_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/sandwich_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/strucchange_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/pixmap_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/Hmisc_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/combinat*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/snow_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/R2HTML*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/rsprng_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/cgh_*.gz 


/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/IRanges*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/Biobase*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/DynDoc*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/widgetTools*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/tkWidgets*.gz  
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/Biostrings_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/AnnotationDbi_* 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/preprocessCore*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/affyio_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/affy_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/matchprobes_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/gcrma_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/affydata_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/affyPLM*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/annotate_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/genefilter_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/limma_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/multtest*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/marray_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/vsn_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/annotate_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/geneplotter_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/tilingArray_*.gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/DNAcopy*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/aCGH_*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/convert_*gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/arrayQuality_*gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/cghMCR_*.gz 
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/GLAD_*gz
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/snapCGH_*gz

/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/ADaCGH_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/pomelo2/imagemap_current.tar.gz

/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/mvtnorm_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/modeltools_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/vcd_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/coin_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/zoo_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/sandwich_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/strucchange_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/party_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/combinat_*
/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/mboost_*

/Part-ramon/sources.programs/R-tests/R-patched/bin/R CMD INSTALL ~/bzr-local-repos/signs/R-packages/SignS2







## just in case a few things did not work

# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/DynDoc*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/widgetTools*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/marray*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/genefilter*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/tilingArray*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/arrayQuality*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/cghMCR*.gz
# sudo R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs-2009-06/snapCGH*.gz




### For BioC 2.4 do


