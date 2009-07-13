#!/bin/bash
#cd /http/R-custom/bin
## cd /home/sources.programs/R-patched/bin
#cd /home/sources.programs/R-devel/bin


### Instal Rmpi seàrately, specifyin lam directory

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/Rmpi_* --configure-args=--with-mpi=/usr/lib/lam

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/aws*.gz
##/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/Rmpi*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/papply*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/polynom*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/sma_*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/GDD*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/gridBase*.gz
##/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/convert*.gz   
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/colorspace*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/hexbin*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/RColorBrewer*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/rlecuyer_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/DBI_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/RSQLite_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/xtable_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/waveslim_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/zoo_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/sandwich_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/strucchange_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/pixmap_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/Hmisc_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/combinat*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/snow_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/R2HTML*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/rsprng_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/cgh_*.gz 


/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/IRanges*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/Biobase*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/DynDoc*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/widgetTools*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/tkWidgets*.gz  
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/Biostrings_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/AnnotationDbi_* 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/preprocessCore*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/affyio_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/affy_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/matchprobes_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/gcrma_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/affydata_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/affyPLM*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/annotate_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/genefilter_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/limma_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/multtest*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/marray_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/vsn_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/annotate_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/geneplotter_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/tilingArray_*.gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/DNAcopy*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/aCGH_*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/convert_*gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/arrayQuality_*gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/cghMCR_*.gz 
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/GLAD_*gz
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/snapCGH_*gz

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/pomelo2/imagemap_current.tar.gz

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/mvtnorm_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/modeltools_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/vcd_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/coin_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/zoo_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/sandwich_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/strucchange_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/party_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/combinat_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/mboost_*

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/signs2/R-packages/SignS2

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/ADaCGH_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/randomForest_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/varSelRF










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


