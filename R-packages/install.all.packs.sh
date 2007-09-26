#!/bin/bash
#cd /http/R-custom/bin
## cd /home/sources.programs/R-patched/bin
#cd /home/sources.programs/R-devel/bin


./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/aws*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/Rmpi*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/papply*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/sma*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/GDD*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/gridBase*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/convert*.gz   
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/colorspace*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/hexbin*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/RColorBrewer*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/rlecuyer_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/DBI_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/RSQLite_*.gz


./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/Biobase*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/DynDoc*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/widgetTools*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/tkWidgets*.gz  
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/Biostrings*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/affyio_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/affy_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/matchprobes_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/gcrma_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/affydata_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/affyPLM*.gz


./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/genefilter_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/limma_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/multtest*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/marray_*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/vsn_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/waveslim_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/xtable_*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/zoo_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/sandwich_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/strucchange_*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/annotate_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/geneplotter_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/pixmap_*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/tilingArray_*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/DNAcopy*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/cgh_*.gz 
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/aCGH_*.gz 

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/arrayQuality_*gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/cghMCR_*.gz 

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/GLAD_*gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/snapCGH_*gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/Hmisc_*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/rsprng_*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/combinat*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/snow*.gz

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/R2HTML*.gz

## just in case a few things did not work

./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/DynDoc*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/widgetTools*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/marray*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/genefilter*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/tilingArray*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/arrayQuality*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/cghMCR*.gz
./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/snapCGH*.gz







# #### If you need BioC-2.1 uncomment lines below and comment all above

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/aws*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/Rmpi*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/papply*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/sma*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/GDD*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/gridBase*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/convert*.gz   
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/colorspace*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/hexbin*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/RColorBrewer*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/rlecuyer_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/DBI_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/RSQLite_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/xtable_*.gz


# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/Biobase*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/DynDoc*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/widgetTools*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/tkWidgets*.gz  
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/Biostrings*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/AnnotationDbi_* 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/preprocessCore*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/affyio_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/affy_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/matchprobes_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/gcrma_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/affydata_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/affyPLM*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/annotate_*

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/genefilter_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/limma_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/multtest*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/marray_*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/vsn_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/waveslim_*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/zoo_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/sandwich_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/strucchange_*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/annotate_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/geneplotter_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/pixmap_*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/tilingArray_*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/DNAcopy*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/cgh_*.gz 
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/aCGH_*.gz 

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/arrayQuality_*gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/cghMCR_*.gz 

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/GLAD_*gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/snapCGH_*gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/Hmisc_*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/rsprng_*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/combinat*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/snow*.gz

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/R2HTML*.gz

# ## just in case a few things did not work

# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/DynDoc*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/widgetTools*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/marray*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/genefilter*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/tilingArray*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/arrayQuality*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/cghMCR*.gz
# ./R CMD INSTALL ~/bzr-local-repos/adacgh-required-packs/BioC-2.1/snapCGH*.gz


