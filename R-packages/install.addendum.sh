#!/bin/bash

echo "******************"
echo "start in "
hostname

/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/bitops_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/caTools_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/gtools_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/gdata_*
/var/www/bin/R-local-7-LAM-MPI/bin/R CMD INSTALL /http/adacgh-required-packs/gplots_*

echo "end in "
hostname
echo "+++++++++++++++"
