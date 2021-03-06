#!/bin/sh

## change the R we are using, if needed
RDIR="/Part-ramon/sources.programs/R-devel-2010-04-05"
#RDIR="/Part-ramon/sources.programs/R-tests/R-patched"
alias R=$RDIR/bin/R


VERSION=$(grep Version ./ADaCGH/DESCRIPTION | sed 's/Version: //')
SRCDIR=$(pwd)

rm -r -f ADaCGH.Rcheck
rm -r -f ADaCGH.Rcheck.all.run

rm ./ADaCGH/src/*.so
rm ./ADaCGH/src/*.o


## we need to install the package, so that mpiInit works
## when multiple node tests
R CMD build ADaCGH
R CMD INSTALL ADaCGH_$VERSION.tar.gz


## Run all, including the dontrun, as it should work in my machine
cd ./ADaCGH/man
sed -i 's/\\dontrun{/{%\\dontrun/' *
cd ../../


echo " "
echo " "
echo "*********** Boot a lam universe with 4 nodes "
lamboot -b lamb-host.def
echo " "
echo " "
echo "*********** STARTING CHECK WITHOUT DONTRUN "
echo " "
echo " "

# Beware: the function "cleanEx" that is run during the R CMD check 
# breaks the Rmpi part: require(Rmpi) no longer works, as
# library has been detached and asks for 
# Error in f(libname, pkgname) : 
#   Probably Rmpi has been detached. Please quit R.


R CMD check ADaCGH
## in case we want to check
mv ADaCGH.Rcheck ADaCGH.Rcheck.all.run

echo " "
echo " "
echo "  Halt lamuniverse"
lamhalt
lamwipe

## delete the directory of library, so we can
## emulate CRAN's behavior properly

rm -r -f $RDIR/library/ADaCGH


echo " "
echo " "
echo "*********** STARTING CHECK FOR CRAN "
echo " "
echo " "

## Now, standard check
cd ./ADaCGH/man
sed -i 's/{%\\dontrun/\\dontrun{/' *
cd ../../
R CMD build ADaCGH
R CMD check ADaCGH_$VERSION.tar.gz

echo " "
echo " "
echo "Any errors in full run?"
grep ERROR ./ADaCGH.Rcheck.all.run/00check.log

echo " "
echo " "
echo "Any errors in run for CRAN?"
grep ERROR ./ADaCGH.Rcheck/00check.log
echo " "
echo " "
