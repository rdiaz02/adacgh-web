#!/bin/sh

VERSION=$(grep Version ./ADaCGH/DESCRIPTION | sed 's/Version: //')
SRCDIR=$(pwd)

rm -r -f ADaCGH.Rcheck
rm -r -f ADaCGH.Rcheck.all.run

rm ./ADaCGH/src/*.so
rm ./ADaCGH/src/*.o

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

echo "\n\nAny errors in full run?"
grep ERROR ./ADaCGH.Rcheck.all.run/00check.log
echo "\n\nAny errors in run for CRAN?"
grep ERROR ./ADaCGH.Rcheck/00check.log
