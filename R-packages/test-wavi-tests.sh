#!/bin/sh

## Runs the testing code in wavi

## You should have launched the LAM/MPI environment the way you
## like it first.

## Of course, you need to install ADaCGH first!!

##alias RR=/var/www/bin/R-local-7-LAM-MPI/bin/R
alias RR="/Part-ramon/sources.programs/R-tests/R-patched-2009-07-09/bin/R"

rm -r -f /tmp/wavi-tests
mkdir /tmp/wavi-tests
cp ../../adacgh-server/f1.R /tmp/wavi-tests/.

cp -a ../../adacgh-server/test-cases/dnacopy-ok /tmp/wavi-tests/.
cp -a ../../adacgh-server/test-cases/biohmm-ok /tmp/wavi-tests/.
cp -a ../../adacgh-server/test-cases/hmm-ok /tmp/wavi-tests/.
cp -a ../../adacgh-server/test-cases/glad-ok /tmp/wavi-tests/.
cp -a ../../adacgh-server/test-cases/wavelets-ok /tmp/wavi-tests/.
cp -a ../../adacgh-server/test-cases/cghseg-ok /tmp/wavi-tests/.
cp -a ../../adacgh-server/test-cases/haarseg-ok /tmp/wavi-tests/.

sed -i 's/.__ADaCGH_SERVER_APPL", TRUE)/.__ADaCGH_SERVER_APPL", FALSE)/' /tmp/wavi-tests/f1.R 

cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/dnacopy-ok/.
cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/biohmm-ok/.
cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/hmm-ok/.
cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/glad-ok/.
cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/wavelets-ok/.
cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/cghseg-ok/.
cp /tmp/wavi-tests/f1.R /tmp/wavi-tests/haarseg-ok/.


cd /tmp/wavi-tests
cd dnacopy-ok
RR CMD BATCH f1.R
cd ..
cd biohmm-ok
RR CMD BATCH f1.R
cd ..
cd hmm-ok
RR CMD BATCH f1.R
cd ..
cd glad-ok
RR CMD BATCH f1.R
cd ..
cd wavelets-ok
RR CMD BATCH f1.R
cd ..
cd cghseg-ok
RR CMD BATCH f1.R
cd ..
cd haarseg-ok
RR CMD BATCH f1.R
cd ..


cat ./*/R_Status.txt