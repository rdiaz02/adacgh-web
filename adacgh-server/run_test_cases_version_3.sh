#!/bin/bash

#### BEWARE: if you compare output wth version 1,
#### recall the "ok" fr wavelets and cghseg use
#### MAD, whereas those for version 1
#### use mergelevels. 

rm -r -f ./runs-tmp/tmp-3/*
rm -r -f ./runs-tmp/tmp-3-older-format/*

for i in $(ls ./test-cases-3/); do cp -a ./test-cases-3/$i ./runs-tmp/tmp-3/.; done


for i in $(ls ./runs-tmp/tmp-3/);
do echo "*********"; echo; echo; echo $i; 
    /http/adacgh-server/runADaCGHserver-3.py /http/adacgh-server/runs-tmp/tmp-3/$i;
done





### Now use the older format, but test everything!!!
## the following is ugly. could improve it

cp -a ./test-cases/0-two-chroms ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/dnacopy-ok ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/biohmm-ok ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/hmm-ok ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/glad-ok ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/wavelets-ok* ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/cghseg-ok* ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/haarseg-ok ./runs-tmp/tmp-3-older-format/.
cp -a ./test-cases/140-one-master ./runs-tmp/tmp-3-older-format/140-one
cp -a ./test-cases/140-two-master ./runs-tmp/tmp-3-older-format/140-two
cp -a ./test-cases/several-large ./runs-tmp/tmp-3-older-format/several-large

cp -a ./test-cases/dnacopy-ok ./runs-tmp/tmp-3-several-identical/.
cp -a ./test-cases/biohmm-ok ./runs-tmp/tmp-3-several-identical/.
cp -a ./test-cases/hmm-ok ./runs-tmp/tmp-3-several-identical/.
cp -a ./test-cases/glad-ok ./runs-tmp/tmp-3-several-identical/.
cp -a ./test-cases/wavelets-ok* ./runs-tmp/tmp-3-several-identical/.
cp -a ./test-cases/cghseg-ok* ./runs-tmp/tmp-3-several-identical/.
cp -a ./test-cases/haarseg-ok ./runs-tmp/tmp-3-several-identical/.

for i in $(ls ./runs-tmp/tmp-3-several-identical/); do cp ./test-cases/inputData.RData.9.arrays.some.identical ./runs-tmp/tmp-3-several-identical/$i/inputData.RData; done



echo 
echo
echo
echo "Running older format"

for i in $(ls ./runs-tmp/tmp-3-older-format/);
do echo "*********"; echo; echo; echo $i; 
    /http/adacgh-server/runADaCGHserver-3.py /http/adacgh-server/runs-tmp/tmp-3-older-format/$i;
done

echo 
echo
echo
echo "Running Several identical"

for i in $(ls ./runs-tmp/tmp-3-several-identical/);
do echo "*********"; echo; echo; echo $i; 
    /http/adacgh-server/runADaCGHserver-3.py /http/adacgh-server/runs-tmp/tmp-3-several-identical/$i;
done



echo 
echo
echo
echo "Verify output"
echo  "****************"


for i in $(ls ./runs-tmp/tmp-3/); 
do cd ./runs-tmp/tmp-3/$i; 
    echo "*********"; echo; echo; echo $i; 
    cat R_Status.txt; 
    cd /http/adacgh-server; 
done


echo
echo "Verify output older format"
echo  "****************"


for i in $(ls ./runs-tmp/tmp-3-older-format/); 
do cd ./runs-tmp/tmp-3-older-format/$i; 
    echo "*********"; echo; echo; echo $i; 
    cat R_Status.txt; 
    cd /http/adacgh-server; 
done



echo
echo "Verify output several identical"
echo  "****************"


for i in $(ls ./runs-tmp/tmp-3-several-identical/); 
do cd ./runs-tmp/tmp-3-several-identical/$i; 
    echo "*********"; echo; echo; echo $i; 
    cat R_Status.txt; 
    cd /http/adacgh-server; 
done




