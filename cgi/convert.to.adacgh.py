#!/usr/bin/python

## when ready, turn adacgh2 into adacgh
## or viceversa


import sys
import os

direction = sys.argv[1]

if direction=='adacgh2':
    os.system("sed 's/adacgh.bioinfo/adacgh2.bioinfo/g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/adacgh.bioinfo/adacgh2.bioinfo/g' adacghR.cgi > tmp; mv tmp adacghR.cgi")
    os.system("sed 's/adacgh.bioinfo/adacgh2.bioinfo/g' ace.cgi > tmp; mv tmp ace.cgi")
    os.system("sed 's/adacgh.bioinfo/adacgh2.bioinfo/g' ../www/adacgh.html > tmp; mv tmp ../www/adacgh.html")
    os.system("sed 's/adacgh.bioinfo/adacgh2.bioinfo/g' results-pre.html > tmp; mv tmp results-pre.html")
    os.system("sed 's/adacgh.bioinfo/adacgh2.bioinfo/g' results-pre-ace.html > tmp; mv tmp results-pre-ace.html")    
    os.system("sed 's/\/http\/adacgh\//\/http\/adacgh2\//g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/\/http\/adacgh\//\/http\/adacgh2\//g' adacghR.cgi > tmp; mv tmp adacghR.cgi")
    os.system("sed 's/\/http\/adacgh\//\/http\/adacgh2\//g' ../www/adacgh.html > tmp; mv tmp ../www/adacgh.html")
    os.system("sed 's/\/http\/adacgh\//\/http\/adacgh2\//g' ace.cgi > tmp; mv tmp ace.cgi")
    os.system("sed 's/\/http\/adacgh\//\/http\/adacgh2\//g' results-pre.html > tmp; mv tmp results-pre.html")
    os.system("sed 's/\/http\/adacgh\//\/http\/adacgh2\//g' results-pre-ace.html > tmp; mv tmp results-pre-ace.html")    
    
    
if direction=='adacgh':
    os.system("sed 's/adacgh2.bioinfo/adacgh.bioinfo/g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/adacgh2.bioinfo/adacgh.bioinfo/g' adacghR.cgi > tmp; mv tmp adacghR.cgi")
    os.system("sed 's/adacgh2.bioinfo/adacgh.bioinfo/g' ace.cgi > tmp; mv tmp ace.cgi")
    os.system("sed 's/adacgh2.bioinfo/adacgh.bioinfo/g' ../www/adacgh.html > tmp; mv tmp ../www/adacgh.html")
    os.system("sed 's/adacgh2.bioinfo/adacgh.bioinfo/g' results-pre.html > tmp; mv tmp results-pre.html")
    os.system("sed 's/adacgh2.bioinfo/adacgh.bioinfo/g' results-pre-ace.html > tmp; mv tmp results-pre-ace.html")    
    os.system("sed 's/\/http\/adacgh2\//\/http\/adacgh\//g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/\/http\/adacgh2\//\/http\/adacgh\//g' adacghR.cgi > tmp; mv tmp adacghR.cgi")
    os.system("sed 's/\/http\/adacgh2\//\/http\/adacgh\//g' ../www/adacgh.html > tmp; mv tmp ../www/adacgh.html")
    os.system("sed 's/\/http\/adacgh2\//\/http\/adacgh\//g' ace.cgi > tmp; mv tmp ace.cgi")
    os.system("sed 's/\/http\/adacgh2\//\/http\/adacgh\//g' results-pre.html > tmp; mv tmp results-pre.html")
    os.system("sed 's/\/http\/adacgh2\//\/http\/adacgh\//g' results-pre-ace.html > tmp; mv tmp results-pre-ace.html")    

os.system('chmod u+x /asterias-web-apps/adacgh2/cgi/*.cgi')    
os.system('chown -R www-data /asterias-web-apps/adacgh2')
