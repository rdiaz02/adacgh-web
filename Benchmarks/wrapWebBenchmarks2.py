#!/usr/bin/env python

import os
import sys
import time

NUM_USERS = (1, 1, 1, 1, 1,
             2, 2, 2,
             5,
             10,
             25)

TESTS = (
    'CBS_small',
    'HMM_small',
##    'BioHMM_small',
##    'CGHseg_small',
    'GLAD_small',
    'Wavelets_small',
##    'ACE_small',
##    'PSW_small',
    'CBS_large',
    'HMM_large',
##    'BioHMM_large',
##    'CGHseg_large',
    'GLAD_large',
    'Wavelets_large')
##    'ACE_large',
##    'PSW_large')
         

def launchUTests(test, users):
    t = [-99999 for i in range(users)]
    timef = [-99999 for i in range(users)]
    for uu in range(users):
        iin, t[uu] = os.popen2('fl-run-test benchmarkADaCGH2.py ADaCGH.test' + test)
        time.sleep(15)
    for uu in range(users):
        timef[uu] = float(t[uu].readlines()[1].strip())

    return timef

def launchAll(test, lusers):
    outall = []
    for users in lusers:
        print 'test :' + test + '.   users: ' + str(users)
        outall = outall + launchUTests(test, users)
    return outall



def writeFile(testout, name):
    fout = open(name, mode = 'w')
    for result in testout:
        fout.write(str(result))
        fout.write('\t')
    fout.flush()
    fout.close()
    


SUFFIX = (1, 2, 3, 4, 5, 1, 2, 3, 1, 1, 1)

kwd = zip(NUM_USERS, SUFFIX)


for test in TESTS:
    for ks in kwd:
        nu = ks[0]
        timings = launchUTests(test, nu)
        writeFile(timings, 'web.bnchmk2.' + test + '_' + \
                  str(nu) + '.' + str(ks[1]) + '.txt')

