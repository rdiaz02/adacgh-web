#!/usr/bin/env python

import os
import sys
import time



def launchUTests(test, users):
    t = [-99999 for i in range(users)]
    timef = [-99999 for i in range(users)]
    for uu in range(users):
        iin, t[uu] = os.popen2('fl-run-test benchmarkADaCGH2.py ADaCGH.test' + test)
        time.sleep(5)
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




###########Sizes:
#         small: 2270 x 10 arrays
#         medium: 15000 genes x 40 subjects
#         large:  42325 genes x 40 subjects

#### Medium and small

TESTS = (
   'CGHseg_medium',
   'Wavelets_medium',
   'CBS_medium',
   'HMM_medium',
   'GLAD_medium',
   'CGHseg_small',
   'Wavelets_small',
   'CBS_small',
   'HMM_small',
   'GLAD_small')

NUM_USERS = (20, 10, 1, 1, 1, 1, 1, 5) ## used for medium and small
SUFFIX = (1, 1, 1, 2, 3, 4, 5, 1)


kwd = zip(NUM_USERS, SUFFIX)

for test in TESTS:
    for ks in kwd:
        nu = ks[0]
        timings = -99
        try:
            timings = launchUTests(test, nu)
        except:
            None
        writeFile(timings, 'web.bnchmk7.' + test + '_' + \
                  str(nu) + '.' + str(ks[1]) + '.txt')



#### Large:

TESTS = (
   'CGHseg_large',
   'Wavelets_large',
   'CBS_large',
   'HMM_large',
   'GLAD_large'
)

SUFFIX = (1, 1, 1, 2, 3, 4, 5, 1)

kwd = zip(NUM_USERS, SUFFIX)

for test in TESTS:
    for ks in kwd:
        nu = ks[0]
        timings = -99
        try:
            timings = launchUTests(test, nu)
        except:
            None
        writeFile(timings, 'web.bnchmk7.' + test + '_' + \
                  str(nu) + '.' + str(ks[1]) + '.txt')




