#!/usr/bin/env python

import os
import sys
import time

NUM_USERS = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             5, 5,
             10,
             20,
             50)


TESTS = (
    'CBS_small',
    'HMM_small',
    'GLAD_small',
    'Wavelets_small',
    'CBS_large',
    'HMM_large',
    'GLAD_large',
    'Wavelets_large')

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
    


SUFFIX = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 1, 1, 1)
kwd = zip(NUM_USERS, SUFFIX)

TESTS = (
    'CBS_small',
    'HMM_small',
    'GLAD_small',
    'Wavelets_small')

for test in TESTS:
    for ks in kwd:
        nu = ks[0]
        timings = -99
        try:
            timings = launchUTests(test, nu)
        except:
            None
        writeFile(timings, 'web.bnchmk3.' + test + '_' + \
                  str(nu) + '.' + str(ks[1]) + '.txt')



TESTS = (
    'CBS_large',
    'HMM_large',
    'GLAD_large',
    'Wavelets_large')

for test in TESTS:
    for ks in kwd:
        nu = ks[0]
        timings = -99
        try:
            timings = launchUTests(test, nu)
        except:
            None
        writeFile(timings, 'web.bnchmk3.' + test + '_' + \
                  str(nu) + '.' + str(ks[1]) + '.txt')



test = 'Wavelets_small'
for ks in kwd:
    nu = ks[0]
    timings = -99
    try:
        timings = launchUTests(test, nu)
    except:
        None
    writeFile(timings, 'web.bnchmk3.' + test + '_' + \
              str(nu) + '.' + str(ks[1]) + '.txt')





### Next lines are a hack, because the connection was going down, and I had to restart
## the test

# timings = -99
# try:
#     timings = launchUTests('CBS_medium', 5)
# except:
#     None
# writeFile(timings, 'web.bnchmk2.CBS_medium_5.1.txt')

# timings = -99
# try:
#     timings = launchUTests('HMM_medium', 5)
# except:
#     None
# writeFile(timings, 'web.bnchmk2.HMM_medium_5.1.txt')
# timings = -99
# try:
#     timings = launchUTests('GLAD_medium', 5)
# except:
#     None
# writeFile(timings, 'web.bnchmk2.GLAD_medium_5.1.txt')

# timings = -99
# try:
#     timings = launchUTests('CGHseg_medium', 5)
# except:
#     None
# writeFile(timings, 'web.bnchmk2.CGHseg_medium_5.1.txt')


# TESTS = (
# ##    'CBS_medium',
#     'HMM_medium',
#     'GLAD_medium',
#     'CGHseg_medium')

# NUM_USERS = (1, 1, 1, 1, 1,
#              2, 2, 2)
# SUFFIX = (1, 2, 3, 4, 5, 1, 2, 3)
# kwd = zip(NUM_USERS, SUFFIX)

# for test in TESTS:
#     for ks in kwd:
#         nu = ks[0]
#         timings = -99
#         try:
#             timings = launchUTests(test, nu)
#         except:
#             None
#         writeFile(timings, 'web.bnchmk2.' + test + '_' + \
#                   str(nu) + '.' + str(ks[1]) + '.txt')

# NUM_USERS = (1, 1, 1, 1, 1,
#              2, 2, 2, 5)
# SUFFIX = (1, 2, 3, 4, 5, 1, 2, 3, 1)
# kwd = zip(NUM_USERS, SUFFIX)

# for ks in kwd:
#     nu = ks[0]
#     timings = -99
#     try:
#         timings = launchUTests('CGHseg_small', nu)
#     except:
#         None
#     writeFile(timings, 'web.bnchmk2.CGHseg_small_' + \
#               str(nu) + '.' + str(ks[1]) + '.txt')


    
