#!/usr/bin/python

## All this code is copyright Ramon Diaz-Uriarte.
## Released under the Affero GPL.


## try to use add_to_MPIErrorLog and similar mechanisms for logging
## MPI attempts

### FIXME: MAX_DURATION_TRY: set a constant, and use that. now, we set (two) different
##         values inside the code itself


### FIXME: get rid of check_tping!!!



### For R: inputData.RData
### We return to caller a directory name
### Files:
##        Status_Server_Run.msg
##        Error_msg.txt


import sys
import os
import cgi 
import time
import shutil
import glob
import random
import socket
import fcntl

# sys.path = sys.path + ['/http/mpi.log']
# import counterApplications


R_MAX_time = 196 * 3600 ## 12 hours is max duration allowd for any process
TIME_BETWEEN_CHECKS = 23
MAX_MPI_CRASHES = 20


#tmpDir = sys.argv[1]
inputDir = sys.argv[1]

#ROOT_TMP_DIR = "/http/adacgh2/www/tmp/"
#newDir = tmpDir.replace(ROOT_TMP_DIR, "")

appDir       = "/http/adacgh-server"
appProcs     = appDir + '/runs-tmp'
runningProcs = appProcs + '/R.running.procs'
tmpDD        = appProcs + '/tmp'
logsDir      = '/http/adacgh-server/runs-tmp/logs/'
Rexec        = '/var/www/bin/R-2.10'



#########################################################

####  These are copied from /http/mpi.log/counterApplications.py
####     but I put code and paths here, for better self containment


########################################################

def add_to_log(applicacion, tmpDir, hostname):
    date_time = time.strftime('%Y\t%m\t%d\t%X')    
    outstr = '%s\t%s\t%s\t%s\n' % (applicacion, date_time, hostname, tmpDir)
    cf = open(logsDir + 'ApplicationCounter', mode = 'a')
    fcntl.flock(cf.fileno(), fcntl.LOCK_SH)
    cf.write(outstr)
    fcntl.flock(cf.fileno(), fcntl.LOCK_UN)
    cf.close()

    
def add_to_MPIErrorLog(application, tmpDir, hostname, message = 'MPI crash'):
#     if not os.path.exists(logsDir + application + 'MPIErrorLog'):
#         os.system('touch ' + logsDir + application + 'MPIErrorLog')
    outlog = open(logsDir + application + 'MPIErrorLog', mode = 'a')
    fcntl.flock(outlog.fileno(), fcntl.LOCK_SH)
    outlog.write(message + time.ctime(time.time()) +
                 '   Directory: ' + tmpDir +
                 '   Hostname: ' + hostname + '\n')
    fcntl.flock(outlog.fileno(), fcntl.LOCK_UN)
    outlog.close()


## the following much better from R, since we know the pid of R
## but R could fail, and we would want to see where is the lam suffix
## comming from
    
def add_to_LAM_SUFFIX_LOG(lamSuffix, application, tmpDir, hostname,
                          Rprocess = 'RprocessPid'):
#     if not os.path.exists(logsDir + 'LAM_SUFFIX_Log'):
#         os.system('touch /http/mpi.log/LAM_SUFFIX_Log')
    outlog = open(logsDir + 'LAM_SUFFIX_Log', mode = 'a')
    fcntl.flock(outlog.fileno(), fcntl.LOCK_SH)
    outlog.write(Rprocess + '\t' + lamSuffix + '\t' +
                 hostname + '\t' + 
                 tmpDir + '\t' +
                 '\t' + application +  '\t' + 
                 str(time.time()) + '\t' +
                 time.ctime(time.time()) + '\n')
    fcntl.flock(outlog.fileno(), fcntl.LOCK_UN)
    outlog.close()

#########################################################
########################################################



def writeStatusServer(message, outfile = 'Status_Server_Run.msg'):
    ## we take tmpDir as a given! from local envir
    statusfile = open(tmpDir + '/' + outfile, mode = 'w')
    statusfile.write(message)
    statusfile.close()


def writeErrorMessage(message = '', infile = 'R_Error_msg.txt',
                      outfile = 'Error_msg.txt'):
    """ Write the error message for users by 
    concatenating the message and the infile."""
    if (message == '') and (infile != 'NULL'):
        shutil.copyfile(tmpDir + '/' + infile,
                        tmpDir + '/Error_msg.txt')
    elif (message != '') and (infile == 'NULL'):
        errorfile = open(tmpDir + '/' + outfile, mode = 'w')
        errorfile.write(message)
        errorfile.close()
    elif (message != '') and (infile != 'NULL'):
        errorfile = open(tmpDir + '/' + outfile, mode = 'w')
        infr = open(tmpDir + '/' + infile, mode = 'r').read()
        errorfile.write(message +'\n\n\n' + infr)
        errorfile.close()
    elif (message == '') and (infile == 'NULL'):
        pass

    
def create_tmpDir(baseDir = tmpDD):
    """ Create a new directory with appropriate permissions"""
    newDir = str(random.randint(1, 10000000)) + str(int(time.time())) + str(os.getpid())
    tmpDir = baseDir + "/" + newDir
    os.mkdir(tmpDir)
    os.chmod(tmpDir, 0700)
    return (newDir, tmpDir)

def set_defaults_lam(tmpDir):
    """ Set defaults for lamboot and Rslaves and number procs
    based on the size of the data file. This is all heuristic,
    but works for us with 6 GB RAM per node. The key idea is to
    prevent swapping. ncpu are the Rslaves spawned by lamd, or the cpu=ncpu
    in the lamb-host file. max_num_procs is the maximum number of simultaneous
    adacgh processes running at any time.
    We return the tuple ncpu, max_num_procs"""
#     datsize1 = 0
#     datsize2 = 0
#     if os.path.exists(tmpDir + '/acghData'):
#         datsize1 = int(os.popen('ls ' + tmpDir + '/acghData -sk').read().split()[0])
#     if os.path.exists(tmpDir + '/acghAndPosition'):
#         datsize2 = int(os.popen('ls ' + tmpDir + '/acghAndPosition -sk').read().split()[0])
    # datsize = int(os.popen('ls ' + tmpDir + '/inputData.RData -sk').read().split()[0])
    # if datsize < 2000:
    #     return (2, 4)
    # elif datsize < 6000:
    #     return (2, 3)
    # elif datsize < 14000:
    #     return (1, 3)
    # else:
    #     return (1, 2) 
    ## With ff, no longer necessary all that mess.
    return (2, 4)


def issue_echo(fecho, tmpDir):
    """Silly function to output small tracking files"""
    timeHuman = '##########   ' + \
                str(time.strftime('%d %b %Y %H:%M:%S')) 
    os.system('echo "' + timeHuman + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "' + fecho + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "    " >> ' + tmpDir + '/checkdone.echo')


def kill_pid_machine(pid, machine):
    'as it says: to kill somehting somewhere and do not get error messages'
    fi, fo, fu = os.popen3('ssh ' + machine + ' "kill -s 9 ' + pid + '"')
    fi.close()
    fo.close()
    fu.close()


def printErrorRun(errorfile):
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    errormsg = open(errorfile).read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>ADaCGH results </title></head><body>\n")
    outf.write("<h1> ERROR: There was a problem with the R code </h1> \n")
    outf.write("<p>  This could be a bug on our code, or a problem  ")
    outf.write("with your data (that we hadn't tought of). Below is all the output from the execution ")
    outf.write("of the run. Unless it is obvious to you that this is a fault of your data ")
    outf.write("(and that there is no way we could have avoided the crash) ")
    outf.write("please let us know so we can fix the problem. ")
    outf.write("Please sed us this URL and the output below</p>")
    outf.write("<p> This is the results file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("<p> And this is the error message:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(errormsg))
    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")



def printRKilled():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>ADaCGH results </title></head><body>\n")
    outf.write("<h1> ERROR: R process killed </h1> \n")
    outf.write("<p>  The R process lasted longer than the maximum  allowed time, ")
    outf.write(str(R_MAX_time))
    outf.write(" seconds,  and was killed.")
###     outf.write("<p> This is the output from the R run:<p>")
###     outf.write("<pre>")
###     outf.write(cgi.escape(soFar))
###     outf.write("</pre>")
    outf.write("<p> This is the results file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


def logMPIerror(tmpDir, numtries, application = 'ADaCGH-server'):
    if not os.path.exists(logsDir + application + 'ErrorLog'):
        os.system('touch ' + logsDir  + application + 'ErrorLog')
    outlog = open(logsDir + application + 'ErrorLog', mode = 'a')
    outlog.write('MPI fails more than ' + numtries + 'numtries on ' +
                 time.ctime(time.time()) +
                 ' Directory: ' + tmpDir + '\n')
    outlog.close()
#     out1 = open(tmpDir + "/natural.death.pid.txt", mode = "w")
#     out2 = open(tmpDir + "/kill.pid.txt", mode = "w")
#     out1.write('MPI initialization error!!')
#     out2.write('MPI initialization error!!')
#     out1.close()
#     out2.close()
#     outf = open(tmpDir + "/pre-results.html", mode = "w")
#     outf.write("<html><head><title> MPI initialization problem.</title></head><body>\n")
#     outf.write("<h1> MPI initialization problem.</h1>")
#     outf.write("<p> After " + str(numtries) + " attempts we have been unable to ")
#     outf.write(" initialize MPI.</p>")
#     outf.write("<p> We will be notified of this error, but we would also ")
#     outf.write("appreciate if you can let us know of any circumstances or problems ")
#     outf.write("so we can diagnose the error.</p>")
#     outf.write("</body></html>")
#     outf.close()
#     shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


def logMPITooBusy(tmpDir, MAX_DURATION_TRY, application = 'ADaCGH-server'):
    if not os.path.exists(logsDir + application + 'ErrorLog'):
        os.system('touch ' + logsDir + application + 'ErrorLog')
    outlog = open(logsDir + application + 'ErrorLog', mode = 'a')
    outlog.write('MPI too busy on ' + time.ctime(time.time()) +
                 ' Directory: ' + tmpDir + '\n')
    outlog.close()
#     out1 = open(tmpDir + "/natural.death.pid.txt", mode = "w")
#     out2 = open(tmpDir + "/kill.pid.txt", mode = "w")
#     out1.write('Cannot start!!')
#     out2.write('Cannot start!!')
#     out1.close()
#     out2.close()
#     outf = open(tmpDir + "/pre-results.html", mode = "w")
#     outf.write("<html><head><title> Cannot start application.</title></head><body>\n")
#     outf.write("<h1> Cannot start application.</h1>")
#     outf.write("<p> After " + str(MAX_DURATION_TRY) + " seconds we have been unable to ")
#     outf.write(" start the application.</p>")
#     outf.write("<p> Most likely this means the servers are too busy and many ")
#     outf.write("are running ahead of yours. ")
#     outf.write("Please try again later. You can also get in touch with us ")
#     outf.write("if you think this is our error.</p>")
#     outf.write("</body></html>")
#     outf.close()
#     shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")



def lamboot(lamSuffix, ncpu, runningProcs = runningProcs):
    'Boot a lam universe and leave a sentinel file behind'
    issue_echo('before sentinel inside lamboot', tmpDir)
    issue_echo('newDir is ' + newDir, tmpDir)
    issue_echo('lamSuffix ' + lamSuffix, tmpDir)
    issue_echo('runningProcs ' + runningProcs, tmpDir)
    sentinel = os.open(''.join([runningProcs, '/sentinel.lam.', newDir, '.', lamSuffix]),
                       os.O_RDWR | os.O_CREAT | os.O_NDELAY)
    issue_echo('before fullCommand inside lamboot', tmpDir)
    fullCommand = 'export LAM_MPI_SESSION_SUFFIX="' + lamSuffix + \
                  '"; /http/mpi.log/tryBootLAM2.py ' + lamSuffix + \
                  ' ' + str(ncpu)
    issue_echo('before os.system inside lamboot', tmpDir)
    lboot = os.system(fullCommand)
    issue_echo('after lboot ---os.system--- inside lamboot. Exiting lamboot', tmpDir)


def check_tping(lamSuffix, tmpDir, tsleep = 15, nc = 2):
    """ Use tping to verify LAM universe OK.
    tsleep is how long we wait before checking output of tping.
    Verify also using 'lamexec C hostname' """
    
    tmp2 = os.system('export LAM_MPI_SESSION_SUFFIX="' +\
                     lamSuffix + '"; cd ' + tmpDir + \
                     '; tping C N -c' + str(nc) + \
                     ' > tping.out & ')
    time.sleep(tsleep)
    tmp = int(os.popen('cd ' + tmpDir + \
                       '; wc tping.out').readline().split()[0])
    os.system('rm ' + tmpDir + '/tping.out')
    timeHuman = '##########   ' + \
                str(time.strftime('%d %b %Y %H:%M:%S')) 
    os.system('echo "' + timeHuman + \
              '" >> ' + tmpDir + '/checkTping.out')
    if tmp == 0:
        os.system('echo "tping fails" >> ' + \
                  tmpDir + '/checkTping.out')
        return 0
    elif tmp > 0:
        os.system('echo "tping OK" >> ' + \
                  tmpDir + '/checkTping.out')
        lamexec = os.system('export LAM_MPI_SESSION_SUFFIX="' +\
                            lamSuffix + '"; lamexec C hostname')
        if lamexec == 0:
            os.system('echo "lamexec OK" >> ' + \
                      tmpDir + '/checkTping.out')
            return 1
        else:
            os.system('echo "lamexec fails" >> ' + \
                      tmpDir + '/checkTping.out')
            return 0
    else:
        os.system('echo "tping weird ' + str(tmp) + '" >> ' + \
                  tmpDir + '/checkTping.out')
        return 0



def lam_crash_log(tmpDir, value):
    """ Write to the lam crash log, 'recoverFromLAMCrash.out' """
    timeHuman = str(time.strftime('%d %b %Y %H:%M:%S')) 
    os.system('echo "' + value + '  at ' + timeHuman + \
              '" >> ' + tmpDir + '/recoverFromLAMCrash.out')
    
def recover_from_lam_crash(tmpDir, NCPU, MAX_NUM_PROCS, lamSuffix,
                           runningProcs= runningProcs,
                           machine_root = 'karl'):
    """ lam crashed during R run. Restart R
    after possibly rebooting the lam universe.
    Leave a trace of what happened."""
    
    os.remove(''.join([runningProcs, '/sentinel.lam.', newDir, '.', lamSuffix]))
    del_mpi_logs(tmpDir, machine_root)
    lam_crash_log(tmpDir, 'Crashed')
    ## We need to halt the universe, or else we can keep a lamd with no R hanging from
        ## it, but that leads to too many lamds, so it cannot start. Like a vicious circle
    lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    try:
        os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                  '; lamhalt -H; lamwipe -H')
    except:
        None
    issue_echo('inside recover_from_lam_crash: lamhalting', tmpDir)
    try:
        os.system('mv ' + tmpDir + '/mpiOK ' + tmpDir + '/previous_mpiOK')
    except:
        None
    check_room = my_queue(MAX_NUM_PROCS)
    if check_room == 'Failed':
        printMPITooBusy(tmpDir, MAX_DURATION_TRY = 5 * 3600)

    ## how could the next be anything but 0? we have done a lamhalt!! FIXME
    lam_ok = check_tping(lamSuffix, tmpDir)
    if lam_ok == 0:
        lboot = lamboot(lamSuffix, NCPU)
    Rrun(tmpDir, lamSuffix)
    lam_crash_log(tmpDir, '..... recovering')


def Rrun(tmpDir, lamSuffix):
    """ Launch R, after setting the lam stuff."""
    Rcommand = 'export LAM_MPI_SESSION_SUFFIX="' + lamSuffix + \
               '"; cd ' + tmpDir + \
               '; sleep 1; ' + Rexec + \
               ' --no-readline --no-save --slave <f3.R >>f3.Rout 2>> R_Status.txt &'
    Rtorun = os.system(Rcommand)
    
### FIXME: something potentially weird here:
    ## R redirects standard error to R_Status.txt
    ## but R also writes to that very file when it exits. (caughtOurError, .Last, etc).
    ## this makes life simpler for "status_run" but ...
    ## but then, a normal execution should leave nothing in stderr.



def status_run(tmpDir):
    """ Read R_Status.txt and return status."""
    status_r = open(tmpDir + '/R_Status.txt').read()
    if status_r.find('Execution halted\n') > -1:
        return('R_ExecutionHalted')
    if status_r.find('Other Error\n') > -1:
        return('R_Other_Error')
    if status_r.find('User Error\n') > -1:
        return('R_User_Error')
    if status_r.find('Our Error\n') > -1:
        return('R_Our_Error')
    if status_r.find('Rmpi error\n') > -1:
        return('Error_Rmpi')
    if status_r.find('Running\n') > -1:
        return('Running')
    if status_r.find('Normal termination\n') > -1:
        return('R_NormalTermination')
#     if status_r.find('Run out of time; killed\n') > -1:
#         return('Out_of_time')


def did_R_crash_in_slaves(tmpDir, machine_root = 'karl'):
    """ Verify whether R crashed in any of the slaves by
    checking lam logs."""
    R_LAM_MSGS = 'Error:  Error in'
    lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
    in_lam_logs = 0
    for lam_log in lam_logs:
        # I could do it more pythonically:
        # t1 = open(lam_log).read().find(R_LAM_MSGS)
        # if t1 > -1 ...
        # But use a try to protect for non-existing files
        
        tmp1 = int(os.popen('grep "' + R_LAM_MSGS + '" ' + \
                            lam_log + ' | wc').readline().split()[0])
        if tmp1 > 0:
            in_lam_logs = 1
            break
    if in_lam_logs > 0:
        return True, lam_log
    else:
        return False, 'NA'



def did_lam_crash(tmpDir, machine_root = 'karl'):
    """ Verify whether LAM/MPI crashed by checking logs and f3.Rout
    for single universe lamboot."""
    OTHER_LAM_MSGS = 'Call stack within LAM:'
    lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
    in_error_msg = int(os.popen('grep MPI_Error_string ' + \
                                tmpDir + '/R_Status.txt | wc').readline().split()[0])
#     no_universe = int(os.popen('grep "Running serial version of papply" ' + \
#                                tmpDir + '/f1.Rout | wc').readline().split()[0])
## We do NOT want that, because sometimes a one node universe is legitimate!!!
    if in_error_msg > 0:
        for lam_log in lam_logs:
            os.system('rm ' + lam_log)
#     elif no_universe > 0:
#         os.system("sed -i 's/Running serial version of papply/already_seen:running serial version of papply/g'" + \
#                   tmpDir + "/f1.Rout")
    else: ## look in lam logs
        in_lam_logs = 0
        for lam_log in lam_logs:
            tmp1 = int(os.popen('grep "' + OTHER_LAM_MSGS + '" ' + \
                                lam_log + ' | wc').readline().split()[0])
            if tmp1 > 0:
                in_lam_logs = 1
                break
    if (in_error_msg > 0) or (in_lam_logs > 0):
        return True
    else:
        return False
    
def did_mpi_crash(tmpDir, machine_root = 'karl'):
    """ Either Rmpi or LAM crashed"""
    if (status_run(tmpDir) == 'Error_Rmpi') or \
       did_lam_crash(tmpDir, machine_root):
        return True
    else:
        return False

def del_mpi_logs(tmpDir, machine_root = 'karl'):
    """ Delete logs from LAM/MPI."""
    lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
    try:
        os.system('rm ' + tmpDir + '/R_Status.txt')
    except:
        None
    try:
        for lam_log in lam_logs:
            os.system('rm ' + lam_log)    
    except:
        None

def did_run_out_of_time(tmpDir, R_MAX_time):
    """ Did the process run longer than allowed?"""
    issue_echo('did we run out of time?', tmpDir)
    if not os.path.exists(tmpDir + "/pid.txt"):
        return False
    elif ((time.time() - os.path.getmtime(tmpDir + "/pid.txt")) > R_MAX_time) and \
       (status_run(tmpDir) == 'Running'):
        return True
    else:
        return False
                           

def cleanups(tmpDir, newDir,
             lamSuffix,
             runningProcs= runningProcs,
             newnamepid = 'finished_pid.txt'):
    """ Clean up actions; kill lam, delete running.procs files, clean process table."""
    lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    try:
        rinfo = open(tmpDir + '/current_R_proc_info', mode = 'r').readline().split()
    except:
        None
    try:
        killpidmachine = kill_pid_machine(rinfo[1], rinfo[0])
    except:
        None
    try:
        os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                  '; lamhalt -H; lamwipe -H')
    except:
        None
    try:
        os.system('rm ' + runningProcs + '/R.' + newDir + '*')
    except:
        None
    try:
        os.rename(tmpDir + '/pid.txt', tmpDir + '/' + newnamepid)
    except:
        None
    try:
        os.remove(''.join([runningProcs, '/sentinel.lam.', newDir, '.', lamSuffix]))
    except:
        None


# def finished_ok(tmpDir):
#     """ check ok termination and act accordingly."""
#     if status_run(tmpDir) == 'FinishedOK':
#         return True
#     else:
#         return False

# def halted(tmpDir):
#     """ check halted execution and act accordingly."""
#     if status_run(tmpDir) == 'Halted':
#         return True
#     else:
#         return False


def master_out_of_time(time_start):
    """If this process run longer than allowed, kill it and kill lam and R."""
    if (time.time () - time_start) > R_MAX_time:
        return True
    else:
        return False
        




def del_from_proc_table(del_procs = 1):
    """Decrease count in the process table."""
    fo = open(procTable, mode = 'r+')
    currentProcs = int(fo.read())
    fo.seek(0)
    fo.write(str(currentProcs - del_procs))
    fo.close()
    return 'OK'


def my_queue(MAX_NUM_PROCS,
             runningProcs = runningProcs,
             ADD_PROCS = 1,
             CHECK_QUEUE = 23,
             MAX_DURATION_TRY = 25 * 3600):
    """ Wait here until the number of processes is smaller than
    MAX_NUM_PROCS and number of slaves smaller than MAX_NUM_PROCS + ADD_PROCS
    (so we allow for other apps. launching lamd).
    But only wait for MAX_DURATION. Check
    every CHECK_QUEUE seconds. If able to find an opening, return
    OK, otherwise return Failed"""
    out_value = 'OK'
    startTime = time.time()
## We no longer check for room!!
    return out_value
#     while True:
# ##        killedlamandr = os.system('/http/mpi.log/killOldLamAllMachines.py')
#         issue_echo('     inside my_queue ', tmpDir)
#         if (time.time() - startTime) > MAX_DURATION_TRY:
#             out_value = 'Failed'
#             break
#         num_lamd = int(os.popen('pgrep -u www-data lamd | wc').readline().split()[0])
#         num_sentinel = int(len(glob.glob(''.join([runningProcs, '/sentinel.lam.*']))))
#         if (num_lamd < (MAX_NUM_PROCS + ADD_PROCS)) and (num_sentinel < MAX_NUM_PROCS):
#             issue_echo('     OK; num_lamd = ' + str(num_lamd) + \
#                        '; num_sentinel = ' + str(num_sentinel), tmpDir)
#             break
#         else:
# 	    issue_echo('     wait:  num_lamd = ' + str(num_lamd) + \
#                        '; num_sentinel = ' + str(num_sentinel), tmpDir)
#             time.sleep(CHECK_QUEUE + random.uniform(0.1, 5))
 

def generate_lam_suffix(tmpDir):
    """As it says. Generate and write it out"""
    lamSuffix = str(int(time.time())) + \
                str(os.getpid()) + str(random.randint(10, 999999))
    lamenvfile = open(tmpDir + '/lamSuffix', mode = 'w')
    lamenvfile.write(lamSuffix)
    lamenvfile.flush()
    lamenvfile.close()
    return lamSuffix

######################################################################
######################################################################
######################################################################

## Starting. First, the very first run.


## Clean up lam
os.system('/http/mpi.log/killOldLam.py')


##(newDir, tmpDir) = create_tmpDir()

tmpDir = inputDir
tmp1 = tmpDir.split('/')
newDir = tmp1[len(tmp1) - 1]


print tmpDir

os.system('echo "' + str(os.getpid()) + ' ' + socket.gethostname() +\
           '"> ' + tmpDir + '/run_and_checkPID')


issue_echo('starting', tmpDir)
writeStatusServer('Running')

## put data from inputDir and f1.R and prepare other initial files
# shutil.copyfile(inputDir + '/inputData.RData',
#                 tmpDir + '/inputData.RData')
# shutil.copyfile(inputDir + '/options.txt',
#                 tmpDir + '/options.txt')

shutil.copyfile(appDir + '/f3.R',
                tmpDir + '/f3.R')
os.system('/bin/touch ' + tmpDir + '/R_Error_msg.txt')
touchRout = os.system("/bin/touch " + tmpDir + "/f3.Rout") 
touchRrunning = os.system('/bin/touch ' +
                          runningProcs + '/R.' + newDir +
                          "@" + socket.gethostname())
checkpoint = os.system("echo 0 > " + tmpDir + "/checkpoint.num")

       
NCPU, MAX_NUM_PROCS = set_defaults_lam(tmpDir)

try:
    add_to_log('ADaCGH-server', tmpDir, socket.gethostname())
except:
    None

issue_echo('at 2', tmpDir)

lamSuffix = generate_lam_suffix(tmpDir)

issue_echo('at 3', tmpDir)

time.sleep(random.uniform(0.1, 3)) ## Break ties if starting at identical times

issue_echo('at 4', tmpDir)


# check_room = my_queue(MAX_NUM_PROCS)
issue_echo('after check_room', tmpDir)

# if check_room == 'Failed':
#     printMPITooBusy(tmpDir, MAX_DURATION_TRY = 5 * 3600)
#     writeStatusServer('ERROR!!!\n')
#     writeErrorMessage('MPI too busy. Too many jobs in the servers.' +
#                       '\nPlease try later', 'NULL')
#     sys.exit()

issue_echo('before lamboot', tmpDir)
lamboot(lamSuffix, NCPU)
issue_echo('after lamboot', tmpDir)

add_to_LAM_SUFFIX_LOG(lamSuffix, 'ADaCGH-server', tmpDir,
                      socket.gethostname())

## start R


Rrun(tmpDir, lamSuffix)
issue_echo('         at          A1', tmpDir)
        
time_start = time.time()
time.sleep(TIME_BETWEEN_CHECKS + random.uniform(0.1, 3))

count_mpi_crash = 0

### FIXME: maybe we can get rid of the issue_echo calls


while True:
    if (status_run(tmpDir) == 'R_NormalTermination'):
        issue_echo(status_run(tmpDir), tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('Normal termination\n')
        ### FIXME: what about PaLS, etc? See old "printOKRun()"
        break
    elif (status_run(tmpDir) == 'R_User_Error'):
        
        issue_echo(status_run(tmpDir), tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('User ERROR\n')
        writeErrorMessage('', 'R_Error_msg.txt')
        break
    elif (status_run(tmpDir) == 'R_ExecutionHalted'):
        issue_echo(status_run(tmpDir), tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('ERROR!!!\n')
        writeErrorMessage('', 'R_Error_msg.txt')
        break
    elif (status_run(tmpDir) == 'R_Other_Error'):
        issue_echo(status_run(tmpDir), tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('ERROR!!!\n')
        writeErrorMessage('', 'R_Error_msg.txt')
        break
    elif (status_run(tmpDir) == 'R_Our_Error'):
        issue_echo(status_run(tmpDir), tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('ERROR!!!\n')
        writeErrorMessage('', 'R_Error_msg.txt')
        break

    elif master_out_of_time(time_start):
        issue_echo('Master out of time', tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('ERROR!!!\n')
        writeErrorMessage('Master out of time', 'NULL')
        break

    elif did_run_out_of_time(tmpDir, R_MAX_time):
        issue_echo('R run out of time', tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('ERROR!!!\n')
        writeErrorMessage('R run out of time', 'NULL')
        break
    
    elif did_R_crash_in_slaves(tmpDir, machine_root = 'karl')[0]:
        issue_echo('R crash in slaves', tmpDir)
        cleanups(tmpDir, newDir, lamSuffix)
        writeStatusServer('ERROR!!!\n')
        writeErrorMessage('R crash in slaves\n ' +
                          'See ' +
                          did_R_crash_in_slaves(tmpDir, machine_root = 'karl')[1],
                          'NULL')
        break

    elif did_mpi_crash(tmpDir, machine_root = 'karl'):
        count_mpi_crash += 1
        add_to_MPIErrorLog('ADaCGH-server',
                           tmpDir, socket.gethostname(),
                           message = 'MPI crash')
        if count_mpi_crash > MAX_MPI_CRASHES:
            logMPIerror(tmpDir, MAX_MPI_CRASHES)
            issue_echo('count_mpi_crash > MAX_MPI_CRASHES', tmpDir)
            cleanups(tmpDir, newDir, lamSuffix)
            writeStatusServer('ERROR!!!\n')
            writeErrorMessage('More MPI crashes than allowed', 'NULL')
            break
        else:
            recover_from_lam_crash(tmpDir, NCPU, MAX_NUM_PROCS,
                                   lamSuffix,
                                   machine_root = 'karl')
    else:
        lam_crash_log(tmpDir, 'NoCrash') ## if we get here, this much we know
    time.sleep(TIME_BETWEEN_CHECKS)


issue_echo('before LAM clean up', tmpDir)

## Clean up lam systemwide (a favor for future runs)
os.system('/http/mpi.log/killOldLamAllMachines.py')

issue_echo('at the very end!', tmpDir)





