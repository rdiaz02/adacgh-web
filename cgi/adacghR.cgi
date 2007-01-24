#!/usr/bin/python

####  Copyright (C)  2005-2006, Ramon Diaz-Uriarte <rdiaz02@gmail.com>

#### This program is free software; you can redistribute it and/or
#### modify it under the terms of the Affero General Public License
#### as published by the Affero Project, version 1
#### of the License.

#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### Affero General Public License for more details.

#### You should have received a copy of the Affero General Public License
#### along with this program; if not, you can download if
#### from the Affero Project at http://www.affero.org/oagpl.html




import socket
import sys
import os
import cgi 
import types
import time
import shutil
import dircache
import string
import whrandom
import re
from stat import ST_SIZE
import cgitb
cgitb.enable() ## zz: eliminar for real work?
sys.stderr = sys.stdout

MAX_adacgh =  25 ## MAX_genesrf + 1 = Maximum number of R processes running 
## at same time.
MAX_time = 3600 * 24 * 5 ## 5 is days until deletion of a tmp directory
R_MAX_time = 3600 * 12 ## 8 hours is max duration allowed for any process
MAX_covariate_size = 363948523L ## a 500 * 40000 array of floats
#MAX_time_size = 61897L
##  f5 <- rep(paste(paste(letters, collapse = ""),
##                  paste(LETTERS, collapse="")), 1000)
## so each of 1000 labels has 48 chars.

acceptedMethodaCGH = ('Wavelets', 'PSW', 'DNAcopy', 'ACE', 'GLAD', 'HMM', 'BioHMM', 'CGHseg')
acceptedMethodCentering = ('Median', 'Mean', 'None')
acceptedIDTypes = ('None', 'cnio', 'affy', 'clone', 'acc', 'ensembl', 'entrez', 'ug')
acceptedOrganisms = ('None', 'Hs', 'Mm', 'Rn')


def commonOutput():
    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    <title>ADaCGH</title>
    </head>
    <body>
    """

## For redirections, from Python Cookbook
def getQualifiedURL(uri = None):
    """ Return a full URL starting with schema, servername and port.
        *uri* -- append this server-rooted uri (must start with a slash)
    """
    schema, stdport = ('http', '80')
    host = os.environ.get('HTTP_HOST')
    if not host:
        host = os.environ.get('SERVER_NAME')
        port = os.environ.get('SERVER_PORT', '80')
        if port != stdport: host = host + ":" + port
    result = "%s://%s" % (schema, host)
    if uri: result = result + uri
    return result

def getScriptname():
    """ Return te scriptname part of the URL."""
    return os.environ.get('SCRIPT_NAME', '')

def getBaseURL():
    """ Return a fully qualified URL to this script. """
    return getQualifiedURL(getScriptname())

def fileUpload(fieldName):
    """Upload and get the files and do some checking. We assume there is an existing call
    to fs = cgi.FieldStorage()"""
## we don't deal with OS specific "\n"
## because R does not have a problem (at least with Windows files)
## no problem in R either with empty carriage returns at end of file
    
    if fs.has_key(fieldName):
        fileClient = fs[fieldName].file
        if not fileClient:
            shutil.rmtree(tmpDir)
            commonOutput()
            print "<h1> ADaCGH ERROR </h1>"    
            print "<p> The ", fieldName, "file you entered is not a file </p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()
    else:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", fieldName, "file required </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
            
    # transferring files to final destination;

    fileInServer = tmpDir + "/" + fieldName
    srvfile = open(fileInServer, mode = 'w')
    fileString = fs[fieldName].value
    srvfile.write(fileString)
    srvfile.close()

    os.chmod(fileInServer, 0666)
        
    if os.path.getsize(fileInServer) == 0:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"
        print "<p>", fieldName, " file has size 0 </p>"
        print "<p> Please enter a file with something in it.</p>"
        print "</body></html>"
        sys.exit()



def valueNumUpload(fieldName, testNumber = 'float', minValue = 0, errorMessageName = ''):
    """Upload and get the values and do some checking. For text and radio selections
    with positive numeric data.
    We assume there is an existing call to fs = cgi.FieldStorage()"""
    if errorMessageName == '': errorMessageName = fieldName
    if not fs.has_key(fieldName):
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", errorMessageName, "value required </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    if fs[fieldName].filename:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", errorMessageName, "should not be a file. </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    if type(fs[fieldName]) == type([]):
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", errorMessageName, "should be a single value.</p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    else:
        tmp = fs[fieldName].value
    ## Accept only numeric values that can be turned to floats or ints
    if testNumber == 'float':
        try:
            tmpn = float(tmp)
        except:
            commonOutput()
            print "<h1> ADaCGH ERROR </h1>"    
            print "<p> ", errorMessageName, "is not a valid numeric value.</p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()
    else:
        try:
            tmpn = int(tmp)
        except:
            commonOutput()
            print "<h1> ADaCGH ERROR </h1>"    
            print "<p> ", errorMessageName, "is not a valid numeric value.</p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()
    if tmpn < minValue:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", errorMessageName, "smaller than smallest accepted value (", minValue, "). </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    ## transferring files to final destination;
    fileInServer = tmpDir + "/" + fieldName
    srvfile = open(fileInServer, mode = 'w')
    srvfile.write(str(tmpn))
    srvfile.close()
    os.chmod(fileInServer, 0666)
    return tmpn




def radioUpload(fieldName, acceptedValues):
    """Upload and get the values and do some checking. For radio selections
    with text data; check those are in acceptedValues.
    We assume there is an existing call to fs = cgi.FieldStorage()"""

    if not fs.has_key(fieldName):
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p>", fieldName, "required </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    if fs[fieldName].filename:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", fieldName, "should not be a file. </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    if type(fs[fieldName]) == type([]):
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p>", fieldName, "should be a single value.</p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    else:
        tmp = fs[fieldName].value
            
    if tmp not in acceptedValues:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> The", fieldName, "choosen is not valid.</p>"
        print "<p> Please fill up the required fields and try again.</p>"
        print "</body></html>"
        sys.exit()

    fileInServer = tmpDir + "/" + fieldName
    srvfile = open(fileInServer, mode = 'w')
    fileString = tmp
    srvfile.write(fileString)
    srvfile.close()
    os.chmod(fileInServer, 0666)

    return tmp




#########################################################
#########################################################

####          Execution starts here      ################

#########################################################
#########################################################



## Deleting tmp directories older than MAX_time
currentTime = time.time()
currentTmp = dircache.listdir("/http/adacgh/www/tmp")
for directory in currentTmp:
    tmpS = "/http/adacgh/www/tmp/" + directory
    if (currentTime - os.path.getmtime(tmpS)) > MAX_time:
        shutil.rmtree(tmpS)


### Creating temporal directories
newDir = str(whrandom.randint(1, 10000)) + str(os.getpid()) + str(whrandom.randint(1, 100000)) + str(int(currentTime)) + str(whrandom.randint(1, 10000))
redirectLoc = "/tmp/" + newDir
tmpDir = "/http/adacgh/www/tmp/" + newDir
os.mkdir(tmpDir)
os.chmod(tmpDir, 0700)

### File and parameter upload
fs = cgi.FieldStorage()

idtype = radioUpload('idtype', acceptedIDTypes)
organism = radioUpload('organism', acceptedOrganisms)
tmp = valueNumUpload('MCR.gapAllowed', 'float', 1)
tmp1 = valueNumUpload('MCR.alteredLow', 'float', 0)
tmp2 = valueNumUpload('MCR.alteredHigh', 'float', 0)
tmp = valueNumUpload('MCR.recurrence', 'float', 1)

if tmp1 >= tmp2:
    commonOutput()
    print "<h1> ERROR </h1>"
    print "<p> It makes no sense to have an alteredLow larger or equal to \
    an alteredHigh </p>"
    print "</body></html>"
    sys.exit()
    






methodaCGH = radioUpload('methodaCGH', acceptedMethodaCGH)
if methodaCGH == 'Wavelets':
    tmp = valueNumUpload('Wave.minDiff', 'float', 0)
    tmp = radioUpload('Wave.merge', ('Yes','No'))
if methodaCGH == 'ACE':
    tmp = valueNumUpload('ACE.fdr', 'float', 0)
if methodaCGH == 'PSW':
    tmp = valueNumUpload('PSW.nIter', 'int', 10)
    tmp = valueNumUpload('PSW.p.crit', 'float', 0)
centering = radioUpload('centering', acceptedMethodCentering)
twofiles = radioUpload('twofiles', ('Two.files', 'One.file'))


if twofiles == 'Two.files':
    ##check if file coming from preP

    if(fs.getfirst("covariate2")!= None):
        prep_tmpdir = fs.getfirst("covariate2")
        shutil.copy("/http/prep/www/tmp/" + prep_tmpdir +"/outdata.txt",tmpDir + "/acghData")
    else:
        fileUpload('acghData')
        if os.stat(tmpDir + '/acghData')[ST_SIZE] > MAX_covariate_size:
            shutil.rmtree(tmpDir)
            commonOutput()
            print "<h1> ADaCGH ERROR </h1>"
            print "<p> aCGH file way too large </p>"
            print "<p> aCGH files this size not allowed.</p>"
            print "</body></html>"
            sys.exit()
    fileUpload('positionInfo')
    if os.stat(tmpDir + '/positionInfo')[ST_SIZE] > MAX_covariate_size:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"
        print "<p> Position info. file way too large </p>"
        print "<p> Position info. files this size not allowed.</p>"
        print "</body></html>"
        sys.exit()
elif twofiles == 'One.file':
    fileUpload('acghAndPosition')
    if os.stat(tmpDir + '/acghAndPosition')[ST_SIZE] > MAX_covariate_size:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"
        print "<p> aCGH file way too large </p>"
        print "<p> aCGH files this size not allowed.</p>"
        print "</body></html>"
        sys.exit()





## Upload worked OK. We store the original names of the files in the
## browser for later report:
# fileNamesBrowser = open(tmpDir + '/fileNamesBrowser', mode = 'w')
# fileNamesBrowser.write(fs['covariate'].filename + '\n')
# fileNamesBrowser.write(fs['time'].filename + '\n')
# fileNamesBrowser.write(fs['event'].filename + '\n')
# fileNamesBrowser.close()


## current number of processes > max number of processes?
## and we do it here, not before, so that we have the most
## current info about number of process right before we launch R.


## First, delete any R file left (e.g., from killing procs, etc).
RrunningFiles = dircache.listdir("/http/adacgh/www/R.running.procs")
for Rtouchfile in RrunningFiles:
    tmpS = "/http/adacgh/www/R.running.procs/" + Rtouchfile
    if (currentTime - os.path.getmtime(tmpS)) > R_MAX_time:
        os.remove(tmpS)

## Now, verify any processes left
numRadacgh = len(dircache.listdir("/http/adacgh/www/R.running.procs"))
if numRadacgh > MAX_adacgh:
    shutil.rmtree(tmpDir)
    commonOutput()
    print "<h1> ADaCGH problem: The servers are too busy </h1>"
    print "<p> Because of the popularity of the application "
    print " the maximum number of simultaneous runs of ADaCGH has been reached.</p>"
    print "<p> Please try again later.</p>"
    print "<p> We apologize for the inconvenience.</p>"    
    print "</body></html>"
    sys.exit()
    

################        Launching R   ###############


## what would I need this for????
    ## To have a simple way  of loading array names


# prepare the arrayNames file:
if twofiles == 'Two.files':
    covarInServer = tmpDir + "/acghData"
elif twofiles == 'One.file':
    covarInServer = tmpDir + "/acghAndPosition"
geneNames = tmpDir + "/geneNames"
arrayNames = tmpDir + "/arrayNames"

srvfile = open(covarInServer, mode = 'rU')
arrayfile = open(arrayNames, mode = 'w')
genenfile = open(geneNames, mode = 'w')

covarR = open(tmpDir + '/covarR', mode = 'w')

all_covar_lines = srvfile.readlines()
srvfile.close()

read_at_least_one = False

if twofiles == 'One.file':
    positions_to_remove = 4
    positionInfo = tmpDir + "/positionInfo"
    positionfile = open(positionInfo, mode = 'w')
else:
    positions_to_remove = 1

arraynames =[]
num_name_lines = 0
for nr in range(0, len(all_covar_lines)):
    line = all_covar_lines[nr]
    if (line.find("#name") == 0) or (line.find("#NAME") == 0) or (line.find("#Name") == 0) \
           or (line.find('"#name"') == 0) or (line.find('"#NAME"') == 0) or (line.find('"#Name"') == 0):
        num_name_lines = num_name_lines + 1
        if num_name_lines > 1:
            commonOutput()
            print """ You have more than one line with #Name (or #NAME or #name), in the data matrix \
                   but only one is allowed."""
            sys.exit()
        arraynames = line.split('\t', positions_to_remove)[-1]
    elif (line.find("#", 0, 1) == 0):
        continue
    else:
        line_splitted = line.split('\t', positions_to_remove)   
        if line_splitted[-1].find(',') >= 0:
            commonOutput()
            print "<h1> ERROR </h1>"
            print """<p> You have commas in the data matrix.
            You probably are using commas instead of "." for the
            decimal separator. Please, use a "." instead of a
            ",".</p>
            """
            covarR.close()
            genenfile.close()
            sys.exit()
        else:
            gname = line_splitted[0]
            if twofiles == 'One.file':
                if len(line_splitted) < 5:
                    commonOutput()
                    print "<h1> ERROR </h1>"
                    print """<p> There are no data after the four position fields
                    (name, chromosome, start and end). Please check your data.</p>
                    """
                    covarR.close()
                    genenfile.close()
                    sys.exit()

                positionfile.write(''.join([line_splitted[0], '\t',
                                            line_splitted[1], '\t',
                                            line_splitted[2], '\t',
                                            line_splitted[3], '\n']))
            line_splitted = line_splitted[-1]   
            if line_splitted.rstrip():
                covarR.write(line_splitted)
                line_splitted = line_splitted.split('\t')   
                num_subjects = len(line_splitted)
                if read_at_least_one == False:
                    previous_num_subjects = num_subjects
                    read_at_least_one = True
                if num_subjects != previous_num_subjects:
                    commonOutput()
                    print "<h1>  <i>pre</i>P ERROR </h1>"
                    print """<p> Not all rows have the same number of
                    columns (subjects or arrays). We first detect this
                    discrepancy in rows """
                    print nr - 1
                    print "and"
                    print nr
                    print "</p>"
                    covarR.close()
                    genenfile.close()
                    sys.exit()
                else:
                    genenfile.write(gname + '\n')
                    previous_num_subjects = num_subjects

if not len(arraynames):
    arraynames = '\t'.join([''.join(['array.', str(i)]) for i in range(1, num_subjects + 1)])
arrayfile.write(arraynames)
arrayfile.write("\n\n")
arrayfile.close()
os.chmod(arrayNames, 0600)


covarR.close()
genenfile.close()
arrayfile.close()
if twofiles == "One.file":
    positionfile.close()

## Get rid of quotes and double quotes:
tmp = os.system("cd " + tmpDir +  """; sed 's/"//g' geneNames > tmp234; mv tmp234 geneNames""")
tmp = os.system("cd " + tmpDir + """; sed "s/'//g" geneNames > tmp234; mv tmp234 geneNames""")
tmp = os.system("cd " + tmpDir +  """; sed 's/"//g' positionInfo > tmp234; mv tmp234 positionInfo""")
tmp = os.system("cd " + tmpDir + """; sed "s/'//g" positionInfo > tmp234; mv tmp234 positionInfo""")
tmp = os.system("cd " + tmpDir +  """; sed 's/"//g' arrayNames > tmp234; mv tmp234 arrayNames""")
tmp = os.system("cd " + tmpDir + """; sed "s/'//g" arrayNames > tmp234; mv tmp234 arrayNames""")




## It would be good to use spawnl or similar instead of system,
## but I have no luck with R. This, I keep using system.
## Its safety depends crucially on the newDir not being altered,
## but newDir is not passed from any other user-reachable place
## (it is created here).


## recall to include in R
    ##pid <- Sys.getpid()
    ##write.table(file = "pid.txt", pid, row.names = FALSE, col.names = FALSE)

## touch Rout, o.w. checkdone can try to open a non-existing file
touchRout = os.system("/bin/touch " + tmpDir + "/f1.Rout") 
touchRrunning = os.system("/bin/touch /http/adacgh/www/R.running.procs/R." + newDir +
                          "@" + socket.gethostname())
shutil.copy("/http/adacgh/cgi/f1.R", tmpDir)
## we add the 2> error.msg because o.w. if we kill R we get a server error as standard
## error is sent to the server
# Rcommand = "cd " + tmpDir + "; " + "/usr/bin/R CMD BATCH --no-restore --no-readline --no-save -q f1.R 2> error.msg &"
# Rrun = os.system(Rcommand)
##os.system('cp /http/js/js.adacgh.squeleton1.html ' + tmpDir + '/.')
##os.system('cp /http/js/js.adacgh.squeleton2.html ' + tmpDir + '/.')
##os.system('cp /http/js/final.adacgh.squeleton.html ' + tmpDir + '/.')

tryrrun = os.system('/http/mpi.log/tryRrun2.py ' + tmpDir +' 10 ' + 'ADaCGH &')
createResultsFile = os.system("/bin/touch " + tmpDir + "/results.txt")



###########   Creating a results.hmtl   ###############

## Copy to tmpDir a results.html that redirects to checkdone.cgi
## If communication gets broken, there is always a results.html
## that will do the right thing.
shutil.copy("/http/adacgh/cgi/results-pre.html", tmpDir)
os.system("cd " + tmpDir + "; /bin/sed 's/sustituyeme/" +
          newDir + "/g' results-pre.html > results.html; rm results-pre.html")

##############    Redirect to checkdone.cgi    ##################
print "Location: "+ getQualifiedURL("/cgi-bin/checkdone.cgi") + "?newDir=" + newDir, "\n\n"


