#!/usr/bin/python

import sys
import os
import cgi 
import types
import time
import shutil
import re
import tarfile
import glob
from stat import ST_SIZE

import cgitb
cgitb.enable() ## zz: eliminar for real work?
sys.stderr = sys.stdout


ACE_MAX_time = 2 * 3600 ## 4 hours is max duration allowd for any process


### There are lots of ugly hacks here. Could turn into a module many of the
### functions below _if_ they took parameters (instead of just getting things
### from the global environment.

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

def relaunchCGI():
    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    """
    print '<meta http-equiv="Refresh"'
    print 'content="30; URL=' + getBaseURL() + '?newDir=' + newDir + '">'
    print '<title>ADaCGH results</title>'
    print '</head> <body>'
    print '<p> This is an autorefreshing page; your results will eventually be displayed here.\n'
    print 'If your browser does not autorefresh, the results will be kept for five days at</p>'
    print '<p><a href="' + getBaseURL() + '?newDir=' + newDir + '">', 'http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html</a>.' 
    print '</p> </body> </html>'



def valueNumUpload(fieldName, testNumber = 'float', minValue = 0):
    """Upload and get the values and do some checking. For text and radio selections
    with positive numeric data.
    We assume there is an existing call to fs = cgi.FieldStorage()"""
    if not fs.has_key(fieldName):
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", fieldName, "value required </p>"
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
        print "<p> ", fieldName, "should be a single value.</p>"
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
            print "<p> ", fieldName, "is not a valid numeric value.</p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()
    else:
        try:
            tmpn = int(tmp)
        except:
            commonOutput()
            print "<h1> ADaCGH ERROR </h1>"    
            print "<p> ", fieldName, "is not a valid numeric value.</p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()
    if tmpn < minValue:
        shutil.rmtree(tmpDir)
        commonOutput()
        print "<h1> ADaCGH ERROR </h1>"    
        print "<p> ", fieldName, "smaller than smallest accepted value (", minValue, "). </p>"
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



def printErrorRun():
    if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
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
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


def printOKRun():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>ADaCGH results</title></head><body>\n")

    if os.path.exists(tmpDir + "/ErrorFigure.png"):
        outf.write('<IMG WIDTH="500" HEIGHT="500" BORDER="0" SRC="ErrorFigure.png">')
        outf.write("<br /><br /> <hr>")
        outf.write("<pre>")
        outf.write('<br /><br /><h2> Results <a href="http://adacgh2.bioinfo.cnio.es/help/adacgh-help.html#outputText">(help)</a></h2> \n')
        outf.write("<br /><br /> <hr>")
        outf.write(cgi.escape(resultsFile))
        outf.write("</pre>")
        outf.write("</body></html>")
        outf.close()
        Rresults.close()
        shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
    else:
        outf.write('<br /><br /><center><h2> ADaCGH Results <a href="http://adacgh2.bioinfo.cnio.es/help/adacgh-help.html#outputText">(help)</a></center></h2> \n')
        outf.write(open(tmpDir + "/results.for.html").read())
        outf.write('<br />')
        outf.write('<hr>')
        outf.write('<h3> Genes/clones per chromosome </h3>')
        outf.write(open(tmpDir + "/clones.per.chrom.html").read())
        outf.write('<br />')
        outf.write('<h3> Summary statistics before centering </h3>')
        outf.write(open(tmpDir + "/stats.before.centering.html").read())
        outf.write('<br />')
        outf.write('<h3> Summary statistics: subject/array by chromosome before centering </h3>')
        outf.write('<h4> Mean </h4>')
        outf.write(open(tmpDir + "/stats.subj.by.chrom.mean.BEFORE.html").read())
        outf.write('<br />')
        outf.write('<h4> Median </h4>')
        outf.write(open(tmpDir + "/stats.subj.by.chrom.median.BEFORE.html").read())
        outf.write('<br />')
        outf.write('<h4> MAD </h4>')
        outf.write(open(tmpDir + "/stats.subj.by.chrom.mad.BEFORE.html").read())
        outf.write('<hr>')
        outf.write('<br />')
        outf.write('<h3> Summary statistics after centering </h3>')
        outf.write(open(tmpDir + "/stats.after.centering.html").read())
        outf.write('<br />')
        outf.write('<h3> Summary statistics: subject/array by chromosome after centering </h3>')
        outf.write('<h4> Mean </h4>')
        outf.write(open(tmpDir + "/stats.subj.by.chrom.mean.AFTER.html").read())
        outf.write('<br />')
        outf.write('<h4> Median </h4>')
        outf.write(open(tmpDir + "/stats.subj.by.chrom.median.AFTER.html").read())
        outf.write('<br />')
        outf.write('<h4> MAD </h4>')
        outf.write(open(tmpDir + "/stats.subj.by.chrom.mad.AFTER.html").read())
        outf.write('<br />')
        ## The following is common to all
        allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
        allResults.add(tmpDir + '/results.txt', 'summary.statistics.txt')

        outf.write('<h2>FDR table</h2>')
        acefdrtable = open(tmpDir + "/ace.fdrtable.html")
        acefdr = acefdrtable.read()
        acefdrtable.close()
        outf.write(acefdr)
        outf.write('<form action="http://adacgh2.bioinfo.cnio.es/cgi-bin/ace.cgi" method="GET">\n')
        outf.write('<input type="hidden" NAME="newDir" VALUE="' + newDir + '">')
        currentfdr = str(open(tmpDir + '/aceFDR').readline())
        outf.write('<br><input type="TEXT" name="fdrace" value="' +
                   currentfdr + '" size="10"  maxlength="10">\n')
        outf.write('<input value="Submit" type = "SUBMIT"> (Change the desired FDR and Press "Submit" to obtain figures with new FDR)')
        outf.write('<h2>Segmented plots</h2><p>Click on thumbnails to expand.</p>')
        thumb(tmpDir,  open(tmpDir + '/arrayNames', mode = 'r').read().split('\n')[0].split('\t'), outf, maxsthumb = 350)
        thumb(tmpDir, ['All_arrays'], outf, maxsthumb = 350)
        outf.write('<p>Inferred gains and losses available from file' +
                   '<a href="./ACE.results.FDR=' + currentfdr + '.txt">' +
                   '"ACE.results.FDR=' + currentfdr + '.txt"</a></p>')
        if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
        if os.path.exists(tmpDir + '/f1.R'): os.remove(tmpDir + '/f1.R')
        if os.path.exists(tmpDir + '/rerunACE.R'): os.remove(tmpDir + '/rerunACE.R')
        if os.path.exists(tmpDir + '/f1.Rout'): os.remove(tmpDir + '/f1.Rout')
        #if os.path.exists(tmpDir + '/.RData'): os.remove(tmpDir + '/.RData')
        allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
        os.chdir(tmpDir)
        ll1 = glob.glob('*.log')
        for dname in ll1:
            os.remove(dname)
        lll = glob.glob('*')
        for flname in lll:
            allResults.add(flname)
        allResults.close()
        outf.write('<hr> <a href="http://adacgh2.bioinfo.cnio.es/tmp/' +
                   newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
        outf.write("</body></html>")
        outf.close()
        Rresults.close()
        shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


def printRKilled():
    if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
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


def thumb(tmpDir, fnames, outf, maxsthumb = 350):
    """ From a set of pngs, obtain thumbnails and
    add links to html. The thumbnails are
    inserted in the webpage, the large png viewed on clik.
    tmpDir is the directory where the files live,
    fnames is a list with the base file names to process
    maxsthumb   max size of thumbnail and outf the html file.
    """
    mst = str(maxsthumb)
    os.chdir(tmpDir)
    for bname in fnames:
        os.system(''.join(['/usr/bin/convert ', bname, '.png',
                           ' -resize ', mst, 'x', mst, ' thumb.', 
                           bname, '.jpeg']))
        outf.write(''.join(['<a href="', bname, '.html"> <img alt="',
	                    bname, '" title="', bname, '" src="thumb.',
                            bname, '.jpeg"></a>']))
    os.chdir('/http/adacgh2/cgi')





############################################################

######     Starting execution here

############################################################

## Logic: it is launched from results.html. Then,, it can be relaunched
    ## indefinite number of time. When launched, verifies if a
    ## rerunACE.Rout exists. If iy does, then an R process to do
    ## figures is running, so verify if done. Otherwise,
    ## launch it, and call back to myself.

fs = cgi.FieldStorage()
if fs.has_key('newDir'):
    value = fs['newDir']
    if type(value) is types.ListType:
        commonOutput()
        print "<h1> ERROR </h1>"    
        print "<p> newDir should not be a list. </p>"
        print "<p> Anyone trying to mess with it?</p>"
        print "</body></html>"
        sys.exit()
    else:
        newDir = value.value
else:
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir is empty. </p>"
    print "</body></html>"
    sys.exit()
            
if re.search(r'[^0-9]', str(newDir)):
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir does not have a valid format. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()
                
##    redirectLoc = "/tmp/" + tmpDir
tmpDir = "/http/adacgh2/www/tmp/" + newDir
if not os.path.isdir(tmpDir):
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> tmpDir is not a valid directory. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()



if not fs.has_key('fdrace'): ## this will only happen when coming from the refresh calls
    fdrace = -99
else:
    fdrace = valueNumUpload('fdrace', testNumber = 'float', minValue = 0)

if fdrace > 1:
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> The FDR cannot be larger than 1. </p>"
    print "</body></html>"
    sys.exit()
                           
if not os.path.exists(tmpDir + '/rerunACE.Rout'):
    shutil.copy("/http/adacgh2/cgi/rerunACE.R", tmpDir)
    Rcommand = "cd " + tmpDir + "; " + "/http/R-custom/bin/R CMD BATCH --no-restore --no-readline --no-save -q rerunACE.R 2> errorACE.msg &"
    Rrun = os.system(Rcommand)

###########   Creating a results.hmtl   ###############
    
## Copy to tmpDir a results.html that redirects to checkdone.cgi
## If communication gets broken, there is always a results.html
## that will do the right thing.
    shutil.copy("/http/adacgh2/cgi/results-pre-ace.html", tmpDir)
    os.system("cd " + tmpDir + "; /bin/sed 's/sustituyeme/" +
              newDir + "/g' results-pre-ace.html > results.html; rm results-pre-ace.html")

##############    Redirect to ace.cgi    ##################
    print "Location: "+ getQualifiedURL("/cgi-bin/ace.cgi") + "?newDir=" + newDir, "\n\n"

else: ## i.e., we have already started the R process
    Rrout = open(tmpDir + "/rerunACE.Rout")
    soFar = Rrout.read()
    Rrout.close()
    finishedOK = soFar.endswith("Normal ACE termination\n")
    errorRun = soFar.endswith("Execution halted\n")

    if os.path.exists(tmpDir + "/ACEpid.txt"):
    ## do we need to kill an R process?
        if (time.time() - os.path.getmtime(tmpDir + "/ACEpid.txt")) > ACE_MAX_time:
            try:
                os.kill(int(open(tmpDir + "/ACEpid.txt", mode = "r").readline()),
                        signal.SIGINT) ## maybe sigint is better than sigkill??  
            finally:
                printRKilled()
                os.rename(tmpDir + '/ACEpid.txt', tmpDir + '/killed.ACEpid.txt')
                if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
                if os.path.exists(tmpDir + '/rerunACE.R'): os.remove(tmpDir + '/rerunACE.R')
                print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'
                sys.exit()

    if errorRun > 0:
        printErrorRun()
        if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
        if os.path.exists(tmpDir + '/rerunACE.R'): os.remove(tmpDir + '/rerunACE.R')
        try: os.rename(tmpDir + '/ACEpid.txt', tmpDir + '/natural.death.ACEpid.txt')
	except: None
        try: os.remove(tmpDir + '/rerunACE.R')
        except: None
        print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'

    elif finishedOK > 0:
        printOKRun()
        if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
        if os.path.exists(tmpDir + '/rerunACE.R'): os.remove(tmpDir + '/rerunACE.R')
        try: os.rename(tmpDir + '/ACEpid.txt', tmpDir + '/natural.death.ACEpid.txt')
	except: None
        try: os.remove(tmpDir + '/f1.R')
        except: None
        print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n' 
    
    else:
    ## we only end up here if: we were not done in a previous run AND no process was overtime 
    ## AND we did not just finish. So we must continue.
        relaunchCGI()
    





