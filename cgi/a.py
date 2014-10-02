#!/usr/bin/python

## All this code is copyright Ramon Diaz-Uriarte. For security reasons, this is for
## now confidential. No license is granted to copy, distribute, or modify it.
## Once everything is OK, it will be distributed under the GPL.


import sys
import os
import cgi 
import types
import time
import shutil
import string
import signal
import re
import glob
import tarfile

import cgitb
cgitb.enable() ## zz: eliminar for real work?
sys.stderr = sys.stdout ## eliminar?

R_MAX_time = 8 * 3600 ## 4 hours is max duration allowd for any process

## For redirections, from Python Cookbook



def pdf2html(rootname, tmpDir, outf, compressedFile, maxsthumb = 350):
    """ From a multipage pdf obtain jpegs; the thumbnails are
    inserted in the webpage, the large jpeg viewed on clik.
    rootname is all the stuff before the `pdf',
    tmpDir is the directory where the files live, maxsthumb
    max size of thumbnail and outf the html file.
    We also decrease the size of the jpeg for showing.
    Finally, we add the generated jpeg to the compressed file"""

    mst = str(maxsthumb)
    mst2 = str(1600)
    os.chdir(tmpDir)
    os.system('/usr/bin/pdftoppm ' + rootname + '.pdf tmpppms')
    tmps = glob.glob('tmpppms*.ppm')
    for fignum in range(len(tmps)):
        os.system('/usr/bin/ppmtojpeg ' + tmps[fignum] + ' > ' + rootname + '.'
                  + str(fignum + 1) + '.jpeg')
        os.system('/usr/bin/convert -size ' + mst + 'x' +
                  mst + ' ' + rootname + '.' + str(fignum + 1) + '.jpeg' +
                  ' -resize ' + mst + 'x' + mst + ' thumb.' + rootname + '.'
                  + str(fignum + 1) + '.jpeg')
        os.system('/usr/bin/convert ' + rootname + '.' + str(fignum + 1) + '.jpeg' +
                  ' -resize ' + mst2 + 'x' + mst2 + ' ' + rootname + '.'
                  + str(fignum + 1) + '.jpeg')
        outf.write('<a href="' + rootname + '.' + str(fignum + 1) +
                   '.jpeg"> <img src="' + 'thumb.' + rootname + '.'
                  + str(fignum + 1) + '.jpeg"></a>\n')
        compressedFile.add(rootname + '.' + str(fignum + 1) + '.jpeg',
                           rootname + '.' + str(fignum + 1) + '.jpeg')
    os.chdir('/asterias-web-apps/adacgh2/cgi')


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


# def getPathinfo():
#     """ Return the remaining part of the URL. """
#     pathinfo = os.environ.get('PATH_INFO', '')
#     return pathinfo

def getBaseURL():
    """ Return a fully qualified URL to this script. """
    return getQualifiedURL(getScriptname())


def commonOutput():
    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    <title>ADaCGH results</title>
    </head>
    <body>
    """
    
## to keep executing myself:
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
    print '<p><a href="' + getBaseURL() + '?newDir=' + newDir + '">', 'http://adacgh2.iib.uam.es/tmp/'+ newDir + '/results.html</a>.' 
    print '</p> </body> </html>'
    

## Output-generating functions
def printErrorRun():
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
        outf.write('<br /><br /><h2> Results <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputText">(help)</a></h2> \n')
        outf.write("<br /><br /> <hr>")
        outf.write(cgi.escape(resultsFile))
        outf.write("</pre>")
        outf.write("</body></html>")
        outf.close()
        Rresults.close()
        shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
    else:
        outf.write('<br /><br /><center><h2> ADaCGH Results <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputText">(help)</a></center></h2> \n')
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


        methodUsed = open(tmpDir + '/methodaCGH').read()
        if (methodUsed == 'CBS') or (methodUsed == 'CBS\n'):
            outf.write('<h2>Diagnostic plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputCBS">(help)</a></h2> \n')
            outf.write('<h3>One plot per array/sample</h3>')
            outf.write('<a href="CBS.diagnostic.plots.pdf">View/save</a> the (multipage) pdf')
            outf.write('<h3>One plot per array/sample and chromosome</h3>')
            outf.write('<a href="CBS.diagnostic.plots.per.array.and.chromosome.pdf">View/save</a> the (multipage) pdf')
            outf.write('<p>(For the rest of the figures,  click on thumbnails to expand.' +
                       ' The larger figures shown are still jpegs; you can get the' +
                       ' much higher quality ones downloading the results, ' +
                       'with output also in pdf format.)</p>')

            outf.write('<h2>Segmented data plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputCBS">(help)</a></h2> \n')
	    pdf2html('CBS.segmented.plots', tmpDir, outf, allResults, 350)
            outf.write('<p>Smoothed values for all genes/clones are available from file "Wavelets.results.txt", below</p>')
            outf.write('<br />')
            outf.write('<h2>Plateau plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputCBS">(help)</a></h2> \n')
            pdf2html('CBS.plateau.plots', tmpDir, outf, allResults, 150)
            outf.write('<br />')
            allResults.add(tmpDir + '/CBS.results.txt', 'CBS.results.txt')
            allResults.add(tmpDir + '/CBS.diagnostic.plots.pdf',
                           'CBS.diagnostic.plots.pdf')
            allResults.add(tmpDir + '/CBS.diagnostic.plots.per.array.and.chromosome.pdf', 
                            'CBS.diagnostic.plots.per.array.and.chromosome.pdf')

            allResults.close()
            outf.write('<hr> <a href="http://adacgh2.iib.uam.es/tmp/' +
                       newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
            outf.write("</body></html>")
            outf.close()
            Rresults.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")

        if (methodUsed == 'WS') or (methodUsed == 'WS\n'):
            outf.write('<h2>Diagnostic autocorrelation plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputWS">(help)</a></h2> \n')
            outf.write('<a href="Autocorrelation.plots.pdf">View/save</a> the (multipage) pdf') ##zz: later, provide thumbnails
            ## and true images
            outf.write('<br />')
            outf.write('<p>(For the rest of the figures,  click on thumbnails to expand.' +
                       ' The larger figures shown are still jpegs; you can get the' +
                       ' much higher quality ones downloading the results, ' +
                       'with output also in pdf format.)</p>')
            outf.write('<h2>Segmented data plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputWS">(help)</a></h2> \n')
            pdf2html('WS.segmented.plots', tmpDir, outf, allResults, 350)
            outf.write('<p>Smoothed values for all genes/clones are available from file "Wavelets.results.txt", below</p>')
            outf.write('<br />')
            outf.write('<h2>Plateau plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputWS">(help)</a></h2> \n')
            pdf2html('WS.plateau.plots', tmpDir, outf, allResults, 150)
            outf.write('<br />')

                
            allResults.add(tmpDir + '/Wavelets.results.txt', 'Wavelets.results.txt')
            allResults.add(tmpDir + '/Autocorrelation.plots.pdf', 
                            'Autocorrelation.plots.pdf')
                
            allResults.close()
            outf.write('<hr> <a href="http://adacgh2.iib.uam.es/tmp/' +
                       newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
            outf.write("</body></html>")
            outf.close()
            Rresults.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
  

        if (methodUsed == 'PSW') or (methodUsed == 'PSW\n'):
            outf.write('<p>(Here are the figures; click on thumbnails to expand.' +
                       ' The larger figures shown are still jpegs; you can get the' +
                       ' much higher quality ones downloading the results, ' +
                       'with output also in pdf format.)</p>')


            outf.write('<h2>Island plots, gains <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputPSW">(help)</a></h2> \n')
            pdf2html('PSW.island.plots.gains', tmpDir, outf, allResults, 350)
            outf.write('<br />')
            
            outf.write('<h2>Island plots, losses <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputPSW">(help)</a></h2> \n')
            pdf2html('PSW.island.plots.losses', tmpDir, outf, allResults, 350)
            outf.write('<br />')

            outf.write('<p>Smith-Waterman results for all genes/clones are available from files "Gains.Price.Smith.Waterman.output.txt" and "Losses.Price.Smith.Waterman.output.txt", below</p>')
                
            allResults.add(tmpDir + '/Gains.Price.Smith.Waterman.output.txt', 'Gains.Price.Smith.Waterman.output.txt')
            allResults.add(tmpDir + '/Losses.Price.Smith.Waterman.output.txt', 'Losses.Price.Smith.Waterman.output.txt')

            allResults.close()
            outf.write('<hr> <a href="http://adacgh2.iib.uam.es/tmp/' +
                       newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
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


    
## Changing to the appropriate directory
    
form = cgi.FieldStorage()
if form.has_key('newDir'):
    value = form['newDir']
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
## newDir can ONLY contain digits.
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir does not have a valid format. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()
    
redirectLoc = "/tmp/" + newDir
tmpDir = "/asterias-web-apps/adacgh2/www/tmp/" + newDir

if not os.path.isdir(tmpDir):
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir is not a valid directory. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()
    

## Were we already done in a previous execution?
## No need to reopen files or check anything else. Return url with results
## and bail out.
if os.path.exists(tmpDir + "/natural.death.pid.txt") or os.path.exists(tmpDir + "/killed.pid.txt"):
    print 'Location: http://adacgh2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'
    sys.exit()

## No, we were not done. Need to examine R output
Rrout = open(tmpDir + "/f1.Rout")
soFar = Rrout.read()
Rrout.close()
finishedOK = soFar.endswith("Normal termination\n")
errorRun = soFar.endswith("Execution halted\n")

if os.path.exists(tmpDir + "/pid.txt"):
    ## do we need to kill an R process?
    if (time.time() - os.path.getmtime(tmpDir + "/pid.txt")) > R_MAX_time:
        try:
            os.kill(int(open(tmpDir + "/pid.txt", mode = "r").readline()),
        	             signal.SIGINT) ## maybe sigint is better than sigkill??  
        finally:
            printRKilled()
            os.rename(tmpDir + '/pid.txt', tmpDir + '/killed.pid.txt')
            os.remove(tmpDir + '/f1.R')
            try:
                os.remove("/asterias-web-apps/adacgh2/www/R.running.procs/R." + newDir)
            finally:
                print 'Location: http://adacgh2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'
                sys.exit()

if errorRun > 0:
    printErrorRun()
    os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    os.remove(tmpDir + '/f1.R')
    try:
        os.remove("/asterias-web-apps/adacgh2/www/R.running.procs/R." + newDir)
    finally:
        print 'Location: http://adacgh2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'

elif finishedOK > 0:
    printOKRun()
    os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    os.remove(tmpDir + '/f1.R')
    try:
        os.remove("/asterias-web-apps/adacgh2/www/R.running.procs/R." + newDir)
    finally:
        print 'Location: http://adacgh2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n' 
    
else:
    ## we only end up here if: we were not done in a previous run AND no process was overtime 
    ## AND we did not just finish. So we must continue.
    relaunchCGI()
    



