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

R_MAX_time = 12 * 3600 ## 12 hours is max duration allowd for any process


def issue_echo(fecho, tmpDir):
    """Silly function to output small tracking files"""
    timeHuman = '##########   ' + \
                str(time.strftime('%d %b %Y %H:%M:%S')) 
    os.system('echo "' + timeHuman + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "' + fecho + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "    " >> ' + tmpDir + '/checkdone.echo')


def kill_lamcheck(pid, machine):
    'as it says: to kill lamcheck; actually, anything'
    os.system('ssh ' + machine + ' "kill -s 9 ' + pid + '"')

## For redirections, from Python Cookbook

def results_print_general(outf, tmpDir, newDir, Rresults):
    outf.write('<h2>Segmented data plots <a href="http://adacgh2.bioinfo.cnio.es/help/adacgh-help.html#output">(help)</a></h2> \n')
    thumb(tmpDir, open(tmpDir + '/arrayNames', mode = 'r').read().split('\n')[0].split('\t'),
          outf, maxsthumb = 350)
    thumb(tmpDir, ['All_arrays'], outf, maxsthumb = 350)
    outf.flush()
    output_name = glob.glob(tmpDir + '/*.output.txt')[0].split('/')[-1]
    outf.write('<p>Smoothed values for all genes/clones are available from file' +
               ' <a href="./' + output_name + '">"' + output_name + '".</a></p>')
    outf.write('<br />')
    if os.path.exists(tmpDir + "/mcr.results.html"):
        outf.write('<h2>Minimal common regions</h2>\n')
        outf.write(open(tmpDir + "/mcr.results.html").read())
    outf.write('<br />')
#     allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
    os.chdir(tmpDir)
    ll1 = glob.glob('*.log')
    for dname in ll1:
        os.remove(dname)
#     lll = glob.glob('*')
#     for flname in lll:
#         try: allResults.add(flname)
#         except: None
    allResults.close() ## this is done externally
    os.system('tar -czvf all.results.tar.gz * &')
    outf.write('<hr> <a href="http://adacgh2.bioinfo.cnio.es/tmp/' +
               newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
   
    try:
        outf.write(printPalsURLADaCGH(newDir))
    except:
        None
         
    outf.write("</body></html>")
    outf.flush()
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")





def printPalsURLADaCGH(newDir, application_url = "http://adacgh2.bioinfo.cnio.es"):
    """ Based on Pomelo II's Send_to_Pals.cgi."""
    f=open("idtype")
    idtype = f.read().strip()
    f.close()
    f=open("organism")
    organism = f.read().strip()
    f.close()
    if (idtype != "None" and organism != "None"):
        url_org_id = "org=" + organism + "&idtype=" + idtype + "&"
    else:
        url_org_id = ""
    gl_base = application_url + '/tmp/' + newDir + '/'
    gl1 = gl_base + 'Lost_for_PaLS.txt'
    gl2 = gl_base + 'Gained_for_PaLS.txt'
    gl3 = gl_base + 'Gained_or_Lost_for_PaLS.txt'
   
    outstr0 = '<br /> <hr> ' + \
              '<h3> Send results to <a href = "http://pals.bioinfo.cnio.es">' + \
              '<IMG BORDER="0" SRC="../../palsfavicon40.png" align="middle"></a></h3>'
   
    outstr = outstr0 + \
             '<p> Send set of <a href="http://pals.bioinfo.cnio.es?' + \
             url_org_id + 'datafile=' + gl1 + '"> genes with copy number LOSS to PaLS</a></p>' + \
             '<p> Send set of <a href="http://pals.bioinfo.cnio.es?' + \
             url_org_id + 'datafile=' + gl2 + '"> genes with copy number GAIN to PaLS</a></p>' + \
             '<p> Send set of <a href="http://pals.bioinfo.cnio.es?' + \
             url_org_id + 'datafile=' + gl3 + '"> genes with copy number ALTERATION (either gain or loss) to PaLS</a></p>'
    return(outstr)




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
#         compressedFile.add(rootname + '.' + str(fignum + 1) + '.jpeg',
#                            rootname + '.' + str(fignum + 1) + '.jpeg')
    os.chdir('/http/adacgh2/cgi')



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
    print '<p><a href="' + getBaseURL() + '?newDir=' + newDir + '">', 'http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html</a>.' 
    print '</p> </body> </html>'
    

## Output-generating functions
def printErrorRun():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    errormsg = open(tmpDir + "/error.msg").read()
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


def printOKRun():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>ADaCGH results</title></head><body>\n")

    if os.path.exists(tmpDir + "/ErrorFigure.png"):
        outf.write('<IMG BORDER="0" SRC="ErrorFigure.png">')
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
        outf.flush()

        methodUsed = open(tmpDir + '/methodaCGH').read()
        if (methodUsed == 'PSW') or (methodUsed == 'PSW\n'):
            arrayNames = open(tmpDir + '/arrayNames', mode = 'r').read().split('\n')[0].split('\t')
            outf.write('<h2>Island plots, gains <a href="http://adacgh2.bioinfo.cnio.es/help/adacgh-help.html#outputPSW">(help)</a></h2> \n')
            outf.write('<p>Click on thumbnails to expand.</p>')
            gains_fig_list = [''.join(['Gains.', aname]) for aname in arrayNames]
            thumb(tmpDir, gains_fig_list, outf, maxsthumb = 350)
            outf.write('<br />')
            
            outf.write('<h2>Island plots, losses <a href="http://adacgh2.bioinfo.cnio.es/help/adacgh-help.html#outputPSW">(help)</a></h2> \n')
            outf.write('<p>Click on thumbnails to expand.</p>')
            loss_fig_list = [''.join(['Losses.', aname]) for aname in arrayNames]
            thumb(tmpDir, loss_fig_list, outf, maxsthumb = 350)
            outf.write('<br />')

            outf.write('<p>Smith-Waterman results for all genes/clones are available from files ' +
                       '<a href="./Gains.Price.Smith.Waterman.results.txt">"Gains.Price.Smith.Waterman.results.txt"</a>' +
                       ' <a href="./Losses.Price.Smith.Waterman.results.txt">"Losses.Price.Smith.Waterman.results.txt."</a></p>')

            ##if os.path.exists(tmpDir + '/f1.R'): os.remove(tmpDir + '/f1.R')
            if os.path.exists(tmpDir + '/ace-figs.R'): os.remove(tmpDir + '/ace-figs.R')
            ##if os.path.exists(tmpDir + '/f1.Rout'): os.remove(tmpDir + '/f1.Rout')
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

            outf.write(printPalsURLADaCGH(newDir))

            outf.write("</body></html>")
            outf.close()
            Rresults.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")

        elif (methodUsed == 'ACE') or (methodUsed == 'ACE\n'):
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
            ##if os.path.exists(tmpDir + '/f1.R'): os.remove(tmpDir + '/f1.R')
            if os.path.exists(tmpDir + '/rerunACE.R'): os.remove(tmpDir + '/rerunACE.R')
            ##if os.path.exists(tmpDir + '/f1.Rout'): os.remove(tmpDir + '/f1.Rout')
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

            outf.write(printPalsURLADaCGH(newDir))

            outf.write("</body></html>")
            outf.close()
            Rresults.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
        else:
            results_print_general(outf, tmpDir, newDir, Rresults)



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
    
##redirectLoc = "/tmp/" + newDir
tmpDir = "/http/adacgh2/www/tmp/" + newDir

if not os.path.isdir(tmpDir):
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir is not a valid directory. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()
    


lam_check = open(tmpDir + '/lamCheckPID', mode = 'r'). readline().split()
lam_check_machine = lam_check[1]
lam_check_pid = lam_check[0]

## Were we already done in a previous execution?
## No need to reopen files or check anything else. Return url with results
## and bail out.
if os.path.exists(tmpDir + "/natural.death.pid.txt") or os.path.exists(tmpDir + "/killed.pid.txt"):
    print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'
    try:
        kill_lamcheck(lam_check_pid, lam_check_machine)
    except:
        None
    sys.exit()

## No, we were not done. Need to examine R output
   
Rrout = open(tmpDir + "/f1.Rout")
soFar = Rrout.read()
Rrout.close()
finishedOK = soFar.endswith("Normal termination\n")
errorRun = soFar.endswith("Execution halted\n")

if os.path.exists(tmpDir + '/RterminatedOK'):
    finishedOK = True

issue_echo('right after checking Rterminated', tmpDir)

if (not finishedOK) and (not errorRun) and (os.path.exists(tmpDir + "/pid.txt")):
    issue_echo('did we run out of time?', tmpDir)
    if (time.time() - os.path.getmtime(tmpDir + "/pid.txt")) > R_MAX_time:
	lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
        try:
            os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                      '; lamhalt -H; lamwipe -H')
        except:
            None
        kill_lamcheck(lam_check_pid, lam_check_machine)
 	printRKilled()
	os.rename(tmpDir + '/pid.txt', tmpDir + '/killed.pid.txt')
	try:
	    os.system("rm /http/adacgh2/www/R.running.procs/R." + newDir + "*")
	except:
	    None
	print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'
	sys.exit()

if errorRun > 0:
    issue_echo('errorRun is 1', tmpDir)
    kill_lamcheck(lam_check_pid, lam_check_machine)
    printErrorRun()
    os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    try:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    except:
        None
    try:
        lamkill = os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv + '; lamhalt -H; lamwipe -H')
    except:
        None
    try:
        os.remove("/http/adacgh2/www/R.running.procs/R." + newDir + "*")
    except:
	None
    print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'


elif finishedOK > 0:
    issue_echo('finishedOK is 1', tmpDir)
    kill_lamcheck(lam_check_pid, lam_check_machine)
    try:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    except:
        None
    try:
        lamkill = os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv + '; lamhalt -H; lamwipe -H')
    except:
        None
    try: os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    except: None
    printOKRun()
#     try: os.remove(tmpDir + '/f1.R')
#     except: None
    try:
        os.system("rm /http/adacgh2/www/R.running.procs/R." + newDir + "*")
    finally:
        print 'Location: http://adacgh2.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n' 
    
else:
    ## we only end up here if: we were not done in a previous run AND no process was overtime 
    ## AND we did not just finish. So we must continue.
    relaunchCGI()
    
issue_echo('END of checkdone.cgi', tmpDir)


