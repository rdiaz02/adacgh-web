#!/usr/bin/python

## All this code is copyright Ramon Diaz-Uriarte.
## Released under the Affero GPL.


## try to use add_to_MPIErrorLog and similar mechanisms for logging
## MPI attempts

import sys
import os
import cgi 
import time
import shutil
import glob
import random
import socket
##import fcntl


sys.path = sys.path + ['/asterias-web-apps/web-apps-common']
import counterApplications
from web_apps_config import *
import toMapMod

tmpDir = sys.argv[1]
ROOT_TMP_DIR = "/asterias-web-apps/adacgh2/www/tmp/"
newDir = tmpDir.replace(ROOT_TMP_DIR, "")
runningProcs = tmpDir.split('/tmp/')[0] + '/R.running.procs/'



def issue_echo(fecho, tmpDir):
    """Silly function to output small tracking files"""
    timeHuman = '##########   ' + \
                str(time.strftime('%d %b %Y %H:%M:%S')) 
    os.system('echo "' + timeHuman + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "' + fecho + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "    " >> ' + tmpDir + '/checkdone.echo')

issue_echo('at 00', tmpDir)    



## I think we no longer check tmpDir is OK, because this is not launched
## by the user, byt by the signsR.cgi file.




## procTable = tmpDir.split('/tmp/')[0] + '/R.running.procs/procTable'

## Must ensure the procTable exists and has a valid value
## No longer used
# if not os.path.exists(procTable):
#     fo = open(procTable, mode = 'w')
#     fo.write('0')
#     fo.close()

# def set_defaults_lam(tmpDir):
#     """ Set defaults for lamboot and Rslaves and number procs
#     based on the size of the data file. This is all heuristic,
#     but works for us with 6 GB RAM per node. The key idea is to
#     prevent swapping. ncpu are the Rslaves spawned by lamd, or the cpu=ncpu
#     in the lamb-host file. max_num_procs is the maximum number of simultaneous
#     adacgh processes running at any time.
#     We return the tuple ncpu, max_num_procs"""
#     datsize1 = 0
#     datsize2 = 0
#     if os.path.exists(tmpDir + '/acghData'):
#         datsize1 = int(os.popen('ls ' + tmpDir + '/acghData -sk').read().split()[0])
#     if os.path.exists(tmpDir + '/acghAndPosition'):
#         datsize2 = int(os.popen('ls ' + tmpDir + '/acghAndPosition -sk').read().split()[0])
#     datsize = max(datsize2, datsize1)
#     if datsize < 2000:
#         return (2, 3)
#     elif datsize < 6000:
#         return (2, 2)
#     elif datsize < 14000:
#         return (2, 1)
#     else:
#         return (1, 1) 


def Rrun(tmpDir, R_bin):
    """ Launch R."""
    issue_echo(' inside Rrun ', tmpDir)
    Rcommand = 'cd ' + tmpDir + \
               '; sleep 1; ' + \
               R_bin + ' --no-readline --no-save --slave <f1.R >>f1.Rout 2>> Status.msg &'
    issue_echo('the Rcommand is ' + Rcommand, tmpDir)
    Rtorun = os.system(Rcommand)



def collectZombies(k = 10):
    """ Make sure there are no zombies in the process tables.
    This is probably an overkill, but works.
    """
    for nk in range(k):
        try:
            tmp = os.waitpid(-1, os.WNOHANG)
        except:
            None

def kill_pid_machine(pid):
    'as it says: to kill somehting and not get error messages'
    fi, fo, fu = os.popen3('kill -s 9 ' + pid )
    fi.close()
    fo.close()
    fu.close()
            

# def kill_pid_machine(pid):
#     'as it says: to kill somehting somewhere'
#     os.system("kill -s 9 " + pid )

# def kill_pid_machine(pid, machine):
#     'as it says: to kill somehting somewhere'
#     os.system('ssh ' + machine + ' "kill -s 9 ' + pid + '"')


def results_print_general(outf, tmpDir, newDir, Rresults):
    outf.write('<h2>Segmented data plots <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#output">(help)</a></h2> \n')
    imagemaps2(tmpDir)
    thumb2(tmpDir, outf, maxsthumb = 350)
    # These are no longer created
    # thumb(tmpDir, ['All_arrays'], outf, maxsthumb = 350)
    outf.flush()
##    output_name = glob.glob(tmpDir + '/*.output.txt')[0].split('/')[-1]
    output_name = glob.glob(tmpDir + '/ADaCGH.results.txt')[0].split('/')[-1]
    outf.write('<p>Smoothed values for all genes/clones are available from file' +
               ' <a href="./' + output_name + '">"' + output_name + '".</a></p>')
    outf.write('<br />')
    # if os.path.exists(tmpDir + "/mcr.results.html"):
    #     outf.write('<h2>Minimal common regions</h2>\n')
    #     outf.write(open(tmpDir + "/mcr.results.html").read())
    outf.write('<br />')
    os.chdir(tmpDir)
    ll1 = glob.glob('*.log')
    for dname in ll1:
        os.remove(dname)
    outf.write('<hr> <a href="http://adacgh2.iib.uam.es/tmp/' +
               newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
    # try:
    #     outf.write(printPalsURLADaCGH(newDir))
    # except:
    #     None
         
    outf.write("</body></html>")
    outf.flush()
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
    os.system('tar -czf all.results.tar.gz * &')



# def printPalsURLADaCGH(newDir, application_url = "http://adacgh2.iib.uam.es"):
#     """ Based on Pomelo II's Send_to_Pals.cgi."""
#     f=open("idtype")
#     idtype = f.read().strip()
#     f.close()
#     f=open("organism")
#     organism = f.read().strip()
#     f.close()
#     if (idtype != "None" and organism != "None"):
#         url_org_id = "org=" + organism + "&idtype=" + idtype + "&"
#     else:
#         url_org_id = ""
#     gl_base = application_url + '/tmp/' + newDir + '/'
#     gl1 = gl_base + 'Lost_for_PaLS.txt'
#     gl2 = gl_base + 'Gained_for_PaLS.txt'
#     gl3 = gl_base + 'Gained_or_Lost_for_PaLS.txt'
   
#     outstr0 = '<br /> <hr> ' + \
#               '<h3> Send results to <a href = "http://pals.iib.uam.es">' + \
#               '<IMG BORDER="0" SRC="../../palsfavicon40.png" align="middle"></a></h3>'
   
#     outstr = outstr0 + \
#              '<p> Send set of <a href="http://pals.iib.uam.es?' + \
#              url_org_id + 'datafile=' + gl1 + '"> genes with copy number LOSS to PaLS</a></p>' + \
#              '<p> Send set of <a href="http://pals.iib.uam.es?' + \
#              url_org_id + 'datafile=' + gl2 + '"> genes with copy number GAIN to PaLS</a></p>' + \
#              '<p> Send set of <a href="http://pals.iib.uam.es?' + \
#              url_org_id + 'datafile=' + gl3 + '"> genes with copy number ALTERATION (either gain or loss) to PaLS</a></p>'
#     return(outstr)




# def pdf2html(rootname, tmpDir, outf, maxsthumb = 350):
#     """ From a multipage pdf obtain jpegs; the thumbnails are
#     inserted in the webpage, the large jpeg viewed on clik.
#     rootname is all the stuff before the `pdf',
#     tmpDir is the directory where the files live, maxsthumb
#     max size of thumbnail and outf the html file.
#     We also decrease the size of the jpeg for showing.
#     Finally, we add the generated jpeg to the compressed file"""

#     mst = str(maxsthumb)
#     mst2 = str(1600)
#     os.chdir(tmpDir)
#     os.system('/usr/bin/pdftoppm ' + rootname + '.pdf tmpppms')
#     tmps = glob.glob('tmpppms*.ppm')
#     for fignum in range(len(tmps)):
#         os.system('/usr/bin/ppmtojpeg ' + tmps[fignum] + ' > ' + rootname + '.'
#                   + str(fignum + 1) + '.jpeg')
#         os.system('/usr/bin/convert -size ' + mst + 'x' +
#                   mst + ' ' + rootname + '.' + str(fignum + 1) + '.jpeg' +
#                   ' -resize ' + mst + 'x' + mst + ' thumb.' + rootname + '.'
#                   + str(fignum + 1) + '.jpeg')
#         os.system('/usr/bin/convert ' + rootname + '.' + str(fignum + 1) + '.jpeg' +
#                   ' -resize ' + mst2 + 'x' + mst2 + ' ' + rootname + '.'
#                   + str(fignum + 1) + '.jpeg')
#         outf.write('<a href="' + rootname + '.' + str(fignum + 1) +
#                    '.jpeg"> <img src="' + 'thumb.' + rootname + '.'
#                   + str(fignum + 1) + '.jpeg"></a>\n')
#     os.chdir('/asterias-web-apps/adacgh2/cgi')

# FIXME: this is REALLY an ugly hack. Come on! Create a proper function and call it!
# def imagemaps(tmpDir):
#     os.chdir(tmpDir)
#     pngs = glob.glob("Chr*.png")
#     for pngfile in pngs:
#         os.system(''.join(['/asterias-web-apps/adacgh2/cgi/toMap.py ',
#                            os.path.splitext(pngfile)[0], " None", " None "]))

def imagemaps2(tmpDir):
    os.chdir(tmpDir)
    pngs = glob.glob("Chr*.png")
    for pngfile in pngs:
        toMapMod.html_image_maps(os.path.splitext(pngfile)[0])
        
        

        
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
    os.chdir('/asterias-web-apps/adacgh2/cgi')

def thumb2(tmpDir, outf, maxsthumb = 350):
    """ From a set of pngs, obtain thumbnails and
    add links to html. The thumbnails are
    inserted in the webpage, the large png viewed on clik.
    tmpDir is the directory where the files live,
    fnames is a list with the base file names to process
    maxsthumb   max size of thumbnail and outf the html file.
    """
    mst = str(maxsthumb)
    os.chdir(tmpDir)
    fnames = glob.glob("Chr*.png")
    for bname in fnames:
        bn2 = os.path.splitext(bname)[0]
        os.system(''.join(['/usr/bin/convert ', bname,
                           ' -resize ', mst, 'x', mst, ' thumb.', 
                           bn2, '.jpeg']))
        outf.write(''.join(['<a href="', bn2, '.html"> <img alt="',
	                    bn2, '" title="', bn2, '" src="thumb.',
                            bn2, '.jpeg"></a>']))
    os.chdir('/asterias-web-apps/adacgh2/cgi') ## FIXME
    

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
    issue_echo('inside relaunchCGI', tmpDir)

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
    issue_echo('end of relaunchCGI', tmpDir)
    


    

## Output-generating functions
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


def printOKRun():
    issue_echo('starting printOKRun', tmpDir)
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>ADaCGH results</title></head><body>\n")

    if os.path.exists(tmpDir + "/ErrorFigure.png"):
        outf.write('<IMG BORDER="0" SRC="ErrorFigure.png">')
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
        outf.flush()

        methodUsed = open(tmpDir + '/methodaCGH').read()
        # if (methodUsed == 'PSW') or (methodUsed == 'PSW\n'):
        #     arrayNames = open(tmpDir + '/arrayNames', mode = 'r').read().split('\n')[0].split('\t')
        #     outf.write('<h2>Island plots, gains <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputPSW">(help)</a></h2> \n')
        #     outf.write('<p>Click on thumbnails to expand.</p>')
        #     gains_fig_list = [''.join(['Gains.', aname]) for aname in arrayNames]
        #     thumb(tmpDir, gains_fig_list, outf, maxsthumb = 350)
        #     outf.write('<br />')
            
        #     outf.write('<h2>Island plots, losses <a href="http://adacgh2.iib.uam.es/help/adacgh-help.html#outputPSW">(help)</a></h2> \n')
        #     outf.write('<p>Click on thumbnails to expand.</p>')
        #     loss_fig_list = [''.join(['Losses.', aname]) for aname in arrayNames]
        #     thumb(tmpDir, loss_fig_list, outf, maxsthumb = 350)
        #     outf.write('<br />')

        #     outf.write('<p>Smith-Waterman results for all genes/clones are available from files ' +
        #                '<a href="./Gains.Price.Smith.Waterman.results.txt">"Gains.Price.Smith.Waterman.results.txt"</a>' +
        #                ' <a href="./Losses.Price.Smith.Waterman.results.txt">"Losses.Price.Smith.Waterman.results.txt."</a></p>')

        #     ##if os.path.exists(tmpDir + '/f1.R'): os.remove(tmpDir + '/f1.R')
        #     if os.path.exists(tmpDir + '/ace-figs.R'): os.remove(tmpDir + '/ace-figs.R')
        #     ##if os.path.exists(tmpDir + '/f1.Rout'): os.remove(tmpDir + '/f1.Rout')
        #     #if os.path.exists(tmpDir + '/.RData'): os.remove(tmpDir + '/.RData')
        #     ## allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
        #     os.chdir(tmpDir)
        #     ll1 = glob.glob('*.log')
        #     for dname in ll1:
        #         os.remove(dname)
        #     outf.write('<hr> <a href="http://adacgh2.iib.uam.es/tmp/' +
        #                newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  

        #     outf.write(printPalsURLADaCGH(newDir))

        #     outf.write("</body></html>")
        #     outf.close()
        #     Rresults.close()
        #     shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
        #     os.system('tar -czf all.results.tar.gz * &')
        # elif (methodUsed == 'ACE') or (methodUsed == 'ACE\n'):
        #     outf.write('<h2>FDR table</h2>')
        #     acefdrtable = open(tmpDir + "/ace.fdrtable.html")
        #     acefdr = acefdrtable.read()
        #     acefdrtable.close()
        #     outf.write(acefdr)
        #     outf.write('<form action="http://adacgh2.iib.uam.es/cgi-bin/ace.cgi" method="GET">\n')
        #     outf.write('<input type="hidden" NAME="newDir" VALUE="' + newDir + '">')
        #     currentfdr = str(open(tmpDir + '/aceFDR').readline())
        #     outf.write('<br><input type="TEXT" name="fdrace" value="' +
        #                currentfdr + '" size="10"  maxlength="10">\n')
        #     outf.write('<input value="Submit" type = "SUBMIT"> (Change the desired FDR and Press "Submit" to obtain figures with new FDR)')
        #     outf.write('<h2>Segmented plots</h2><p>Click on thumbnails to expand.</p>')
        #     thumb(tmpDir,  open(tmpDir + '/arrayNames', mode = 'r').read().split('\n')[0].split('\t'), outf, maxsthumb = 350)
        #     thumb(tmpDir, ['All_arrays'], outf, maxsthumb = 350)
        #     outf.write('<p>Inferred gains and losses available from file' +
        #                '<a href="./ACE.output.FDR=' + currentfdr + '.txt">' +
        #                '"ACE.output.FDR=' + currentfdr + '.txt"</a></p>')
        #     if os.path.exists(tmpDir + '/rerunACE.Rout'): os.remove(tmpDir + '/rerunACE.Rout')
        #     ##if os.path.exists(tmpDir + '/f1.R'): os.remove(tmpDir + '/f1.R')
        #     if os.path.exists(tmpDir + '/rerunACE.R'): os.remove(tmpDir + '/rerunACE.R')
        #     ##if os.path.exists(tmpDir + '/f1.Rout'): os.remove(tmpDir + '/f1.Rout')
        #     #if os.path.exists(tmpDir + '/.RData'): os.remove(tmpDir + '/.RData')
        #     ## allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
        #     os.chdir(tmpDir)
        #     ll1 = glob.glob('*.log')
        #     for dname in ll1:
        #         os.remove(dname)
        #     outf.write('<hr> <a href="http://adacgh2.iib.uam.es/tmp/' +
        #                newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  

        #     #outf.write(printPalsURLADaCGH(newDir))

        #     outf.write("</body></html>")
        #     outf.close()
        #     Rresults.close()
        #     shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
        #     os.system('tar -czf all.results.tar.gz * &')
        # else:
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


# def printMPIerror(tmpDir, numtries, application = 'ADaCGH2'):
#     if not os.path.exists('/asterias-web-apps/mpi.log/' + application + 'ErrorLog'):
#         os.system('touch /asterias-web-apps/mpi.log/' + application + 'ErrorLog')
#     outlog = open('/asterias-web-apps/mpi.log/' + application + 'ErrorLog', mode = 'a')
#     outlog.write('MPI fails more than ' + numtries + 'numtries on ' +
#                  time.ctime(time.time()) +
#                  ' Directory: ' + tmpDir + '\n')
#     outlog.close()
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


def printTooBusy(tmpDir, MAX_DURATION_TRY, application = 'ADaCGH'):
    if not os.path.exists('/asterias-web-apps/log/' + application + 'ErrorLog'):
        os.system('touch /asterias-web-apps/log/' + application + 'ErrorLog')
    outlog = open('/asterias-web-apps/log/' + application + 'ErrorLog', mode = 'a')
    outlog.write('Something fails on ' + time.ctime(time.time()) +
                 ' Directory: ' + tmpDir + '\n')
    outlog.close()
    out1 = open(tmpDir + "/natural.death.pid.txt", mode = "w")
    out2 = open(tmpDir + "/kill.pid.txt", mode = "w")
    out1.write('Cannot start!!')
    out2.write('Cannot start!!')
    out1.close()
    out2.close()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title> Cannot start application.</title></head><body>\n")
    outf.write("<h1> Cannot start application.</h1>")
    outf.write("<p> After " + str(MAX_DURATION_TRY) + " seconds we have been unable to ")
    outf.write(" start the application.</p>")
    outf.write("<p> Most likely this means the servers are too busy and many ")
    outf.write("are running ahead of yours. ")
    outf.write("Please try again later. You can also get in touch with us ")
    outf.write("if you think this is our error.</p>")
    outf.write("</body></html>")
    outf.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")



# def lamboot(lamSuffix, ncpu, runningProcs = runningProcs):
#     'Boot a lam universe and leave a sentinel file behind'
#     issue_echo('before sentinel inside lamboot', tmpDir)
#     issue_echo('newDir is ' + newDir, tmpDir)
#     issue_echo('lamSuffix ' + lamSuffix, tmpDir)
#     issue_echo('runningProcs ' + runningProcs, tmpDir)
#     sentinel = os.open(''.join([runningProcs, 'sentinel.lam.', newDir, '.', lamSuffix]),
#                        os.O_RDWR | os.O_CREAT | os.O_NDELAY)
#     issue_echo('before fullCommand inside lamboot', tmpDir)
#     fullCommand = 'export LAM_MPI_SESSION_SUFFIX="' + lamSuffix + \
#                   '"; /asterias-web-apps/mpi.log/tryBootLAM2.py ' + lamSuffix + \
#                   ' ' + str(ncpu)
#     issue_echo('before os.system inside lamboot', tmpDir)
#     lboot = os.system(fullCommand)
#     issue_echo('after lboot ---os.system--- inside lamboot. Exiting lamboot', tmpDir)


# def check_tping(lamSuffix, tmpDir, tsleep = 15, nc = 2):
#     """ Use tping to verify LAM universe OK.
#     tsleep is how long we wait before checking output of tping.
#     Verify also using 'lamexec C hostname' """
    
#     tmp2 = os.system('export LAM_MPI_SESSION_SUFFIX="' +\
#                      lamSuffix + '"; cd ' + tmpDir + \
#                      '; tping C N -c' + str(nc) + \
#                      ' > tping.out & ')
#     time.sleep(tsleep)
#     tmp = int(os.popen('cd ' + tmpDir + \
#                        '; wc tping.out').readline().split()[0])
#     os.system('rm ' + tmpDir + '/tping.out')
#     timeHuman = '##########   ' + \
#                 str(time.strftime('%d %b %Y %H:%M:%S')) 
#     os.system('echo "' + timeHuman + \
#               '" >> ' + tmpDir + '/checkTping.out')
#     if tmp == 0:
#         os.system('echo "tping fails" >> ' + \
#                   tmpDir + '/checkTping.out')
#         return 0
#     elif tmp > 0:
#         os.system('echo "tping OK" >> ' + \
#                   tmpDir + '/checkTping.out')
#         lamexec = os.system('export LAM_MPI_SESSION_SUFFIX="' +\
#                             lamSuffix + '"; lamexec C hostname')
#         if lamexec == 0:
#             os.system('echo "lamexec OK" >> ' + \
#                       tmpDir + '/checkTping.out')
#             return 1
#         else:
#             os.system('echo "lamexec fails" >> ' + \
#                       tmpDir + '/checkTping.out')
#             return 0
#     else:
#         os.system('echo "tping weird ' + str(tmp) + '" >> ' + \
#                   tmpDir + '/checkTping.out')
#         return 0



# def lam_crash_log(tmpDir, value):
#     """ Write to the lam crash log, 'recoverFromLAMCrash.out' """
#     timeHuman = str(time.strftime('%d %b %Y %H:%M:%S')) 
#     os.system('echo "' + value + '  at ' + timeHuman + \
#               '" >> ' + tmpDir + '/recoverFromLAMCrash.out')
    
def generic_crash_log(tmpDir, value):
    """ Write to the lam crash log, 'recoverFromLAMCrash.out' """
    timeHuman = str(time.strftime('%d %b %Y %H:%M:%S')) 
    os.system('echo "' + value + '  at ' + timeHuman + \
              '" >> ' + tmpDir + '/recoverFromLAMCrash.out')
    
# def recover_from_lam_crash(tmpDir, NCPU, MAX_NUM_PROCS, lamSuffix,
#                            runningProcs = runningProcs,
#                            machine_root = 'karl'):
#     """Check if lam crashed during R run. If it did, restart R
#     after possibly rebooting the lam universe.
#     Leave a trace of what happened."""
    
#     os.remove(''.join([runningProcs, 'sentinel.lam.', newDir, '.', lamSuffix]))
#     del_mpi_logs(tmpDir, machine_root)
#     lam_crash_log(tmpDir, 'Crashed')
#     ## We need to halt the universe, or else we can keep a lamd with no R hanging from
#         ## it, but that leads to too many lamds, so it cannot start. Like a vicious circle
#     lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
#     try:
#         os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
#                   '; lamhalt -H; lamwipe -H')
#     except:
#         None
#     issue_echo('inside recover_from_lam_crash: lamhalting', tmpDir)
#     try:
#         os.system('mv ' + tmpDir + '/mpiOK ' + tmpDir + '/previous_mpiOK')
#     except:
#         None
#     check_room = my_queue(MAX_NUM_PROCS)
#     if check_room == 'Failed':
#         printMPITooBusy(tmpDir, MAX_DURATION_TRY = 5 * 3600)

#     lam_ok = check_tping(lamSuffix, tmpDir)
#     if lam_ok == 0:
#         lboot = lamboot(lamSuffix, NCPU)
#     Rrun(tmpDir, lamSuffix)
#     lam_crash_log(tmpDir, '..... recovering')



def status_run(tmpDir):
    """ Read Status.msg and return status."""
    status_r = open(tmpDir + '/Status.msg').read()
    if status_r.find('Normal termination\n') > -1:
        return('FinishedOK')
    if status_r.find('Execution halted\n') > -1:
        return('Halted')
    if status_r.find('Running\n') > -1:
        return('Running')
    if status_r.find('Rmpi error\n') > -1:
        return('Error_mpi')
    if status_r.find('Run out of time; killed\n') > -1:
        return('Out_of_time')


# def did_R_crash_in_slaves(tmpDir, machine_root = 'karl'):
#     """ Verify whether R crashed in any of the slaves by
#     checking lam logs."""
#     R_LAM_MSGS = 'Error:  Error in'
#     lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
#     in_lam_logs = 0
#     for lam_log in lam_logs:
#         tmp1 = int(os.popen('grep "' + R_LAM_MSGS + '" ' + \
#                             lam_log + ' | wc').readline().split()[0])
#         if tmp1 > 0:
#             in_lam_logs = 1
#             break
#     if in_lam_logs > 0:
#         return True, lam_log
#     else:
#         return False, 'NA'



# def did_lam_crash(tmpDir, machine_root = 'karl'):
#     """ Verify whether LAM/MPI crashed by checking logs and f1.Rout
#     for single universe lamboot."""
#     OTHER_LAM_MSGS = 'Call stack within LAM:'
#     lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
#     in_error_msg = int(os.popen('grep MPI_Error_string ' + \
#                                 tmpDir + '/Status.msg | wc').readline().split()[0])
# #     no_universe = int(os.popen('grep "Running serial version of papply" ' + \
# #                                tmpDir + '/f1.Rout | wc').readline().split()[0])
# ## We do NOT want that, because sometimes a one node universe is legitimate!!!
#     if in_error_msg > 0:
#         for lam_log in lam_logs:
#             os.system('rm ' + lam_log)
# #     elif no_universe > 0:
# #         os.system("sed -i 's/Running serial version of papply/already_seen:running serial version of papply/g'" + \
# #                   tmpDir + "/f1.Rout")
#     else: ## look in lam logs
#         in_lam_logs = 0
#         for lam_log in lam_logs:
#             tmp1 = int(os.popen('grep "' + OTHER_LAM_MSGS + '" ' + \
#                                 lam_log + ' | wc').readline().split()[0])
#             if tmp1 > 0:
#                 in_lam_logs = 1
#                 break
#     if (in_error_msg > 0) or (in_lam_logs > 0):
#         return True
#     else:
#         return False
    
# def did_mpi_crash(tmpDir, machine_root = 'karl'):
#     """ Either Rmpi or LAM crashed"""
#     if (status_run(tmpDir) == 'Error_mpi') or \
#        did_lam_crash(tmpDir, machine_root):
#         return True
#     else:
#         return False

# def del_mpi_logs(tmpDir, machine_root = 'karl'):
#     """ Delete logs from LAM/MPI."""
#     lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
#     try:
#         os.system('rm ' + tmpDir + '/Status.msg')
#     except:
#         None
#     try:
#         for lam_log in lam_logs:
#             os.system('rm ' + lam_log)    
#     except:
#         None

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
                           

def cleanups(tmpDir, newDir, newnamepid,
             runningProcs = runningProcs,
             appl = 'adacgh2'):
    """ Clean up actions; kill lam, delete running.procs files, clean process table."""
##    lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    rinfo = open(tmpDir + '/current_R_proc_info', mode = 'r').readline().split()
    try:
        kill_pid_machine(rinfo[1])
    except:
        None
    # try:
    #     os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
    #               '; lamhalt -H; lamwipe -H')
    # except:
    #     None
    try:
        os.system('rm /asterias-web-apps/' + appl + '/www/R.running.procs/R.' + newDir + '*')
    except:
        None
    try:
        os.rename(tmpDir + '/pid.txt', tmpDir + '/' + newnamepid)
    except:
        None
    # try:
    #     os.remove(''.join([runningProcs, 'sentinel.lam.', newDir, '.', lamSuffix]))
    # except:
    #     None


def finished_ok(tmpDir):
    """ check ok termination and act accordingly."""
    if status_run(tmpDir) == 'FinishedOK':
        return True
    else:
        return False

def halted(tmpDir):
    """ check halted execution and act accordingly."""
    if status_run(tmpDir) == 'Halted':
        return True
    else:
        return False


def master_out_of_time(time_start):
    """If this process run longer than allowed, kill it and kill lam and R."""
    if (time.time () - time_start) > R_MAX_time:
        return True
    else:
        return False
        

# def add_to_proc_table(max_num_procs, add_procs = 1):
#     """Try to add add_procs to the process table. If it can
#     returns OK, otherwise (e.g., too many procs) return Failed.
#     Locking would be great ... but it does not work over NFS. """
    
#     fo = open(procTable, mode = 'r+')
#     fcntl.flock(fo.fileno(), fcntl.LOCK_EX)
#     currentProcs = int(fo.read())
#     if currentProcs >= max_num_procs:
#         fcntl.flock(fo.fileno(), fcntl.LOCK_UN)
#         fo.close()
#         return 'Failed'
#     else:
#         fo.seek(0)
#         fo.write(str(currentProcs + add_procs))
#         fcntl.flock(fo.fileno(), fcntl.LOCK_UN)
#         fo.close()
#         return 'OK'

# def add_to_proc_table(max_num_procs, add_procs = 1):
#     """Try to add add_procs to the process table. If it can
#     returns OK, otherwise (e.g., too many procs) return Failed."""
#     fo = open(procTable, mode = 'r+')
#     currentProcs = int(fo.read())
#     if currentProcs >= max_num_procs:
#         fo.close()
#         return 'Failed'
#     else:
#         fo.seek(0)
#         fo.write(str(currentProcs + add_procs))
#         fo.close()
#         return 'OK'




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
             MAX_DURATION_TRY = MAX_DURATION_TRY_ADaCGH):
    """ Wait here until the number of processes is smaller than
    MAX_NUM_PROCS 
    But only wait for MAX_DURATION. Check
    every CHECK_QUEUE seconds. If able to find an opening, return
    OK, otherwise return Failed"""

    # taken from the smae function in signs
    out_value = 'OK'
    startTime = time.time()
    while True:
#        killedlamandr = os.system('/asterias-web-apps/mpi.log/killOldLamAllMachines.py')
        issue_echo('     inside my_queue ', tmpDir)
        if (time.time() - startTime) > MAX_DURATION_TRY:
            out_value = 'Failed'
            break
 #       num_lamd = int(os.popen('pgrep -u www lamd | wc').readline().split()[0])
        num_sentinel = int(len(glob.glob(''.join([runningProcs, 'sentinel.lam.*']))))
        if (num_sentinel < MAX_NUM_PROCS):
            issue_echo(' num_sentinel = ' + str(num_sentinel), tmpDir)
            break
        else:
	    issue_echo('     wait:  ' + \
                       '; num_sentinel = ' + str(num_sentinel), tmpDir)
            time.sleep(CHECK_QUEUE + random.uniform(0.1, 5))
    return out_value


    
# def my_queue(MAX_NUM_PROCS,
#              runningProcs = runningProcs,
#              ADD_PROCS = 1,
#              CHECK_QUEUE = 23,
#              MAX_DURATION_TRY = 25 * 3600):
#     """ Wait here until the number of processes is smaller than
#     MAX_NUM_PROCS and number of slaves smaller than MAX_NUM_PROCS + ADD_PROCS
#     (so we allow for other apps. launching lamd).
#     But only wait for MAX_DURATION. Check
#     every CHECK_QUEUE seconds. If able to find an opening, return
#     OK, otherwise return Failed"""
#     out_value = 'OK'
#     startTime = time.time()
#     while True:
#         killedlamandr = os.system('/asterias-web-apps/mpi.log/killOldLamAllMachines.py')
#         issue_echo('     inside my_queue ', tmpDir)
#         if (time.time() - startTime) > MAX_DURATION_TRY:
#             out_value = 'Failed'
#             break
#         num_lamd = int(os.popen('pgrep -u www-data lamd | wc').readline().split()[0])
#         num_sentinel = int(len(glob.glob(''.join([runningProcs, 'sentinel.lam.*']))))
#         if (num_lamd < (MAX_NUM_PROCS + ADD_PROCS)) and (num_sentinel < MAX_NUM_PROCS):
#             issue_echo('     OK; num_lamd = ' + str(num_lamd) + \
#                        '; num_sentinel = ' + str(num_sentinel), tmpDir)
#             break
#         else:
# 	    issue_echo('     wait:  num_lamd = ' + str(num_lamd) + \
#                        '; num_sentinel = ' + str(num_sentinel), tmpDir)
#             time.sleep(CHECK_QUEUE + random.uniform(0.1, 5))
#     return out_value

# def generate_lam_suffix(tmpDir):
#     """As it says. Generate and write it out"""
#     lamSuffix = str(int(time.time())) + \
#                 str(os.getpid()) + str(random.randint(10, 999999))
#     lamenvfile = open(tmpDir + '/lamSuffix', mode = 'w')
#     lamenvfile.write(lamSuffix)
#     lamenvfile.flush()
#     lamenvfile.close()
#     return lamSuffix



## Starting. First, the very first run.

issue_echo('starting', tmpDir)

        
# NCPU, MAX_NUM_PROCS = set_defaults_lam(tmpDir)

try:
    counterApplications.add_to_counter_log('ADaCGH2', tmpDir, socket.gethostname())
except:
    None

issue_echo('at 2', tmpDir)

# lamSuffix = generate_lam_suffix(tmpDir)

issue_echo('at 3', tmpDir)

time.sleep(random.uniform(0.1, 15)) ## Break ties if starting at identical times

check_room = my_queue(MAX_adacgh, MAX_DURATION_TRY = MAX_DURATION_TRY_ADaCGH)
issue_echo('after check_room', tmpDir)

if check_room == 'Failed':
    printMPITooBusy(tmpDir, MAX_DURATION_TRY = 5 * 3600)
    sys.exit()

# issue_echo('before lamboot', tmpDir)
# lamboot(lamSuffix, NCPU)
# issue_echo('after lamboot', tmpDir)

# counterApplications.add_to_LAM_SUFFIX_LOG(lamSuffix, 'ADaCGH2', tmpDir,
#                                           socket.gethostname())

issue_echo('before  Rrun', tmpDir)

Rrun(tmpDir, R_bin)

issue_echo('after Rrun', tmpDir)

time_start = time.time()
time.sleep(TIME_BETWEEN_CHECKS + random.uniform(0.1, 3))

count_mpi_crash = 0

while True:
    if did_run_out_of_time(tmpDir, R_MAX_time):
        issue_echo('run out of time', tmpDir)
        cleanups(tmpDir, newDir, 'killed.pid.txt')
        printRKilled()
        break
    elif finished_ok(tmpDir):
        issue_echo('finished OK', tmpDir)
        cleanups(tmpDir, newDir, 'natural.death.pid.txt')
        printOKRun()
        break
    elif halted(tmpDir):
        issue_echo('halted', tmpDir)
        cleanups(tmpDir, newDir, 'natural.death.pid.txt')
        printErrorRun(tmpDir + '/Status.msg')
        break
    # elif did_R_crash_in_slaves(tmpDir, machine_root = 'karl')[0]:
    #     issue_echo('R crash in slaves', tmpDir)
    #     cleanups(tmpDir, newDir, 'natural.death.pid.txt', lamSuffix)
    #     printErrorRun(did_R_crash_in_slaves(tmpDir, machine_root = 'karl')[1])
    #     break
    elif master_out_of_time(time_start):
        issue_echo('master out of time', tmpDir)
        cleanups(tmpDir, newDir, 'killed.pid.txt')
        printRKilled()
        break
    # elif did_mpi_crash(tmpDir, machine_root = 'karl'):
    #     count_mpi_crash += 1
    #     counterApplications.add_to_MPIErrorLog('ADaCGH2',
    #                                            tmpDir, socket.gethostname(),
    #                                            message = 'MPI crash')
    #     if count_mpi_crash > MAX_MPI_CRASHES:
    #         printMPIerror(tmpDir, MAX_MPI_CRASHES)
    #         cleanups(tmpDir, newDir, 'MPIerror.pid.txt', lamSuffix)
    #         break
    #     else:
    #         recover_from_lam_crash(tmpDir, NCPU, MAX_NUM_PROCS,
    #                                lamSuffix,
    #                                machine_root = 'karl')
    else:
        generic_crash_log(tmpDir, 'NoCrash') ## if we get here, this much we know
    time.sleep(TIME_BETWEEN_CHECKS)



issue_echo('at the very end!', tmpDir)

