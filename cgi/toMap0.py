#!/usr/bin/python
import sys


nameMap = sys.argv[1]
idtype = sys.argv[2]
organism = sys.argv[3]

#nameMap = 'Chr1@yetAnother@88'
#idtype = 'acc'
#organism = 'Hs'

nameSrc = nameMap + '.png'
nameHTML = nameMap + '.html'

gene_F = open('geneNamesChr', mode = 'r')
gene_Names = [L.rstrip('\n') for L in gene_F]
## gene_Names = gene_F.read().split('\n') this leaves a trailing 
## empty line
gene_F.close()
map_F = open('pngCoordChr', mode = 'r')
map_coord = [L.rstrip('\n') for L in map_F]
map_F.close()


def create_map():
    global outList
    for numline in range(len(gene_Names)):
        gene_line = gene_Names[numline]
        coords = map_coord[numline].split('\t')
        outstring = ''.join(['<area shape="circle" coords="',
                            coords[0], ' ', coords[1], ' ', 
                            coords[2], 
                            '" title="',gene_line,
                            '" onClick="fixedtooltip(',str(numline + 1),
                            ', \'<a class="tip" href="javascript:hidetip(1)">X</a><a class=\'tip2\' href="http://idclight.iib.uam.es/idclight.prog?',
                            gene_line,'&idtype=', idtype, '&org=', organism, '">',
                            gene_line, '</a>\', this, event, \'90px\'); return false"  >   '])
        outList.append(outstring)

def create_map_none():
    global outList
    """Like create_map, but when there are no known identifiers.
    It changes the javascript call"""
    for numline in range(len(gene_Names)):
        gene_line = gene_Names[numline]
        coords = map_coord[numline].split('\t')
        outstring = ''.join(['<area shape="circle" coords="',
                            coords[0], ' ', coords[1], ' ', 
                            coords[2], 
                            '" title="',gene_line,
                            '" onClick="fixedtooltip(',str(numline + 1),
                            ', \'<a class="tip" href="javascript:hidetip(1)">X</a> <a class=\'tip2\'',
                            gene_line, '</a>\', this, event, \'90px\'); return false"  >   '])
        outList.append(outstring)


def create_div():
    global outList
    for numline in range(len(gene_Names)):
        outstring = ''.join(['document.write(\'<div id="', str(numline + 1),
                            '" class="fixedtipdiv" style="visibility:hidden;width:\'+tipwidth+\';background-color:\'+tipbgcolor+\'" ></div>\')\n'])
        outList.append(outstring)



out_squeleton1 = """
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
<head>
<title>Chromosome view</title>
<style type="text/css">

.fixedtipdiv{
position:absolute;
padding: 0px;
border:1px solid blue;
font-family: Verdana,Arial,Helvetica;
font-size: x-small;
font-style: normal;
color:#000000;
line-height:18px; vertical-align:text-top;
z-index:100;
a:link{color: black;}
}

a.tip:link{color: black; text-align:right; vertical-align:top; display:block; padding-bottom:0px; padding-top:0px; font-size: small; text-decoration: none; font-weight: bolder; background-color:#FFFFFF;}
a.tip2:link{color: black; padding-left:4px; padding-bottom:3px;text-decoration: none;}
a.tip2:active{color: black; padding-left:4px;padding-bottom:3px;text-decoration: underline;font-weight: bolder;}
a.tip2:visited{color: red; padding-left:4px;padding-bottom:3px;text-decoration: none;}



</style>

</head>


<body>
<script type="text/javascript">

/**********************************************
*  Original code taken from ToolTip script, by
*  Dynamic Drive (see below). A few changes, additions, and
*  deletions by us (Oscar Rueda Palacio and
*  Ramon Diaz-Uriarte).  A few ideas taken from 
*  overLIB.4.21, by Erik Bosrup.
*
*  These notices apply to this file (js.squeleton1),
*  as well as js.squeleton2).
**********************************************/


/***********************************************
* Fixed ToolTip script- © Dynamic Drive (www.dynamicdrive.com)
* This notice MUST stay intact for legal use
* Visit http://www.dynamicdrive.com/ for full source code
***********************************************/
        
var tipwidth='150px' //default tooltip width
var tipbgcolor="#CCCCFF"  //tooltip bgcolor
var disappeardelay=250  //tooltip disappear speed onMouseout (in miliseconds)
var vertical_offset="0px" //horizontal offset of tooltip from anchor link
var horizontal_offset="-1px" //horizontal offset of tooltip from anchor link

/////No further editting needed

var ie4=document.all
var ns6=document.getElementById&&!document.all

if (ie4||ns6)
"""

out_squeleton2 = """
function getposOffset(what, offsettype){
var totaloffset=(offsettype=="left")? what.offsetLeft : what.offsetTop;
var parentEl=what.offsetParent;
while (parentEl!=null){
totaloffset=(offsettype=="left")? totaloffset+parentEl.offsetLeft : totaloffset+parentEl.offsetTop;
parentEl=parentEl.offsetParent;
}
return totaloffset;
}


function showhide(obj, e, visible, hidden, tipwidth){
if (ie4||ns6)
dropmenuobj.style.left=dropmenuobj.style.top=-500
if (tipwidth!=""){
dropmenuobj.widthobj=dropmenuobj.style
dropmenuobj.widthobj.width=tipwidth
}
if (e.type=="click" && obj.visibility==hidden || e.type=="mouseover")
obj.visibility=visible
else if (e.type=="click")
obj.visibility=visible
}

function iecompattest(){
return (document.compatMode && document.compatMode!="BackCompat")? document.documentElement : document.body
}

function clearbrowseredge(obj, whichedge){
var edgeoffset=(whichedge=="rightedge")? parseInt(horizontal_offset)*-1 : parseInt(vertical_offset)*-1
if (whichedge=="rightedge"){
var windowedge=ie4 && !window.opera? iecompattest().scrollLeft+iecompattest().clientWidth-15 : window.pageXOffset+window.innerWidth-15
dropmenuobj.contentmeasure=dropmenuobj.offsetWidth
if (windowedge-dropmenuobj.x < dropmenuobj.contentmeasure)
edgeoffset=dropmenuobj.contentmeasure-obj.offsetWidth
}
else{
var windowedge=ie4 && !window.opera? iecompattest().scrollTop+iecompattest().clientHeight-15 : window.pageYOffset+window.innerHeight-18
dropmenuobj.contentmeasure=dropmenuobj.offsetHeight
if (windowedge-dropmenuobj.y < dropmenuobj.contentmeasure)
edgeoffset=dropmenuobj.contentmeasure+obj.offsetHeight
}
return edgeoffset
}

function fixedtooltip(index, menucontents, obj, e, tipwidth){
if (window.event) event.cancelBubble=true
else if (e.stopPropagation) e.stopPropagation()
clearhidetip()
dropmenuobj=document.getElementById? document.getElementById(index) : fixedtipdiv
dropmenuobj.innerHTML=menucontents

if (ie4||ns6){
showhide(dropmenuobj.style, e, "visible", "hidden", tipwidth)
dropmenuobj.x=e.pageX
dropmenuobj.y=e.pageY
dropmenuobj.style.left=dropmenuobj.x-clearbrowseredge(obj, "rightedge")+"px"
dropmenuobj.style.top=dropmenuobj.y-clearbrowseredge(obj, "bottomedge")
}
}

function hidetip(index){
if (typeof dropmenuobj!="undefined"){
if (ie4||ns6)
dropmenuobj=document.getElementById? document.getElementById(index) : fixedtipdiv
dropmenuobj.style.visibility="hidden"
}
}

function delayhidetip(){
if (ie4||ns6)
delayhide=setTimeout("hidetip()",disappeardelay)
}

function clearhidetip(){
if (typeof delayhide!="undefined")
clearTimeout(delayhide)
}

</script>
"""


outList = []
outList.append(out_squeleton1)
create_div()
outList.append(out_squeleton2)

outList.append(''.join(['<h1>Chromosome view:', nameMap, '</h1>\n',
    '<img src="', nameSrc, '"usemap="#', nameMap, '" ISMAP>\n',
    '<map name="', nameMap, '">']))

if idtype == 'None' or organism == 'None':
    create_map_none()
else:
    create_map()

outList.append(['</map> </body> </html>'])


fileout = nameHTML
fout = open(fileout, mode = 'w')
for nl in range(len(outList)):
    fout.writelines(outList[nl])
fout.close()

