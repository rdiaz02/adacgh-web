#!/usr/bin/python
import sys
from toMapMod import *

nameMap = sys.argv[1]
idtype = sys.argv[2]
organism = sys.argv[3]


nameSrc = nameMap + '.png'
nameHTML = nameMap + '.html'

gene_F = open('geneNamesChr', mode = 'r')
gene_Names = [L.rstrip('\n') for L in gene_F]
gene_F.close()
map_F = open('pngCoordChr', mode = 'r')
map_coord = [L.rstrip('\n') for L in map_F]
map_F.close()




outList = []
outList.append(out_squeleton1)
outList.append(create_div(gene_Names))
outList.append(out_squeleton2)

outList.append(''.join(['<h1>Chromosome view: ', nameMap, '</h1>\n',
    '<img src="', nameSrc, '"usemap="#', nameMap, '" ISMAP>\n',
    '<map name="', nameMap, '">\n']))

if idtype == 'None' or organism == 'None':
    outList.append(create_map_none(gene_Names, map_coord, idtype, organism))
else:
    outList.append(create_map(gene_Names, map_coord, idtype, organism))

outList.append(['</map> </body> </html>'])


fileout = nameHTML
fout = open(fileout, mode = 'w')
for nl in range(len(outList)):
    fout.writelines(outList[nl])
fout.close()

