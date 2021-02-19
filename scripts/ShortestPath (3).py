import collections
import xml.etree.ElementTree as ET
import networkx as nx
import matplotlib.pyplot as plt
import numpy
from scipy import linalg
import math

from GermanDrama_new.py import permutation
#from GermanDrama_new.py import permutationSingle
from GermanDrama_new.py import tempus
from GermanDrama_new.py import checkForDoubles
from GermanDrama_new.py import adjustColor
from GermanDrama_new.py import adjustSize
#from GermanDrama_new.py import adjustValue

def FloydWahrschall(graph):
    V = len(graph)
    shortestPath = graph
    for k in range(V):
        shortestPathI = [list(row) for row in shortestPath]
        for l in range(V):
            for j in range(V):
                shortestPathI[j][l] = min(shortestPath[j][l], shortestPath[l][k] + shortestPath[k][j])
        shortestPath = shortestPathI
    return shortestPath


document = 'input\\beer-der-paria.xml'
#   'input\\arnim-das-loch.xml' ##VOLLSTÄNDIGER GRAPH!!
#   'input\\alberti-brot.xml'
#   'input\\beer-der-paria.xml'
#   'input\\test.xml'

root = ET.parse(document).getroot()
NS = {'tei': 'http://www.tei-c.org/ns/1.0'}

# finde pro Szene vorkommende Personen und speichere diese als Liste in einer Liste
sceneList = []
for div_tag in root.findall(".//tei:div", namespaces=NS):
    newlist = []
    for sp_tag in div_tag.findall(".//tei:sp[@who]", namespaces=NS):
        for name in checkForDoubles(sp_tag.get('who')):
            if name not in newlist:
                newlist.append(name)
    if (len(newlist) > 0):
        sceneList.append(newlist)

# finde insgesamt alle vorkommenden Personen
speakerList = []
for sp_tag in root.findall(".//tei:sp[@who]", namespaces=NS):
    value = sp_tag.get('who')
    for name in checkForDoubles(value):
        if name not in speakerList and value is not None:
            speakerList.append(name)

dynamicGraph = nx.DiGraph()
for i in range(len(sceneList)):
    dynamicGraph.add_nodes_from(tempus(i, speakerList))
    tempusI = tempus(i, sceneList[i])
    dynamicGraph.add_edges_from(permutation(tempusI))
for i in range(len(sceneList) - 1):
    for j in speakerList:
        dynamicGraph.add_edge(j + "~@" + str(i), j + "~@" + str(i + 1))

rgb1 = 177
rgb2 = 224
rgb3 = 185

plt.figure(1, dpi=400)
dynamicGraph_Layout = nx.spring_layout(dynamicGraph)
nx.draw_networkx(dynamicGraph, node_size=adjustSize(dynamicGraph), node_color=adjustColor(dynamicGraph),
                 width=0.2, with_labels=False, pos=dynamicGraph_Layout)
nx.draw_networkx_labels(dynamicGraph, pos=dynamicGraph_Layout)
plt.title("aggregate network of " + document.replace("input\\", ""))
#plt.show()

# Matrizen
Matrix = nx.adjacency_matrix(dynamicGraph)
print(Matrix)

adjMatrix1 = nx.to_numpy_matrix(dynamicGraph)
# print(adjMatrix1)

adjMatrix2 = nx.adjacency_matrix(dynamicGraph).todense()
print(adjMatrix2)

# G = (V, E)
# V = Knoten; V = 1, 2, ... , N Gewichte w[i,j]
# E = Kanten
# kürzester Pfad -> shortestPath(i,j,k)
# shortestPath(i,j,0) -> Gewichtsmatrix
# shortestPath(i,j,k+1) = min(shortestPath(i,j,k), shortestPath(i,k+1,k)+shortestPath(k+1,j,k))
INF = 9999

#print(FloydWahrschall(adjMatrix1))

Graph = [[0, 5, INF, 10],
             [INF, 0, 3, INF],
             [INF, INF, 0, 1],
             [INF, INF, INF, 0]
        ]

adjMatrixListTest = adjMatrix1.tolist()
for x in range(len(adjMatrixListTest)):
    for y in range(len(adjMatrixListTest)):
        if adjMatrixListTest[x][y]==0:
            adjMatrixListTest[x][y] = INF
#adjMatrix2 = numpy.array(adjMatrix2)[0].tolist()
print(FloydWahrschall(adjMatrixListTest))
#print(PrintFloydWahrschall(FloydWahrschall(Graph)))