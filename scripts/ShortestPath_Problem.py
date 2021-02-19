import collections
import xml.etree.ElementTree as ET
import networkx as nx
import matplotlib.pyplot as plt
import numpy
from scipy import linalg
import math


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



def permutation(liste):
    # Gibt Permutationen der Listenelemente aus, ohne (name, name)
    tupelListe = []
    for name in liste:
        for name2 in liste:
            if name != name2:
                tupelListe.append((name, name2))
    return tupelListe


def permutationSingle(liste):
    # Wie Permutation, aber falls (name, name2) schon existiert wird (name2, name) nicht noch hinzugefügt
    tupelListe = []
    for name in liste:
        for name2 in liste:
            if name != name2:
                if (name2, name) not in tupelListe:
                    tupelListe.append((name, name2))
    return tupelListe


def tempus(time, person):
    # Liste Person + Zeit(int)
    #
    newlist = []
    for name in person:
        name2 = name + "~@" + str(time)
        newlist.append(name2)
    return newlist


def checkForDoubles(String):
    # Prüft ob mehrere Personen gleichzeitig sprechen;
    # Rückgabe Liste der sprechenden Pers. : "#Anna #Berta" -> ["#Anna", "#Berta"]
    stringList = []
    y = len(String)
    x = 0
    for i in range(len(String) - 1, -1, -1):
        if String[i] == "#":
            x = i
            stringList.append(String[x:y])
            y = x - 1
    return stringList


def adjustColor(graph):
    # Nutzt adjustValue um eine Liste der Farben der Knoten zurück zu geben
    degreeTuple = []
    for node in list(graph):
        degreeTuple.append((node, nx.degree(graph, node)))

    colorTime = []
    referenceColorMax = (255 / 256, 50 / 256, 00 / 256)
    referenceColorMean = (255 / 256, 160 / 256, 00 / 256)
    referenceColorMin = (230 / 256, 220 / 256, 30 / 256)
    colors = [referenceColorMin, referenceColorMean, referenceColorMax]

    degreeMax = 0
    degreeMean = 0
    degreeMin = 4000
    for node, degree in degreeTuple:
        if degreeMax < degree:
            degreeMax = degree
        if degreeMin > degree:
            degreeMin = degree
        degreeMean += degree
    degreeMean = round(degreeMean / len(degreeTuple))
    degreeValues = [degreeMin, degreeMean, degreeMax]
    for node, dvalue in degreeTuple:
        colorTime.append(adjustValue(dvalue, degreeValues, colors))
    return colorTime


def adjustSize(graph):
    # Nutzt adjustValue um eine Liste mit Größen der Knoten zurück zu geben
    degreeTuple = []
    for node in list(graph):
        degreeTuple.append((node, nx.degree(graph, node)))

    sizeTime = []
    sMax = 300,
    sMean = 150,
    sMin = 50,
    # ',' am Ende macht das Objekt zu einem 1-Tupel, so kann in adjustValue per Index drauf zugegriffen werden
    sizes = [sMin, sMean, sMax]
    degreeMax = 0
    degreeMean = 0
    degreeMin = 4000
    for node, degree in degreeTuple:
        if degreeMax < degree:
            degreeMax = degree
        if degreeMin > degree:
            degreeMin = degree
        degreeMean += degree
    degreeMean = round(degreeMean / len(degreeTuple))
    degreeValues = [degreeMin, degreeMean, degreeMax]
    for node, dvalue in degreeTuple:
        sizeTime.append(adjustValue(dvalue, degreeValues, sizes))
    return sizeTime


def adjustValue(dcurrent, degreevalues, values):
    # dcurrent = degree of current node
    # degreevalues = List of degreeMin, degreeMean, degreeMax
    # values = List of valueMin, valueMean, valueMax
    adjustedValues = []
    if dcurrent == degreevalues[1]:
        return values[1]
    elif dcurrent > degreevalues[1]:
        for i in range(len(values[0])):
            adjustedValues.append((values[1][i] + (values[2][i] - values[1][i]) * (
                    dcurrent - degreevalues[1]) / (
                                           degreevalues[2] - degreevalues[1])) / 1)
    else:
        for i in range(len(values[0])):
            adjustedValues.append((values[1][i] - (values[1][i] - values[0][i]) * (
                    dcurrent - degreevalues[0] + 1) / (
                                           degreevalues[1] - degreevalues[0] + 1)) / 1)
    return tuple(adjustedValues)


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