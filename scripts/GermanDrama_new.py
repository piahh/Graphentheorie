import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import networkx as nx


def permutation(liste):
    #Gibt Permutationen der Listenelemente aus, ohne (name, name)
    tupelListe = []
    for name in liste:
        for name2 in liste:
            if name != name2:
                tupelListe.append((name, name2))
    return tupelListe


def permutationSingle(liste):
    #Wie Permutation, aber falls (name, name2) schon existiert wird (name2, name) nicht noch hinzugefügt
    tupelListe = []
    for name in liste:
        for name2 in liste:
            if name != name2:
                if (name2, name) not in tupelListe:
                    tupelListe.append((name, name2))
    return tupelListe


def tempus(time, person):
    #Liste Person + Zeit(int)
    #
    newlist=[]
    for name in person:
        name2 = name+"~@"+str(time)
        newlist.append(name2)
    return newlist


def checkForDoubles(String):
    # Prüft ob mehrere Personen gleichzeitig sprechen;
    # Rückgabe Liste der sprechenden Pers. : "#Anna #Berta" -> ["#Anna", "#Berta"]
    stringList = []
    y = len(String)
    x = 0
    for i in range(len(String)-1, -1, -1):
        if String[i] == "#":
            x = i
            stringList.append(String[x:y])
            y = x-1
    return stringList


def adjustColor(graph):
    #Nutzt adjustValue um eine Liste der Farben der Knoten zurück zu geben
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
    #Nutzt adjustValue um eine Liste mit Größen der Knoten zurück zu geben
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
            adjustedValues.append(( values[1][i] + (values[2][i] - values[1][i]) * (
                dcurrent - degreevalues[1]) / (
                    degreevalues[2] - degreevalues[1])) /1)
    else:
        for i in range(len(values[0])):
            adjustedValues.append(( values[1][i] - (values[1][i] - values[0][i]) * (
                dcurrent - degreevalues[0] + 1) / (
                    degreevalues[1] - degreevalues[0] + 1)) /1)
    return tuple(adjustedValues)


def FreemanIndexDegreeCentrality(graph):
    #Berechnet den FreemanIndex über den Grad des Knoten mit Hilfe der Degree Centrality
    degree = dict(nx.degree_centrality(graph))
    degreeMax = 0
    sum = 0
    for value in degree.values():
        if degreeMax < value:
            degreeMax = value
    for nodes in graph:
        sum +=(degreeMax - degree[nodes])
    return sum/((len(degree)-1)*(len(degree)-2))

def FreemanIndexEigenvectorCentrality(graph):
    #Berechnet den FreemanIndex über den Grad des Knoten mit Hilfe der Eigenvector Centrality
    degree = dict(nx.eigenvector_centrality(graph))
    degreeMax = 0
    sum = 0
    for value in degree.values():
        if degreeMax < value:
            degreeMax = value
    for nodes in graph:
        sum +=(degreeMax - degree[nodes])
    return sum/((len(degree)-1)*(len(degree)-2))

def FreemanVitalityDegreeCentrality(graph, speakersList, numberOfScenes):
    vitalityDic = {}
    if speakersList[0] in graph:
        for node in speakersList:
            neighborss = list(graph.neighbors(node))
            adjList = []
            for neighbor in neighborss:
                adjList.append((node, neighbor))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexDegreeCentrality(graph)
            graph.add_edges_from(adjList)

    else:
        for node in speakersList:
            instanceOfNode = []
            adjList = []
            for timeScene in range(numberOfScenes):
                instanceOfNode.append(node + "~@" + str(i))
            for timeNode in instanceOfNode:
                neighborss = list(graph.neighbors(timeNode))
                inEdges = list(graph.in_edges(timeNode))
                for neighbor in neighborss:
                    adjList.append((timeNode, neighbor))
                for inEdge in inEdges:
                    adjList.append((inEdge[0], timeNode))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexDegreeCentrality(graph)
            graph.add_edges_from(adjList)
    return vitalityDic

def FreemanVitalityMultiGraph(graph, speakersList, numberOfScenes):
    #FreemanVitality für MultiGraph mit FreemandIndex mit Hilfe der Degree Centrality
    vitalityDic = {}
    if speakersList[0] in graph:
        for node in speakersList:
            neighborss = list(graph.neighbors(node))
            adjList = []
            for neighbor in neighborss:
                adjList.append((node, neighbor))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexDegreeCentrality(graph)
            graph.add_edges_from(adjList)

    else:
        for node in speakersList:
            instanceOfNode = []
            adjList = []
            for timeScene in range(numberOfScenes):
                instanceOfNode.append(node + "~@" + str(i))
            for timeNode in instanceOfNode:
                neighborss = list(graph.neighbors(timeNode))
                #inEdges = list(graph.in_edges(timeNode))
                for neighbor in neighborss:
                    adjList.append((timeNode, neighbor))
                #for inEdge in inEdges:
                    #adjList.append((inEdge[0], timeNode))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexDegreeCentrality(graph)
            graph.add_edges_from(adjList)
    return vitalityDic

def FreemanVitalityMultiGraphEigenvectorCentrality(graph, speakersList, numberOfScenes):
    # FreemanVitality für MultiGraph mit FreemandIndex mit Hilfe der Eigenvector Centrality
    vitalityDic = {}
    if speakersList[0] in graph:
        for node in speakersList:
            neighborss = list(graph.neighbors(node))
            adjList = []
            for neighbor in neighborss:
                adjList.append((node, neighbor))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexEigenvectorCentrality(graph)
            graph.add_edges_from(adjList)

    else:
        for node in speakersList:
            instanceOfNode = []
            adjList = []
            for timeScene in range(numberOfScenes):
                instanceOfNode.append(node + "~@" + str(i))
            for timeNode in instanceOfNode:
                neighborss = list(graph.neighbors(timeNode))
                # inEdges = list(graph.in_edges(timeNode))
                for neighbor in neighborss:
                    adjList.append((timeNode, neighbor))
                # for inEdge in inEdges:
                # adjList.append((inEdge[0], timeNode))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexEigenvectorCentrality(graph)
            graph.add_edges_from(adjList)
    return vitalityDic

def FreemanVitalityEigenvectorCentrality(graph, speakersList, numberOfScenes):
    vitalityDic = {}
    if speakersList[0] in graph:
        for node in speakersList:
            neighborss = list(graph.neighbors(node))
            adjList = []
            for neighbor in neighborss:
                adjList.append((node, neighbor))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexEigenvectorCentrality(graph)
            graph.add_edges_from(adjList)

    else:
        for node in speakersList:
            instanceOfNode = []
            adjList = []
            for timeScene in range(numberOfScenes):
                instanceOfNode.append(node + "~@" + str(i))
            for timeNode in instanceOfNode:
                neighborss = list(graph.neighbors(timeNode))
                inEdges = list(graph.in_edges(timeNode))
                for neighbor in neighborss:
                    adjList.append((timeNode, neighbor))
                for inEdge in inEdges:
                    adjList.append((inEdge[0], timeNode))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexEigenvectorCentrality(graph)
            graph.add_edges_from(adjList)
    return vitalityDic

def PlotFreemanVitalityDegreeCentrality(graph, title):
    freemanVitality = FreemanVitalityDegreeCentrality(graph, speakerList, len(sceneList)) ##TODO keine globalen werte
    nodeList = []
    valueList = []
    for node in freemanVitality:
        nodeList.append(node)
        valueList.append(freemanVitality[node])
    fiPlot = plt.figure('FI', dpi=200)
    fiAxis = fiPlot.add_subplot()
    fiAxis.set_ylabel('Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', linewidth=0.5, markersize=4, label='Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', marker='.', markersize=4, linestyle='None')
    fiAxis.plot(nodeList, [FreemanIndexDegreeCentrality(graph)]*len(nodeList), linestyle='--', linewidth=1, color='orange')
    plt.xticks(rotation=90)
    plt.title(title)
    plt.show()
    return

def PlotFreemanVitalityMultiGraphDegreeCentrality(graph, title):
    freemanVitality = FreemanVitalityMultiGraph(graph, speakerList, len(sceneList)) ##TODO keine globalen werte
    nodeList = []
    valueList = []
    for node in freemanVitality:
        nodeList.append(node)
        valueList.append(freemanVitality[node])
    fiPlot = plt.figure('FI', dpi=200)
    fiAxis = fiPlot.add_subplot()
    fiAxis.set_ylabel('Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', linewidth=0.5, markersize=4, label='Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', marker='.', markersize=4, linestyle='None')
    fiAxis.plot(nodeList, [FreemanIndexDegreeCentrality(graph)]*len(nodeList), linestyle='--', linewidth=1, color='orange')
    plt.xticks(rotation=90)
    plt.title(title)
    plt.show()
    return

def PlotFreemanVitalityMultiGraphEigenvectorCentrality(graph, title):
    freemanVitality = FreemanVitalityMultiGraphEigenvectorCentrality(graph, speakerList, len(sceneList)) ##TODO keine globalen werte
    nodeList = []
    valueList = []
    for node in freemanVitality:
        nodeList.append(node)
        valueList.append(freemanVitality[node])
    fiPlot = plt.figure('FI', dpi=200)
    fiAxis = fiPlot.add_subplot()
    fiAxis.set_ylabel('Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', linewidth=0.5, markersize=4, label='Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', marker='.', markersize=4, linestyle='None')
    fiAxis.plot(nodeList, [FreemanIndexEigenvectorCentrality(graph)]*len(nodeList), linestyle='--', linewidth=1, color='orange')
    plt.xticks(rotation=90)
    plt.title(title)
    plt.show()
    return

def PlotFreemanVitalityEigenvectorCentrality(graph, title):
    freemanVitality = FreemanVitalityEigenvectorCentrality(graph, speakerList, len(sceneList)) ##TODO keine globalen werte
    nodeList = []
    valueList = []
    for node in freemanVitality:
        nodeList.append(node)
        valueList.append(freemanVitality[node])
    fiPlot = plt.figure('FI', dpi=200)
    fiAxis = fiPlot.add_subplot()
    fiAxis.set_ylabel('Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', linewidth=0.5, markersize=4, label='Freeman Index')
    fiAxis.plot(nodeList, valueList, color='grey', marker='.', markersize=4, linestyle='None')
    fiAxis.plot(nodeList, [FreemanIndexEigenvectorCentrality(graph)]*len(nodeList), linestyle='--', linewidth=1, color='orange')
    plt.xticks(rotation=90)
    plt.title(title)
    plt.show()
    return


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
    if(len(newlist)>0):
        sceneList.append(newlist)

# finde insgesamt alle vorkommenden Personen
speakerList = []
for sp_tag in root.findall(".//tei:sp[@who]", namespaces=NS):
    value = sp_tag.get('who')
    for name in checkForDoubles(value):
        if name not in speakerList and value is not None:
            speakerList.append(name)


#Dynamischer Graph
dynamicGraph = nx.DiGraph()
for i in range(len(sceneList)):
    dynamicGraph.add_nodes_from(tempus(i, speakerList))
    tempusI = tempus(i, sceneList[i])
    dynamicGraph.add_edges_from(permutation(tempusI))
for i in range(len(sceneList) - 1):
    for j in speakerList:
        dynamicGraph.add_edge(j + "~@" + str(i), j + "~@" + str(i + 1))


#MultiGraph
multiGraph = nx.MultiGraph()
for i in range(len(sceneList)):
    multiGraph.add_nodes_from(tempus(i, speakerList))
    tempusI = tempus(i, sceneList[i])
    multiGraph.add_edges_from(permutation(tempusI))

for i in range(len(sceneList) - 1):
    for j in speakerList:
        multiGraph.add_edge(j + "~@" + str(i), j + "~@" + str(i + 1))

#Aggregations Graph
aggregateGraph = nx.Graph()
aggregateGraph.add_nodes_from(speakerList)
for singleScene in sceneList:
    permutationList = permutationSingle(singleScene)
    for (u, v) in permutationList:
        if aggregateGraph.has_edge(u, v):
            currentWeight = aggregateGraph[u][v]["weight"]
            aggregateGraph.add_edge(u,v, weight=currentWeight+1)
        else:
            aggregateGraph.add_edge(u,v, weight=1)



rgb1=177
rgb2=224
rgb3=185

plt.figure(1, dpi=400)
aggregateGraph_Layout = nx.spring_layout(aggregateGraph)
nx.draw_networkx(aggregateGraph, node_size=adjustSize(aggregateGraph), node_color=adjustColor(aggregateGraph),
                 width=0.2, with_labels=False, pos=aggregateGraph_Layout)
nx.draw_networkx_labels(aggregateGraph, pos=aggregateGraph_Layout)
plt.title("aggregate network of "+document.replace("input\\", ""))
plt.show()

PlotFreemanVitalityDegreeCentrality(aggregateGraph, "Freeman Index with Degree Centrality of Aggregate Graph")
PlotFreemanVitalityDegreeCentrality(dynamicGraph, "Freeman Index with Degree Centrality of Dynamic Graph")
#PlotFreemanVitalityDegreeCentrality(multiGraph, "Freeman Index of multi graph")

PlotFreemanVitalityEigenvectorCentrality(aggregateGraph, "Freeman Index with Eigenvector Centrality of Aggregate Graph")
PlotFreemanVitalityEigenvectorCentrality(dynamicGraph, "Freeman Index Eigenvector Centrality of Dynamic Graph")
PlotFreemanVitalityMultiGraphEigenvectorCentrality(multiGraph, "Freeman Index with Eigenvector Centrality of Multi Graph")
#PlotFreemanVitalityEigenvectorCentrality(multiGraph, "Freeman Index 2 of multi graph")

PlotFreemanVitalityMultiGraphDegreeCentrality(multiGraph, "Freeman Index with Degree Centrality of Multi Graph")
# #Matrizen
# adjMatrix1 = nx.to_numpy_matrix(dynamicGraph)
# print(adjMatrix1)
#
# adjMatrix2 = nx.adjacency_matrix(dynamicGraph).todense()
# print(adjMatrix2)
#
# aggMatrix = nx.to_numpy_matrix(aggregateGraph, nodelist=speakerList)
# print(aggMatrix)




