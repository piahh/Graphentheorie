import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import networkx as nx
from operator import itemgetter


def DocumentInput():
    """Gewünschtes Dokument wird abgefragt und kann aufgerufen werden.
    Fehlermeldung, wenn Pfad falsch angegeben."""
    # auf das Dokument kann dann beim Erstellen und Plotten des Graphen zugegriffen werden.
    # hieraus werden die Sprecher und Szenen gezogen
    # Grundlage
    try:
        document = input("Pfad zum xml des Dramas: ") # Pfad des XML Dokuments einfügen. Wenn im gleichen Ordner reicht input\dokument.xml
                                                      # sonst kompletten Pfad einfügen.
        #   'input\\arnim-das-loch.xml' ## VOLLSTÄNDIGER GRAPH
        #   'input\\alberti-brot.xml'  ## großes Netzwerk
        #   'input\\beer-der-paria.xml'  ## kleines Netzwerk

        AggregateGraph(document)
        return document
    except:
        print("Ungültiger Dateipfad")
        return ""


def SceneList(document):
    """findet pro Szene vorkommende Personen und speichert diese als Liste in einer Liste.
    Darauf kann dann zugegriffen werden, um den Graphen zu visualisieren."""
    #document ist Pfad zum xml
    root = ET.parse(document).getroot()
    NS = {'tei': 'http://www.tei-c.org/ns/1.0'}

    sceneList = []
    for div_tag in root.findall(".//tei:div", namespaces=NS):
        newlist = []
        for sp_tag in div_tag.findall(".//tei:sp[@who]", namespaces=NS):
            for name in CheckForDoubles(sp_tag.get('who')):
                if name not in newlist:
                    newlist.append(name)
        if (len(newlist) > 0):
            sceneList.append(newlist)
    return sceneList


def SpeakerList(document):
    """findet  alle vorkommenden Personen, unabhängig von der Szene."""
    # document ist Pfad zum xml

    root = ET.parse(document).getroot()
    NS = {'tei': 'http://www.tei-c.org/ns/1.0'}

    speakerList = []
    for sp_tag in root.findall(".//tei:sp[@who]", namespaces=NS):
        value = sp_tag.get('who')
        for name in CheckForDoubles(value):
            if name not in speakerList and value is not None:
                speakerList.append(name)
    return speakerList


def Permutation(liste):
    """Gibt Permutationen der ListeneElemente aus, ohne (name, name). Eine Liste an Tupeln wird erstellt,
    die alle Verbindungen zwischen den Personen angibt."""
    tupelListe = []
    for name in liste:
        for name2 in liste:
            if name != name2:
                tupelListe.append((name, name2))
    return tupelListe


def Tempus(time, person):
    """Liste der Personen + die Zeit/Scenen-Information als Zahl hinter den Namen gehängt,
    damit eindeutig zugewiesen werden kann."""
    newlist = []
    for name in person:
        name2 = name + "~@" + str(time)
        newlist.append(name2)
    return newlist


def CheckForDoubles(String):
    """Prüft ob mehrere Personen gleichzeitig sprechen und teilt diese dann auf einzelne Einträge auf.
    Rückgabe Liste der sprechenden Pers. : "#Anna #Berta" -> ["#Anna", "#Berta"]"""
    stringList = []
    y = len(String)
    x = 0
    for i in range(len(String) - 1, -1, -1):
        if String[i] == "#":
            x = i
            stringList.append(String[x:y])
            y = x - 1
    return stringList


def PermutationSingle(liste):
    """Wie Permutation, aber falls (name, name2) schon existiert, wird (name2, name) nicht noch hinzugefügt"""
    tupelListe = []
    for name in liste:
        for name2 in liste:
            if name != name2:
                if (name2, name) not in tupelListe:
                    tupelListe.append((name, name2))
    return tupelListe


def AggregateGraph(document):
    """Aggregate Graph wird erstellt. """
    #greift auf document (behandeltes Drama) zu, Liste der Sprecher allgemein.
    #Kanten sind gewichtet, um Wichtigkeit einer Person anzuzeigen.
    aggregateGraph = nx.Graph()
    aggregateGraph.add_nodes_from(SpeakerList(document))
    for singleScene in SceneList(document):
        permutationList = PermutationSingle(singleScene)
        for (u, v) in permutationList:
            if aggregateGraph.has_edge(u, v):
                currentWeight = aggregateGraph[u][v]["weight"]
                aggregateGraph.add_edge(u, v, weight=currentWeight + 1)
            else:
                aggregateGraph.add_edge(u, v, weight=1)
    return aggregateGraph


def DynamicGraph(document):
    """dynamischer Graph wird erstellt."""
    #inklusive Zeitinformation, durch Szenen.
    dynamicGraph = nx.Graph()
    for i in range(len(SceneList(document))):
        dynamicGraph.add_nodes_from(Tempus(i, SpeakerList(document)))
        tempusI = Tempus(i, SceneList(document)[i])
        dynamicGraph.add_edges_from(Permutation(tempusI), weight=1)

    for i in range(len(SceneList(document)) - 1):
        for j in SpeakerList(document):
            dynamicGraph.add_edge(j + "~@" + str(i), j + "~@" + str(i + 1), weight=1)
    return dynamicGraph


def PlottedAggregateGraph(aggregateGraph, document):
    """aggregats Graph wird geplottet."""
    plt.figure(1, dpi=400)
    # nur die wichtigsten Personen des Dramas.
    # Anzahl im Verhältnis zur Länge des Dramas.
    mostImportantNodes = list(reversed(sorted(WeightedDegreeTuple(aggregateGraph),
                                              key=itemgetter(1))))[0:round(len(aggregateGraph)**0.67)+1]
    labelDic = {} # nur die wichtigsten Personen werden in das Dictionaryhinzugefügt.
    for x,y in mostImportantNodes:
        labelDic[x] = x.replace("#", "").replace("_", " ") # Raute und Unterstrich werden, im Label der wichtigsten/angezeigten Personen, entfernt

    aggregateGraph_Layout = nx.spring_layout(aggregateGraph)
    nx.draw_networkx(aggregateGraph, node_size=AdjustSize(aggregateGraph), node_color=AdjustColor(aggregateGraph),
                     width=0.2, with_labels=False, pos=aggregateGraph_Layout)
    # nur die wichtigsten Namen (gespeichert im Dictionary) werden im geplotteten Graphen angezeigt
    nx.draw_networkx_labels(aggregateGraph, labels=labelDic,  pos=aggregateGraph_Layout)
    # Name des verwendeten Dramas im Titel der Visualisierung
    plt.title("Aggregate Network of "+ShortenPath(document))
    plt.show()


def AdjustValue(dcurrent, degreevalues, values):
    """dcurrent = degree of current node
    degreevalues = List of degreeMin, degreeMean, degreeMax
    values = List of valueMin, valueMean, valueMax;
    values ist das Spektrum der Eigenschaft(zb. Größe) welche einem bestimmten Knoten(dcurrent) zugeordnet werden soll"""
    adjustedValues = []
    # Entspricht der Grad des Knotens dem Durchschnitt wird die durchschnittliche Ausprägung der Eigenschaft zurückgegeben
    if dcurrent == degreevalues[1]:
        return values[1]
    # Ist der Grad des Knotens größer als der Durchschnitt wird eine, im gleichen Verhältniss stehende, Ausprägung der
    # Eigenschaft zurück gegeben: d.h. die Differenz zwischen maxAusprägung und meanAusprägung
    # wird mit (dcurrent - mean)/(max - mean) multipliziert
    elif dcurrent > degreevalues[1]:
        for i in range(len(values[0])):
            adjustedValues.append((values[1][i] + (values[2][i] - values[1][i]) * ((
                    dcurrent - degreevalues[1]) / (
                                           degreevalues[2] - degreevalues[1]))))
    else:
        for i in range(len(values[0])):
            adjustedValues.append((values[1][i] - (values[1][i] - values[0][i]) * ((
                    degreevalues[1] - dcurrent ) / (
                                           degreevalues[1] - degreevalues[0] + 1))))
    return tuple(adjustedValues)


def AdjustColor(graph):
    """Nutzt AdjustValue um eine Liste der Farben der Knoten zurückzugeben
    Farbe der Knoten ist abhängig von der Relevanz einer Person, gemessen an der Häufigkeit des
    Vorkommens in unterschiedlichen Szenen und somit der Verbindung zu anderen Sprechern. """
    degreeTuple = WeightedDegreeTuple(graph)

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
        colorTime.append(AdjustValue(dvalue, degreeValues, colors))
    return colorTime


def AdjustSize(graph):
    """Nutzt AdjustValue um eine Liste mit Größen der Knoten zurück zu geben
    # Größe der Knoten ist abhängig von der Relevanz einer Person, gemessen an der Häufigkeit des Vorkommens in unterschiedlichen Szenen."""
    degreeTuple = WeightedDegreeTuple(graph)

    sizeTime = []
    sMax = 300,
    sMean = 150,
    sMin = 50,
    # ',' am Ende macht das Objekt zu einem 1-Tupel, so kann in AdjustValue per Index darauf zugegriffen werden
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
        sizeTime.append(AdjustValue(dvalue, degreeValues, sizes))
    return sizeTime


def WeightedDegreeTuple(graph):
    '''Nutzt Summe des Gewichts der umliegenden Kanten zum Berechnen der passenden Größe/Farbe, etc'''
    degreeTuple = []
    graphAttributes = nx.get_edge_attributes(graph, "weight")
    weightedDegreeTuple = {}
    for node in list(graph):
        weightedDegreeTuple[node] = 0
    for (u,v) in graphAttributes:
        weightedDegreeTuple[u] += graphAttributes[(u,v)]
        weightedDegreeTuple[v] += graphAttributes[(u,v)]
    for node in list(graph):
        degreeTuple.append((node, weightedDegreeTuple[node]))
    return degreeTuple


def ShortenPath(string):
    """kürzst den Pfad beim Aufrufen des Dokumentes, um ihn später als Namen in die Überschrift der Graphik einfügen zu können."""
    for i in range(len(string)-1, -1, -1):
        if string[i] == '/' or string[i] == '\\':
            return string[i+1:len(string)]
    return string


def FreemanIndexEigenvectorCentrality(graph):
    """Berechnet den FreemanIndex über den Grad des Knoten mit Hilfe der Eigenvector Centrality"""
    eigen = dict(nx.eigenvector_centrality_numpy(graph, weight="weight"))

    degreeMax = 0
    sum = 0
    for value in eigen.values():
        if value != 'weight' and degreeMax < value:
            degreeMax = value
    for nodes in graph:
        sum += (degreeMax - eigen[nodes])

    # Erzeugt Sterngraph in gleicher größe wie 'graph' um zentralisieren zu können
    # Nach Zweig (2016: 259) und Prado et al. (2016).
    referenzSternGraph = nx.Graph()
    for i in range(1, len(graph)):
        referenzSternGraph.add_edge(0, i)
    eigenRef = dict(nx.eigenvector_centrality_numpy(referenzSternGraph))
    degreeMaxRef = 0
    sumRef = 0
    for valueRef in eigenRef.values():
        if degreeMaxRef < valueRef:
            degreeMaxRef = valueRef
    for nodesRef in referenzSternGraph:
        sumRef += (degreeMaxRef - eigenRef[nodesRef])
    return sum / sumRef


def FreemanVitality(graph, speakersList, numberOfScenes):
    """FreemanVitality für Graph mit FreemanIndex"""
    vitalityDic = {}
    #freemanI = FreemanIndexEigenvectorCentrality(graph)
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
                instanceOfNode.append(node + "~@" + str(timeScene))
            for timeNode in instanceOfNode:
                neighborss = list(graph[timeNode])
                for neighbor in neighborss:
                    adjList.append((timeNode, neighbor))
            graph.remove_edges_from(adjList)
            vitalityDic[node] = FreemanIndexEigenvectorCentrality(graph)
            graph.add_edges_from(adjList)
    return vitalityDic


def PlotFreemanVitality(graph, title, nodelist, numberOfScenes):
    """FreemanVitality wird geplottet"""
    # nodelist = zugrundeliegende Knoten
    freemanVitality = FreemanVitality(graph, nodelist, numberOfScenes)
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
    fiAxis.plot(nodeList, [FreemanIndexEigenvectorCentrality(graph)] * len(nodeList), linestyle='--', linewidth=1,
                color='orange')
    plt.xticks(rotation=90)
    plt.title(title)
    plt.show()
    return


def _main():
    """Zuerst wird das Dokument erfragt und kann eingefügt werden (documentInput).
        Bei falscher Eingabe erscheint eine Fehlermeldung und erneut die Aufforderung den Pfad des Dramas/Dokumentes anzugeben.
    Aus dem Drama werden Sprecher und Szenen extrahiert und jeweils in Listen gespeichert. (SceneList, SpeakerList)
        Die Sprecher stellen Knoten des Netzwerkes dar, die Szenen zeigen die Verbindungen zwischen den Sprechern (Permutation),
        zu unterschiedlichen Zeitpunkten an. (Tempus)
        Personengruppen werden in einzelne Personen aufgeteilt. (CheckForDoubles)
        Wenn eine Verbindung von a zu b besteht, wird die Verbindung von b zu a nicht noch zusätzlich aufgeführt (PermutationSingle)
    Auf Grundlage des ausgewählten Dramas werden ein aggregats und ein dynamischer Graph erstellt, sowie der aggregats Graph geplottet.
                                                                                (AggregateGraph, DynamicGraph, PlottedAggregateGraph).
    Die Visualisierung des Graphen zeigt die Wichtigkeit der Sprecher des Dramas im Verhältnis zueinander an,
        mit Hilfe von Einfärbung und Größe der Knoten. (AdjustSize, AdjustColor)
        Dies funktioniert über Gewichtung (AdjustValue, WeightedDegreeTuple)
        Im Titel der Graphik wird der gekürzte Pfadname des Dramas mitangegeben. (ShortenPath)
    Auf Grundlage der Graphen und somit des Dramas, wird der FreemanIndex des aggregats Graphen berechnet.
        Hierfür wird die Eigenvector Zentralität der Knoten hinzugezogen. Grundlage sind die Paper von Prado et. al. (2016) und Zweig (2016)
        (FreemanIndexEigenvectorCentrality)
        der FreemanIndex allgemein wird als orange Linie in der Graphik dargestellt. (PlotFreemanVitality)
        Über die Vitality wird der Unterschied deutlich, wie sich der Index ändert,
        wenn einzelne Personen/ihre Kanten, aus dem Graphen gelöscht werden. (FreemanVitality)
        """
    document = DocumentInput()
    while document == "":
        document = DocumentInput()

    aggregateGraph = AggregateGraph(document)
    dynamicGraph = DynamicGraph(document)
    PlottedAggregateGraph(aggregateGraph, document)
    # FreemanIndex wird als Graphik geplottet
    PlotFreemanVitality(aggregateGraph, "Freeman Index with Eigenvector Centrality of Aggregate Graph", SpeakerList(document),
                        len(SceneList(document)))
    PlotFreemanVitality(dynamicGraph, "Freeman Index Eigenvector Centrality of Dynamic Graph", SpeakerList(document),
                        len(SceneList(document)))
    return

if __name__ == '__main__':
    _main()


