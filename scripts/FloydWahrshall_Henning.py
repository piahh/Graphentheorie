import networkx as nx
import numpy as np


def FloydWahrshall(matrix):
    """Akzeptiert eine AdjMatrix im Format [[],[]], wobei ein Eintrag die Länge der Kante beschreibt.
    Berechnet den Kürzesten Pfad zwischen zwei Knoten.
    In meinem Beispiel aufgerufen und importiert in GermanDrama_Henning und angewandt auf ein Drama und den dazugehörigen aggregats Graphen,
    kann aber auf beliebigen Graphen oder Matrix angewandt werden."""
    V = len(matrix)
    shortestPath = matrix
    for k in range(V):
        shortestPathI = [list(row) for row in shortestPath]
        for l in range(V):
            for j in range(V):
                shortestPathI[j][l] = min(shortestPath[j][l], shortestPath[l][k] + shortestPath[k][j])
        shortestPath = shortestPathI
    return shortestPath


def GraphToMatrix(graph):
    """Wandelt einen Networkx-Graphen in eine einfache Matrix für den FloydWahrshall-Algorithmus um,
    falls ein Graph eingegeben wird."""
    # 0 wird durch inf ersetzt, da hier keine Verbindung vorliegt.
    adjMatrix = nx.to_numpy_matrix(graph).tolist()
    for x in range(len(adjMatrix)):
        for y in range(len(adjMatrix)):
            if adjMatrix[x][y] == 0:
                adjMatrix[x][y] = np.inf
    return adjMatrix



# Beispiel Matrix
matrix = [[2, 3, 4], [4, 5, 2], [3, 1, 4]]

print(FloydWahrshall(matrix))
print(FloydWahrshall([[2, 3, 4], [4, 5, 2], [3, 1, 4]]))

## Beispiel Graph
Nodes = ["a", "b", "c", "d", "e"]
Edges = [("a", "b"), ("a", "c"), ("c", "d"), ("d", "e"), ("e", "b")]
graph = nx.Graph()
graph.add_edges_from(Edges)
graph.add_nodes_from(Nodes)

print(GraphToMatrix(graph))