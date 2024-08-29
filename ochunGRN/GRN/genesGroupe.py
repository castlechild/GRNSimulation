import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


def findGroupGraph(adj_matrix):
    """
    Identifie le type de sous-graphe à trois nœuds (motif) dans un graphe
    dirigé en utilisant sa matrice d'adjacence.

    Cette fonction prend une matrice d'adjacence de sous-graphe à trois nœuds
    et la classifie en fonction de ses motifs structurels parmi les options
    suivantes :
    "Fan-In", "Fan-Out", "Cascade", "FFL", "Triangular", "C-D", "C-E",
    "SCascade", "dC", "dD", "dE", "5edges", et "6edges".

    Paramètres :
    - adj_matrix (numpy.ndarray): Matrice d'adjacence 3x3 représentant
    le sous-graphe.

    Retourne :
    - str : Le type de motif détecté parmi les options disponibles.
    """
    if adj_matrix.shape != (3, 3):
        raise ValueError("La matrice d'adjacence doit être de taille 3x3.")
    for i in range(3):
        adj_matrix[i, i] = 0
    nEdges = np.sum(adj_matrix)

    if nEdges == 2:
        for i in range(3):
            if np.sum(adj_matrix[:, i]) == 2:
                return "Fan-In"
            elif np.sum(adj_matrix[i, :]) == 2:
                return "Fan-Out"
        return "Cascade"

    elif nEdges == 3:
        for i in range(3):
            for j in range(3):
                if i != j and adj_matrix[i, j] == 0 and adj_matrix[j, i] == 0:
                    if (np.sum(adj_matrix[:, i]) == 1 and
                            np.sum(adj_matrix[:, j]) == 1):
                        return "C-D"
                    else:
                        return "C-E"
        for i in range(3):
            if np.sum(adj_matrix[:, i]) == 2:
                return "FFL"
        return "Triangular"

    elif nEdges == 4:
        for i in range(3):
            for j in range(3):
                if i != j and adj_matrix[i, j] == 0 and adj_matrix[j, i] == 0:
                    return "SCascade"
        for i in range(3):
            if np.sum(adj_matrix[i, :]) + np.sum(adj_matrix[:, i]) == 2:
                if np.sum(adj_matrix[:, i]) == 2:
                    return "dE"
                elif np.sum(adj_matrix[i, :]) == 2:
                    return "dD"
                return "dC"
        raise ValueError("Motif non classifié pour nEdges = 4")

    elif nEdges == 5:
        return "5edges"

    elif nEdges == 6:
        return "6edges"

    raise ValueError("Motif non classifié")


def subgraph3N(Graph):
    resDic = {}
    nodes = list(Graph.nodes)
    edges = list(Graph.edges)
    adj_matrice = nx.to_numpy_array(Graph)
    adj_matrice_undi = np.zeros(np.shape(adj_matrice))
    DictPos = {}
    cache = set()
    for i in range(len(nodes)):
        DictPos[nodes[i]] = i
        for j in range(len(nodes)):
            if adj_matrice[i][j] != 0:
                adj_matrice[i][j] = 1
                adj_matrice_undi[i][j] = 1
                adj_matrice_undi[j][i] = 1
    print(edges)
    for u, v in edges:
        for w in set(Graph[v]).union(Graph.pred[v], Graph[u], Graph.pred[u]):
            if u != v != w != u:
                triplet = tuple(sorted([w, u, v]))
                if triplet not in cache:
                    cache.add(tuple(sorted([w, u, v])))
                    if issubgraphconnected(adj_matrice_undi, DictPos[u],
                                           DictPos[v], DictPos[w]):
                        adj_submatrice = adj_matrice[
                            np.ix_([DictPos[u], DictPos[v], DictPos[w]],
                                   [DictPos[u], DictPos[v], DictPos[w]])]
                        resDic[f"{u}-{v}-{w}"] = findGroupGraph(
                            adj_submatrice)
    return resDic


def issubgraphconnected(adj_matrice, nodeEdgeA, nodeEdgeB, node):
    return bool(adj_matrice[node][nodeEdgeA]) or bool(adj_matrice[node][nodeEdgeB])  # noqa: E501


def main():
    AdjMatrice = np.array([[0, 0, -1, 0, 0, 0, 0, 0],
                           [0, 0, -1, 0, -1, -1, 0, 0],
                           [0, 1, 0, 1, 0, -1, -1, 0],
                           [0, -1, 1, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0, 0, 0],
                           [-1, 0, 1, 0, 0, 0, 0, 1],
                           [0, 0, 0, 0, 1, 0, 0, 0]])
    G = nx.from_numpy_array(AdjMatrice, create_using=nx.DiGraph)
    nx.draw_circular(G, with_labels=True)
    print(subgraph3N(G))

    plt.show()
    TGraph = [nx.DiGraph() for i in range(13)]
    TGraph[0].add_edges_from([(0, 1), (1, 2), (2, 0)])
    TGraph[1].add_edges_from([(1, 0), (1, 2), (2, 0)])
    TGraph[2].add_edges_from([(0, 1), (1, 2)])
    TGraph[3].add_edges_from([(1, 0), (1, 2)])
    TGraph[4].add_edges_from([(0, 1), (2, 1)])
    TGraph[5].add_edges_from([(0, 1), (1, 2), (2, 0), (0, 2)])
    TGraph[6].add_edges_from([(1, 0), (1, 2), (2, 0), (0, 2)])
    TGraph[7].add_edges_from([(0, 1), (2, 1), (2, 0), (0, 2)])
    TGraph[8].add_edges_from([(0, 1), (1, 2), (2, 1)])
    TGraph[9].add_edges_from([(0, 1), (1, 2), (1, 0)])
    TGraph[10].add_edges_from([(0, 1), (1, 2), (1, 0), (2, 1)])
    TGraph[11].add_edges_from([(0, 1), (1, 2), (1, 0), (2, 1), (2, 0)])
    TGraph[12].add_edges_from([(0, 1), (1, 2), (1, 0), (2, 1), (2, 0), (0, 2)])
    for i in range(3):
        plt.subplot(3, 5, i+2)
        nx.draw_circular(TGraph[i], arrowsize=30)
        plt.title(findGroupGraph(nx.to_numpy_array(TGraph[i])))
    for i in range(5):
        plt.subplot(3, 5, i+6)
        nx.draw_circular(TGraph[i+3], arrowsize=30)
        plt.title(findGroupGraph(nx.to_numpy_array(TGraph[i+3])))
    for i in range(5):
        plt.subplot(3, 5, i+11)
        nx.draw_circular(TGraph[i+8], arrowsize=30)
        plt.title(findGroupGraph(nx.to_numpy_array(TGraph[i+8])))
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.savefig("Melvin/Images/allGroups.png")


if __name__ == "__main__":
    main()
