import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import numpy.typing as npt


def findGroupGraph(adj_matrix: npt.ArrayLike) -> str:
    """
    Identifies the type of three-node subgraph (motif) in a directed graph
    using its adjacency matrix.

    This function takes an adjacency matrix of a three-node subgraph
    and classifies it based on its structural motifs.
    The classification includes types such as "Fan-In", "Fan-Out", "Cascade",
    "FFL", "FBL", "Mutual-Out", "Mutual-In", "Bi-Mutual", "Regulated-Mutual",
    "Regulating-Mutual", "Mutual-Cascade", "Semi-Clique" and "Clique".

    Parameters:
    - adj_matrix (np.ndarray): A 3x3 adjacency matrix representing the subgraph

    Returns:
    - str: The detected motif type from the available options.
    """
    # Ensure the matrix is 3x3
    if adj_matrix.shape != (3, 3):
        raise ValueError("The adjacency matrix must be 3x3.")

    # Remove self-loops
    for i in range(3):
        adj_matrix[i, i] = 0

    # Count the number of edges in the subgraph
    nEdges = np.sum(adj_matrix)

    #  Classify the motif based on the number of edges
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
                        return "Mutual-Out"
                    else:
                        return "Mutual-In"
        for i in range(3):
            if np.sum(adj_matrix[:, i]) == 2:
                return "FFL"
        return "FBL"

    elif nEdges == 4:
        for i in range(3):
            for j in range(3):
                if i != j and adj_matrix[i, j] == 0 and adj_matrix[j, i] == 0:
                    return "Bi-Mutual"
        for i in range(3):
            if np.sum(adj_matrix[i, :]) + np.sum(adj_matrix[:, i]) == 2:
                if np.sum(adj_matrix[:, i]) == 2:
                    return "Regulated-Mutual"
                elif np.sum(adj_matrix[i, :]) == 2:
                    return "Regulating-Mutual"
                return "Mutual-Cascade"
        raise ValueError("Motif non classifié pour nEdges = 4")

    elif nEdges == 5:
        return "Semi-Clique"

    elif nEdges == 6:
        return "Clique"

    raise ValueError("Motif non classifié")


def subgraph3N(Graph: nx.digraph.DiGraph) -> dict:
    """
    Identifies all three-node subgraphs (motifs) in a directed graph.

    This function takes a directed graph, finds all three-node subgraphs,
    and classifies each based on its motif type.

    Parameters:
    - Graph (networkx.DiGraph): A directed graph to analyze.

    Returns:
    - dict: A dictionary where keys are triplets of nodes and values are
    the detected motif type.
    """
    resDic = {}
    nodes = list(Graph.nodes)
    edges = list(Graph.edges)
    adj_matrice = nx.to_numpy_array(Graph)
    # the creation of undi matrice is use for finding connected subgraph
    adj_matrice_undi = np.zeros(np.shape(adj_matrice))
    DictPos = {}
    cache = set()

    # Create adjacency matrices for directed and undirected graphs
    for i in range(len(nodes)):
        DictPos[nodes[i]] = i
        for j in range(len(nodes)):
            if adj_matrice[i][j] != 0:
                adj_matrice[i][j] = 1
                adj_matrice_undi[i][j] = 1
                adj_matrice_undi[j][i] = 1

    # Identify all three-node subgraphs
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


def issubgraphconnected(adj_matrice: npt.ArrayLike,
                        nodeEdgeA: int,
                        nodeEdgeB: int,
                        node: int) -> bool:
    """
    Checks if a three-node subgraph is connected. (knowing an edge)

    Parameters:
    - adj_matrice (numpy.ndarray): Adjacency matrix of the graph.
    - nodeEdgeA (int): Index of the first edge node.
    - nodeEdgeB (int): Index of the second edge node.
    - node (int): Index of a third node.

    Returns:
    - bool: True if the subgraph is connected, False otherwise.
    """
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
    plt.suptitle("13 types of 3-nodes subgraphs")
    plt.suptitle("An introduction to systems biology (Alon 2019)",
                 fontsize=20, x=0.8, y=0.08)
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.show()


if __name__ == "__main__":
    main()
