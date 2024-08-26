import networkx as nx
from itertools import combinations
import multiprocessing
# from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd

document = pd.read_excel("ochunGRN/GRN/41598_2021_3625_MOESM5_ESM.xlsx")

arabidopsisThaliana = (document["Supplementary Table S1: Networks. A spreadsheet file with filtered networks"].tolist()[2:], document["Unnamed: 1"].tolist()[2:], "Arabidopsis thaliana")  # noqa: E501
drosophilaMelanogaster = (document["Unnamed: 2"].tolist()[2:], document["Unnamed: 3"].tolist()[2:], "Drosophila Melanogaster")  # noqa: E501
escherichniaColi = (document["Unnamed: 4"].tolist()[2:], document["Unnamed: 5"].tolist()[2:], "Escherichnia coli")  # noqa: E501
homoSapiens = (document["Unnamed: 6"].tolist()[2:], document["Unnamed: 7"].tolist()[2:], "Homo sapiens")  # noqa: E501
saccharomycesCerevisiae = (document["Unnamed: 8"].tolist()[2:], document["Unnamed: 9"].tolist()[2:], "Saccharomyces cerevisiae")  # noqa: E501


def findGroupGraph(threeNGraph):
    threeNGraphWAuto = nx.DiGraph()
    edges = threeNGraph.edges()
    for u, v in edges:
        if u != v:
            threeNGraphWAuto.add_edge(u, v)

    edges = threeNGraphWAuto.edges()
    nEdges = threeNGraphWAuto.number_of_edges()
    if nEdges == 2:
        for u in threeNGraphWAuto:
            if threeNGraphWAuto.in_degree[u] == 2:
                return "Fan-In"
            elif threeNGraphWAuto.out_degree[u] == 2:
                return "Fan-Out"
        return "Cascade"

    elif nEdges == 3:
        for u in threeNGraphWAuto:
            for v in threeNGraphWAuto:
                if u != v and (u, v) not in edges and (v, u) not in edges:
                    if (threeNGraphWAuto.in_degree[u] == 1
                            and threeNGraphWAuto.in_degree[v] == 1):
                        return "C-D"
                    else:
                        return "C-E"
        for u in threeNGraphWAuto:
            if threeNGraphWAuto.in_degree[u] == 2:
                return "FFL"
        return "Triangular"

    elif nEdges == 4:
        for u in threeNGraphWAuto:
            for v in threeNGraphWAuto:
                if u != v and (u, v) not in edges and (v, u) not in edges:
                    return "SCascade"
        for u in threeNGraphWAuto:
            if threeNGraphWAuto.degree[u] == 2:
                if threeNGraphWAuto.in_degree[u] == 2:
                    return "dE"
                elif threeNGraphWAuto.out_degree[u] == 2:
                    return "dD"
                return "dC"
        raise

    elif nEdges == 5:
        return "5edges"

    elif nEdges == 6:
        return "6edges"

    raise ValueError("Unclassified graph pattern")


def findGroupGraph2(adj_matrix):
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
    """
    Extract and analyze all three-node subgraphs from a given graph.

    This function traverses a given graph, extracts all subgraphs containing
    exactly three nodes, and uses the `findGroupGraph` function to classify
    each one.

    Parameters:
    - Graph (networkx.Graph): The graph to analyze.

    Returns:
    - dict: A dictionary where the keys are the types of patterns found and
    the values are the counts of each type.
    """
    resDic = {}
    nodes = list(Graph.nodes)
    total_combinations = len(nodes) * (len(nodes) - 1) * (len(nodes) - 2) // 6
    with tqdm(total=total_combinations) as pbar:
        for i, (u, v, w) in enumerate(combinations(nodes, 3)):
            if issubgraphconnected(Graph, u, v, w):
                Subgraph = nx.subgraph(Graph, [u, v, w])
                resDic[f"{u}-{v}-{w}"] = findGroupGraph(Subgraph)
            pbar.update(1)
    return resDic


def subgraph3N2(Graph):
    resDic = {}
    nodes = list(Graph.nodes)
    adj_matrice = nx.to_numpy_array(Graph)
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if adj_matrice[i][j] != 0:
                adj_matrice[i][j] = 1
    total_combinations = len(nodes) * (len(nodes) - 1) * (len(nodes) - 2) // 6
    with tqdm(total=total_combinations) as pbar:
        for i, (u, v, w) in enumerate(combinations(range(len(nodes)), 3)):
            if issubgraphconnected(adj_matrice, u, v, w):
                adj_submatrice = adj_matrice[np.ix_([u, v, w], [u, v, w])]
                resDic[f"{u}-{v}-{w}"] = findGroupGraph2(adj_submatrice)
            pbar.update(1)
    return resDic


def issubgraphconnected(adj_matrice, nodeA, nodeB, nodeC):
    i = 0
    i += int(adj_matrice[nodeA][nodeC]) | int(adj_matrice[nodeC][nodeA])
    i += int(adj_matrice[nodeB][nodeC]) | int(adj_matrice[nodeC][nodeB])
    i += int(adj_matrice[nodeA][nodeB]) | int(adj_matrice[nodeB][nodeA])
    return i >= 2


def worker(tasks, adj_matrice, resDic):
    for task in tasks:
        u, v, w = task
        if issubgraphconnected(adj_matrice, u, v, w):
            adj_submatrice = adj_matrice[np.ix_([u, v, w], [u, v, w])]
            resDic[f"{u}-{v}-{w}"] = findGroupGraph2(adj_submatrice)


def chunk_tasks(tasks, chunk_size):
    chunk = []
    for task in tasks:
        chunk.append(task)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def subgraph3N_parallel(Graph):
    manager = multiprocessing.Manager()
    resDic = manager.dict()
    nodes = list(Graph.nodes)
    adj_matrice = nx.to_numpy_array(Graph)

    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if adj_matrice[i][j] != 0:
                adj_matrice[i][j] = 1

    tasks = combinations(range(len(nodes)), 3)
    coreNB = 1
    maxTasks = 100
    chunks = chunk_tasks(tasks, maxTasks)
    processes = []
    total_combinations = len(nodes) * (len(nodes) - 1) * (len(nodes) - 2) // 6
    with tqdm(total=total_combinations) as pbar:
        for chunk in chunks:
            if len(processes) >= coreNB:
                processes[0].join()
                pbar.update(1)
                processes.pop(0)
            p = multiprocessing.Process(
                target=worker, args=(chunk, adj_matrice, resDic))
            processes.append(p)
            p.start()
        for p in processes:
            p.join()
            pbar.update(1)
    return resDic


def FFLratio(Gspecies):
    print(Gspecies.number_of_nodes(), Gspecies.number_of_edges())
    print(nx.is_weakly_connected(Gspecies))
    DictFFL = subgraph3N2(Gspecies)
    motifs = {v: 0 for v in DictFFL.values()}
    N = len(DictFFL)
    print(N)
    for i in DictFFL.keys():
        motifs[DictFFL[i]] += 1
    print(sorted(motifs.items(), key=lambda item: item[1]))
    print(motifs['FFL']/N)


def graphCreator(species):
    N = len(species[0])
    Graph = nx.DiGraph()
    for i in range(N):
        Graph.add_edge(species[0][i], species[1][i])
    return Graph


def main():
    print("ratioFFL", FFLratio(graphCreator(drosophilaMelanogaster)))


if __name__ == "__main__":
    main()
