from GRNCreationUtils import meanClustering, BarabasiAlbertAlgorithm
from GRNCreationUtils import adjacenteDiMatriceStaredFromGraph, addColors
from ochunGRN.GRN.genesGroup import subgraph3N
import networkx as nx
import numpy as np


def randomGrn(genesNb: int,
              autoRG: float,
              duoRG: float,
              resDict: dict | None = None) -> dict:
    """
    Generate a random Gene Regulatory Network (GRN).

    This function generates a GRN using the Barabasi-Albert algorithm
    to create a graph, then calculates the adjacency matrix
    and other graph characteristics, such as the average clustering
    coefficient and a 3-node subgraph.

    Parameters:
    - genesNb (int): The number of genes (nodes) in the network.
    - autoRG (float): The rate of self-regulation in the network.
    - duoRG (float): The regulation rate between pairs of genes.
    - resDict (dict, optional): A dictionary to store the results. If no dictionary is provided, a new one is created.
    # noqa: E501

    Returns:
    - dict: A dictionary containing the generated graph, adjacency matrix,
    average clustering coefficient, and a subgraph.
    """
    if resDict is None:
        resDict = {}

    # Generate the initial graph using the Barabasi-Albert model
    Graph = BarabasiAlbertAlgorithm(genesNb, 2)

    # Store results in the provided or new dictionary
    Graph, M = adjacenteDiMatriceStaredFromGraph(Graph, autoRG, duoRG)
    resDict["Graph"] = Graph
    resDict["genesNb"] = genesNb
    resDict["autoRG"], resDict["duoRG"] = autoRG, duoRG
    resDict["AdjMatrice"] = M
    resDict["meanClustering"] = meanClustering(Graph)
    resDict["subGraph"] = subgraph3N(Graph)

    return resDict


def GrnFromAdj(AdjMatrice: np.ndarray,
               resDict: dict = None) -> dict:
    """
    Create a Gene Regulatory Network (GRN) from an adjacency matrix.

    This function constructs a directed graph from a given adjacency matrix,
    then calculates the number of self-regulations and other graph
    characteristics.

    Parameters:
    - AdjMatrice (numpy.ndarray): The adjacency matrix representing the GRN.
    - resDict (dict, optional): A dictionary to store the results. If no dictionary is provided, a new one is created.
    # noqa: E501
    Returns:
    - dict: A dictionary containing the generated graph, the number of genes,
    the number of self-regulations, and other characteristics.
    """
    if resDict is None:
        resDict = {}

    # Create directed graph from adjacency matrix
    genesNb = len(AdjMatrice)
    Graph = nx.from_numpy_array(AdjMatrice, create_using=nx.DiGraph)

    # Calculate self-regulation and pairwise regulation rates
    autoRG = 0
    for i in range(genesNb):
        if AdjMatrice[i][i] != 0:
            autoRG += 1
    N = Graph.number_of_edges()
    N -= autoRG
    duoRG = 0
    for i in range(genesNb):
        for j in range(i+1, genesNb):
            if AdjMatrice[i][j] != 0 and AdjMatrice[j][i] != 0:
                N -= 1
                duoRG += 1
    autoRG /= Graph.number_of_edges()
    duoRG /= N

    # Add colors and store results
    addColors(Graph, AdjMatrice)
    resDict["Graph"] = Graph
    resDict["genesNb"] = genesNb
    resDict["autoRG"], resDict["duoRG"] = autoRG, duoRG
    resDict["AdjMatrice"] = AdjMatrice
    resDict["meanClustering"] = meanClustering(Graph)
    resDict["subGraph"] = subgraph3N(Graph)

    return resDict


def main():
    pass


if __name__ == "__main__":
    main()
