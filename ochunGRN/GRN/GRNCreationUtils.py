#!/usr/bin/env python3

import networkx as nx
import random as rd
import matplotlib.pyplot as plt
import numpy as np

# Barabasi-Albert algorithm
# Ma = np.zeros((genesNb,genesNb))


def BarabasiAlbertAlgorithm(n: int, an: int) -> nx.graph.Graph:
    """
    Generate a random graph using the Barabasi-Albert model.

    This function creates a scale-free network based on preferential
    attachment.

        Parameters:
            - n (int): The total number of nodes in the network.
            - an (int): The number of initial connected nodes in the network.

        Returns:
            - networkx.Graph: The generated graph.
    """
    G = nx.Graph()
    G.add_nodes_from(range(an))

    # Create new nodes with edges following the preferential attachment
    for new_node in range(an, n):
        sum_denominator = 2 * G.number_of_edges() + G.number_of_nodes()
        # probabilistic list creation
        s = 0
        Lprob = []
        for node in G:
            s += (G.degree[node]+1) / sum_denominator
            Lprob.append(s)
        G.add_node(new_node)

        # new edges determination
        for a in range(an):
            random = rd.random()
            final_node = 0
            while random > Lprob[final_node]:
                final_node += 1
            G.add_edge(new_node, final_node)

    # connectivity condition
    if nx.is_connected(G):
        return G
    else:
        return BarabasiAlbertAlgorithm(n, an)


def meanClustering(G: nx.Graph) -> float:
    """
    Calculate the mean clustering coefficient of a graph.

    Parameters:
    - G (networkx.Graph): The input graph.

    Returns:
    - float: The mean clustering coefficient.
    """
    L = list(nx.clustering(G).items())
    return np.mean([L[i][1] for i in range(G.number_of_nodes())])


def createLogVerificationScaleFree(G: nx.Graph) -> tuple:
    """
    Generate data for log-log verification of a scale-free network.

    Parameters:
    - G (networkx.Graph): The input graph.

    Returns:
    - tuple: Two lists containing log(degree) and log(proportion of nodes
    with that degree).
    """
    res = ([], [])
    dicD = {}
    N = G.number_of_nodes()

    # Calculate the degree distribution
    for node in G:
        d = G.degree[node]
        if d not in dicD:
            dicD[d] = 1
        else:
            dicD[d] += 1

    # Calculate log values for scale-free verification
    for d in dicD:
        res[0].append(np.log(d))
        res[1].append(np.log(dicD[d]/N))

    return res


def adjacenteDiMatriceFromGraph(G: nx.Graph,
                                autoRG: float,
                                duoRG: float) -> tuple:
    """
    Generate a directed graph and adjacency matrix from an undirected graph.

    This function assigns directed edges to an undirected graph based on
    auto-regulation and duo-regulation rates.

    Parameters:
    - G (networkx.Graph): The undirected graph.
    - autoRG (float): The self-regulation rate.
    - duoRG (float): The duo-regulation rate.

    Returns:
    - tuple: A directed graph and its adjacency matrix.
    """
    DiG = nx.DiGraph()
    DiG.add_nodes_from(G)

    # Assign directed edges
    for edge in G.edges():
        rdNumber = rd.random()
        if rdNumber < duoRG:
            DiG.add_edges_from((edge, edge[::-1]), color='black')
        else:
            rdNumber = rd.random()
            if rdNumber < 0.5:
                DiG.add_edges_from([edge], color='blue')
            else:
                DiG.add_edges_from([edge[::-1]], color='blue')

    # Assign self-loops
    for node in G:
        rdNumber = rd.random()
        if rdNumber < autoRG:
            DiG.add_edge(node, node, color='gray')

    # Create adjacency matrix and add activations/inhibitions
    M = nx.to_numpy_array(DiG)
    addActivationInhibition(DiG, M)
    return (DiG, M)


def adjacenteDiMatriceStaredFromGraph(G: nx.Graph,
                                      autoRG: float,
                                      duoRG: float) -> tuple:
    """
    Generate a directed graph and adjacency matrix from an undirected graph
    with a specific starting point.

    This function assigns directed edges to an undirected graph based
    on auto-regulation and duo-regulation rates,
    starting from the node with the highest degree.

    Parameters:
    - G (networkx.Graph): The undirected graph.
    - autoRG (float): The self-regulation rate.
    - duoRG (float): The duo-regulation rate.

    Returns:
    - tuple: A directed graph and its adjacency matrix.
    """
    DiG = nx.DiGraph()
    DiG.add_nodes_from(G)
    degree_dict = dict(G.degree())
    motherNode = max(degree_dict, key=degree_dict.get)
    distance = nx.shortest_path_length(G, motherNode)
    cache = set()

    # Assign directed edges from the most connected node
    for nodeA in distance:
        for nodeB in G[nodeA]:
            edge = (nodeA, nodeB)
            if edge not in cache:
                cache.add(edge)
                cache.add(edge[::-1])
                rdNumber = rd.random()
                if rdNumber < duoRG:
                    DiG.add_edges_from((edge, edge[::-1]), color='black')
                else:
                    DiG.add_edges_from([edge], color='blue')
        # Assign self-loops
        rdNumber = rd.random()
        if rdNumber < autoRG:
            DiG.add_edge(nodeA, nodeA, color='gray')

    # Create adjacency matrix and add activations/inhibitions
    M = nx.to_numpy_array(DiG)
    addActivationInhibition(DiG, M)
    return (DiG, M)


def addActivationInhibition(G: nx.Graph,
                            M: np.ndarray) -> None:
    """
    Adds activation or inhibition labels to the graph edges.

    This function randomly assigns activation or inhibition
    to each edge in the graph.

    Parameters:
    - G (networkx.Graph): The directed graph.
    - M (numpy.ndarray): The adjacency matrix of the graph.
    """
    for u, v in G.edges():
        inhibitionBool = rd.random() < 0.5
        if inhibitionBool:
            M[u][v] *= -1
            G[u][v]['acInColor'] = 'r'  # Red for inhibition
        else:
            G[u][v]['acInColor'] = 'g'  # Green for activation


def addColors(G: nx.Graph,
              M: np.ndarray) -> None:
    """
    Adds colors to the graph edges based on their type.

    This function assigns colors to edges in the graph depending
    on whether they are self-loops,
    mutual connections, or one-directional.

    Parameters:
    - G (networkx.Graph): The directed graph.
    - M (numpy.ndarray): The adjacency matrix of the graph.
    """
    for u, v in G.edges():
        if v == u:
            G[u][v]['color'] = "gray"  # Self-loop
        elif (v, u) in G.edges():
            G[u][v]['color'] = "black"  # Mutual connection
        else:
            G[u][v]['color'] = "blue"  # One-directional connection

        if M[u][v] == 1:
            G[u][v]["acInColor"] = 'g'  # Green for activation
        else:
            G[u][v]["acInColor"] = 'r'  # Red for inhibition

##############################################################################


def main():
    genesNb = 5
    autoRG = 0.05
    duoRG = 0.1
    Ma = np.random.randint(2, size=(genesNb, genesNb))
    for i in [1000]:
        for j in [1, 2, 3]:
            G = BarabasiAlbertAlgorithm(i, j)
            # nx.draw(G, with_labels=True, font_weight='bold')
            # for node in G:
            #    print(node, G.degree[node])
            print(i, j)
            meanClustering(G)
            print()
            coord = createLogVerificationScaleFree(G)
            plt.scatter(coord[0], coord[1], label=f"am = {j}")
            plt.xlabel("log d")
            plt.ylabel("log(P(d)/N)")
            plt.legend()
            plt.title("Scale-free verification")
    plt.savefig("Melvin/Images/scaleFreeProof.png")
    G = BarabasiAlbertAlgorithm(20, 2)
    meanClustering(G)

    # G = nx.DiGraph()
    # G.add_nodes_from(range(genesNb))
    #
    # for i in range(genesNb):
    #    for j in range(genesNb):
    #        if Ma[i,j] :
    #            G.add_edge(i,j)
    #

    # print(list(G.edges()))

    plt.subplot(121)
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.subplot(122)
    DiG = adjacenteDiMatriceStaredFromGraph(G, autoRG, duoRG)[0]
    edges = DiG.edges()
    colors = [DiG[u][v]['color'] for u, v in edges]
    nx.draw_kamada_kawai(DiG, with_labels=True,
                         font_weight='bold', edge_color=colors)
    plt.show()

    K = []
    for i in range(genesNb):
        for j in range(genesNb):
            if Ma[i, j]:
                k = rd.random()
                K.append(k)
##############################################################################


if __name__ == "__main__":
    main()
