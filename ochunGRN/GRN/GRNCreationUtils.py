#!/usr/bin/env python3

import networkx as nx
import random as rd
import matplotlib.pyplot as plt
import numpy as np

# Barabasi-Albert algorithm
# Ma = np.zeros((genesNb,genesNb))


def BarabasiAlbertAlgorithm(n, an):
    G = nx.Graph()
    G.add_nodes_from(range(an))
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


def meanClustering(G):
    L = list(nx.clustering(G).items())
    return np.mean([L[i][1] for i in range(G.number_of_nodes())])


def createLogVerificationScaleFree(G):
    res = ([], [])
    dicD = {}
    N = G.number_of_nodes()
    for node in G:
        d = G.degree[node]
        if d not in dicD:
            dicD[d] = 1
        else:
            dicD[d] += 1
    for d in dicD:
        res[0].append(np.log(d))
        res[1].append(np.log(dicD[d]/N))
    return res


def adjacenteDiMatriceFromGraph(G, autoRG, duoRG):
    DiG = nx.DiGraph()
    DiG.add_nodes_from(G)
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
    for node in G:
        rdNumber = rd.random()
        if rdNumber < autoRG:
            DiG.add_edge(node, node, color='gray')
    M = nx.to_numpy_array(DiG)
    addActivationInhibition(DiG, M)
    return (DiG, M)


def addActivationInhibition(G, M):
    for u, v in G.edges():
        inhibitionBool = rd.random() < 0.5
        if inhibitionBool:
            M[u][v] *= -1
            G[u][v]['acInColor'] = 'r'
        else:
            G[u][v]['acInColor'] = 'g'


def addColors(G, M):
    for u, v in G.edges():
        if v == u:
            G[u][v]['color'] = "gray"
        elif (v, u) in G.edges():
            G[u][v]['color'] = "black"
        else:
            G[u][v]['color'] = "blue"

        if M[u][v] == 1:
            G[u][v]["acInColor"] = 'g'
        else:
            G[u][v]["acInColor"] = 'r'

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
    DiG = adjacenteDiMatriceFromGraph(G, autoRG, duoRG)[0]
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
