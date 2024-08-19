#!/usr/bin/env python3

import networkx as nx
import numpy as np
import pandas as pd
import pkg_resources
from genesGroupe import subgraph3N_parallel
if __name__ == "__main__":
    document = pd.read_excel("ochunGRN/GRN/41598_2021_3625_MOESM5_ESM.xlsx")
else:
    excel_path = pkg_resources.resource_filename(
        "ochunGRN", "GRN/41598_2021_3625_MOESM5_ESM.xlsx")
    document = pd.read_excel(excel_path)

arabidopsisThaliana = (
    document["""Supplementary Table S1: Networks.A spreadsheet
             file with filtered networks"""].tolist()[2:],
    document["Unnamed: 1"].tolist()[2:], "Arabidopsis thaliana")
drosophilaMelanogaster = (document["Unnamed: 2"].tolist()[2:],
                          document["Unnamed: 3"].tolist()[2:],
                          "Drosophila Melanogaster")
escherichniaColi = (document["Unnamed: 4"].tolist()[2:],
                    document["Unnamed: 5"].tolist()[2:], "Escherichnia coli")
homoSapiens = (document["Unnamed: 6"].tolist()[2:],
               document["Unnamed: 7"].tolist()[2:], "Homo sapiens")
saccharomycesCerevisiae = (document["Unnamed: 8"].tolist()[2:],
                           document["Unnamed: 9"].tolist()[2:],
                           "Saccharomyces cerevisiae")


def nanDestructor(species):
    while type(species[0][-1]) is not str:
        species[0].pop()
        species[1].pop()


def graphCreator(species):
    N = len(species[0])
    Graph = nx.DiGraph()
    for i in range(N):
        Graph.add_edge(species[0][i], species[1][i])
    return Graph


def findAutoRegulationGenes(species):
    AutoRegulationGenes = []
    N = len(species[0])
    for i in range(N):
        if species[0][i] == species[1][i]:
            AutoRegulationGenes.append(species[0][i])
    return (AutoRegulationGenes, len(AutoRegulationGenes)/len(species[0]))


def findDoubleRegulationGenes(GSpecies):
    DoubleRegulationGenes = []
    for nodeA in GSpecies:
        for nodeB in GSpecies.successors(nodeA):
            if nodeA in GSpecies.successors(nodeB):
                DoubleRegulationGenes.append((nodeA, nodeB))
    return (DoubleRegulationGenes,
            len(DoubleRegulationGenes)/GSpecies.number_of_edges())


def FFLratio(Gspecies):
    print(Gspecies.number_of_nodes(), Gspecies.number_of_edges())
    print(nx.is_weakly_connected(Gspecies))
    DictFFL = subgraph3N_parallel(Gspecies)
    motifs = {v: 0 for v in DictFFL.values()}
    N = len(DictFFL)
    print(N)
    for i in DictFFL.keys():
        motifs[DictFFL[i]] += 1
    print(sorted(motifs.items(), key=lambda item: item[1]))
    print(motifs['FFL']/N)


def main():
    for species in [arabidopsisThaliana, drosophilaMelanogaster,
                    escherichniaColi, homoSapiens, saccharomycesCerevisiae]:
        nanDestructor(species)
        GSpecies = graphCreator(species)
        print(f"\nspecies : {species[2]}")
        print("nodes :", GSpecies.number_of_nodes())
        print("edges :", GSpecies.number_of_edges())
        print("Auto Rg :", findAutoRegulationGenes(species))
        print("Double Rg:", findDoubleRegulationGenes(GSpecies))
        L = list(nx.clustering(GSpecies).items())
        print("moyenne coefficient de clustering", np.mean(
            [L[i][1] for i in range(GSpecies.number_of_nodes())]))
    print("ratioFFL", FFLratio(graphCreator(homoSapiens)))

# options = {
#     'node_color': 'black',
#     'node_size': 0.4,
#     'width': 0.02,
# }
# G = graphCreator(saccharomycesCerevisiae)
# subax1 = plt.subplot(111)
# nx.draw(G, **options)
# plt.show()


if __name__ == "__main__":
    main()
