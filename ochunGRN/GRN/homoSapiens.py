#!/usr/bin/env python3

import networkx as nx
import numpy as np
import pandas as pd
import pkg_resources
import matplotlib.pyplot as plt
from ochunGRN.GRN.genesGroup import subgraph3N

# Load the data from an Excel file
if __name__ == "__main__":
    # If running as a script, load the Excel file directly
    document = pd.read_excel("ochunGRN/GRN/41598_2021_3625_MOESM5_ESM.xlsx")
else:
    # If imported as a module, use pkg_resources to locate the file
    excel_path = pkg_resources.resource_filename(
        "ochunGRN", "GRN/41598_2021_3625_MOESM5_ESM.xlsx")
    document = pd.read_excel(excel_path)

# Load the data for different species from the Excel file
arabidopsisThaliana = (document["""Supplementary Table S1: Networks. A spreadsheet file with filtered networks"""].tolist()[2:], document["Unnamed: 1"].tolist()[2:], "Arabidopsis thaliana")  # noqa : E501
drosophilaMelanogaster = (document["Unnamed: 2"].tolist()[2:], document["Unnamed: 3"].tolist()[2:], "Drosophila Melanogaster")  # noqa : E501
escherichniaColi = (document["Unnamed: 4"].tolist()[2:], document["Unnamed: 5"].tolist()[2:], "Escherichnia coli")  # noqa : E501
homoSapiens = (document["Unnamed: 6"].tolist()[2:], document["Unnamed: 7"].tolist()[2:], "Homo sapiens")  # noqa : E501
saccharomycesCerevisiae = (document["Unnamed: 8"].tolist()[2:], document["Unnamed: 9"].tolist()[2:], "Saccharomyces cerevisiae")  # noqa : E501


def remove_nan_values(species: list) -> None:
    """
    Removes NaN values from the end of the gene and interaction lists of
    a species.

        Args:
            - species (list): A list containing two lists:
            [gene_list, interaction_list]
    """
    # Check if the last element is not a string (i.e., it's NaN)
    while type(species[0][-1]) is not str:
        species[0].pop()  # Remove the last gene entry
        species[1].pop()  # Remove the corresponding interaction entry


def create_graph(species: list) -> nx.DiGraph:
    """
    Creates a directed graph from gene interactions for a given species.

        Args:
            - species (list): A list containing two lists:
            [gene_list, interaction_list]

        Returns:
            nx.DiGraph: A directed graph representing gene interactions.
    """
    N = len(species[0])  # Number of gene interactions
    Graph = nx.DiGraph()  # Create an empty directed graph
    for i in range(N):
        Graph.add_edge(species[0][i], species[1][i])  # Add edges to the graph
    return Graph


def find_autoregulatory_genes(species: list) -> tuple:
    """
    Finds genes that regulate themselves (autoregulatory genes) in a given
    species.

        Args:
            - species (list): A list containing two lists:
            [gene_list, interaction_list]

        Returns:
            tuple: A tuple containing the list of autoregulatory genes
            and their proportion.
    """
    autoregulatory_genes = []
    N = len(species[0])  # Number of gene interactions
    for i in range(N):
        if species[0][i] == species[1][i]:  # Check if a gene regulates itself
            autoregulatory_genes.append(species[0][i])
    # Return the genes and their proportion
    return (autoregulatory_genes, len(autoregulatory_genes) / len(species[0]))


def find_double_regulatory_genes(graph: nx.DiGraph) -> tuple:
    """
    Finds pairs of genes that regulate each other in a given species graph.

        Args:
            - graph (nx.DiGraph): A directed graph representing gene
            interactions.

        Returns:
            tuple: A tuple containing the list of double-regulatory gene pairs
            and their proportion.
    """
    double_regulation_genes = []
    for nodeA in graph:  # Iterate over all nodes (genes)
        # Iterate over successors of each node (genes regulated by nodeA)
        for nodeB in graph.successors(nodeA):
            # Check if nodeA is regulated by nodeB
            if nodeA in graph.successors(nodeB):
                double_regulation_genes.append((nodeA, nodeB))
    # Return the pairs and their proportion
    return (double_regulation_genes,
            len(double_regulation_genes) / graph.number_of_edges())


def calculate_ffl_ratio(graph: nx.DiGraph) -> None:
    """
    Calculates and prints the ratio of feed-forward loops (FFL) in a graph.

        Args:
            - graph (nx.DiGraph): A directed graph representing gene
            interactions.
    """
    print(graph.number_of_nodes(), graph.number_of_edges())
    print(nx.is_weakly_connected(graph))
    dict_ffl = subgraph3N(graph)
    motifs = {v: 0 for v in dict_ffl.values()}  # Initialize motif counts
    N = len(dict_ffl)
    print(N)
    for i in dict_ffl.keys():
        motifs[dict_ffl[i]] += 1
    # Print sorted motifs by count
    print(sorted(motifs.items(), key=lambda item: item[1]))
    print(motifs['FFL'] / N)  # Print the ratio of FFL motifs


def plot_autoregulatory_distribution(graph: nx.DiGraph,
                                     group_count: int) -> None:
    """
    Plots the distribution of autoregulatory genes based on degree groups.

    Args:
    graph (nx.DiGraph): A directed graph representing gene interactions.
    group_count (int): The number of groups to divide the nodes by degree.
    """
    autoregulatory_dict = {}
    total_dict = {}
    degree_list = list(graph.degree())
    degree_values = [degree_list[i][1] for i in range(len(degree_list))]
    degree_values.sort()
    gene_count = graph.number_of_nodes()
    # Create degree group pivot points
    pivot_points = [degree_values[int(i * gene_count / (group_count + 1))]
                    for i in range(1, group_count)]
    print(pivot_points)

    for node in graph:
        group_val = 0
        deg = graph.degree(node)
        while group_val < group_count - 1 and deg > pivot_points[group_val]:
            group_val += 1
        if node in graph.successors(node):  # If a node is autoregulatory
            if group_val in autoregulatory_dict:
                autoregulatory_dict[group_val] += 1
            else:
                autoregulatory_dict[group_val] = 1
        if group_val in total_dict:
            total_dict[group_val] += 1
        else:
            total_dict[group_val] = 1

    sorted_autoregulatory = sorted(autoregulatory_dict.items())
    print(sorted_autoregulatory)

    X, Y = [], []
    for item in sorted_autoregulatory:
        deg = item[0]
        X.append(str(item[0]))
        # Normalize by total number of nodes in each group
        Y.append(item[1] / total_dict[deg])
    plt.bar(X, Y)
    plt.show()


def find_high_degree_nodes(graph: nx.DiGraph, min_degree: int = 10) -> list:
    """
    Finds nodes (genes) in the graph with a degree greater than a specified
    minimum.

        Args:
            - graph (nx.DiGraph): A directed graph representing gene
            interactions.
            - min_degree (int): The minimum degree threshold.

        Returns:
            list: A list of nodes with a degree higher than the specified
            threshold.
    """
    return [node for node in graph if graph.degree(node) > min_degree]


def main():
    for species in [arabidopsisThaliana,
                    drosophilaMelanogaster,
                    escherichniaColi,
                    homoSapiens,
                    saccharomycesCerevisiae]:
        remove_nan_values(species)
        graph = create_graph(species)
        print(f"\nSpecies: {species[2]}")
        print("Nodes:", graph.number_of_nodes())
        print("Edges:", graph.number_of_edges())
        print("Autoregulatory Genes:", find_autoregulatory_genes(species))
        print("Double Regulatory Genes:", find_double_regulatory_genes(graph))
        clustering_coefficients = list(nx.clustering(graph).items())
        print("Average Clustering Coefficient:",
              np.mean([clustering_coefficients[i][1]
                       for i in range(graph.number_of_nodes())]))
        print("Autoregulatory Distribution:",
              plot_autoregulatory_distribution(graph, 100))
        print("Nodes with Degree >= 10:", find_high_degree_nodes(graph))


if __name__ == "__main__":
    main()
