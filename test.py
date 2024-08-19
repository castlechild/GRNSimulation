import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations,islice
import multiprocessing
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from tqdm import tqdm

document = pd.read_excel("ochunGRN/GRN/41598_2021_3625_MOESM5_ESM.xlsx")

arabidopsisThaliana=(document["Supplementary Table S1: Networks. A spreadsheet file with filtered networks"].tolist()[2:],document["Unnamed: 1"].tolist()[2:], "Arabidopsis thaliana" )
drosophilaMelanogaster=(document["Unnamed: 2"].tolist()[2:],document["Unnamed: 3"].tolist()[2:],"Drosophila Melanogaster")
escherichniaColi=(document["Unnamed: 4"].tolist()[2:],document["Unnamed: 5"].tolist()[2:], "Escherichnia coli")
homoSapiens=(document["Unnamed: 6"].tolist()[2:],document["Unnamed: 7"].tolist()[2:], "Homo sapiens")
saccharomycesCerevisiae=(document["Unnamed: 8"].tolist()[2:],document["Unnamed: 9"].tolist()[2:], "Saccharomyces cerevisiae")

def findGroupGraph(threeNGraph):
    """
    Identify the type of three-node subgraph (pattern) in a directed graph.

    This function takes a three-node subgraph and classifies it based on its structural patterns into one of the following types: "Fan-In", "Fan-Out", "Cascade", "FFL", "Triangular", "C-D", "C-E", "SCascade", "dC", "dD", "dE", "5edges", and "6edges".

    Parameters:
    - threeNGraph (networkx.DiGraph): The three-node subgraph to analyze.

    Returns:
    - str: The detected pattern type among the available options.
    """
    threeNGraphWAuto = nx.DiGraph()
    edges = threeNGraph.edges()
    for u,v in edges:
        if u != v:
            threeNGraphWAuto.add_edge(u,v)
    
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
                if u!=v and (u,v) not in edges and (v,u) not in edges:
                    if threeNGraphWAuto.in_degree[u] == 1 and threeNGraphWAuto.in_degree[v] == 1:
                        return "C-D"
                    else :
                        return "C-E"
        for u in threeNGraphWAuto:
            if threeNGraphWAuto.in_degree[u] == 2:
                return "FFL"
        return "Triangular"
    
    elif nEdges == 4:
        for u in threeNGraphWAuto:
            for v in threeNGraphWAuto:
                if u!=v and (u,v) not in edges and (v,u) not in edges:
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

def subgraph3N(Graph):
    """
    Extract and analyze all three-node subgraphs from a given graph.

    This function traverses a given graph, extracts all subgraphs containing exactly three nodes, and uses the `findGroupGraph` function to classify each one.

    Parameters:
    - Graph (networkx.Graph): The graph to analyze.

    Returns:
    - dict: A dictionary where the keys are the types of patterns found and the values are the counts of each type.
    """
    resDic = {}
    nodes = list(Graph.nodes)
    total_combinations = len(nodes) * (len(nodes) - 1) * (len(nodes) - 2) // 6
    with tqdm(total=total_combinations) as pbar:
        for i, (u, v, w) in enumerate(combinations(nodes, 3)):
            if issubgraphconnected(Graph, u,v,w):
                Subgraph = nx.subgraph(Graph, [u,v,w])
                resDic[f"{u}-{v}-{w}"] = findGroupGraph(Subgraph)
            pbar.update(1)
    return resDic

def issubgraphconnected(Graph, nodeA, nodeB, nodeC):
    i=0
    if nodeA in Graph[nodeC] or nodeC in Graph[nodeA]:
        i+=1
    if nodeA in Graph[nodeB] or nodeB in Graph[nodeA]:
        i+=1
    if nodeB in Graph[nodeC] or nodeC in Graph[nodeB]:
        i+=1
    return i>=2

def process_combination(combo,Graph,nodes):
    u, v, w = combo
    if issubgraphconnected(Graph, nodes[u], nodes[v], nodes[w]):
        Subgraph = nx.subgraph_view(Graph, filter_node=lambda n: n in [nodes[u], nodes[v], nodes[w]])
        return f"{u}-{v}-{w}", findGroupGraph(Subgraph)
    return None    
    
def subgraph3N_parallel(Graph):
    print("Début de subgraph3N_parallel")
    resDic = {}
    nodes = list(Graph.nodes)
    print(f"Nombre de nœuds : {len(nodes)}")
    total_combinations = len(nodes) * (len(nodes) - 1) * (len(nodes) - 2) // 6
    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        print("Début de l'exécution parallèle")
        batch_size = 1000
        combo_iter = combinations(range(len(nodes)), 3)
        with tqdm(total=total_combinations) as pbar:
        
            while True:
                batch = list(islice(combo_iter, batch_size))
                if not batch:
                    break
                
                for result in enumerate(executor.map(lambda x: process_combination(x, Graph, nodes), batch)):
                    if result:
                        key, value = result
                        resDic[key] = value
                    pbar.update(1)
    
    print("Fin de subgraph3N_parallel")
    return resDic

def FFLratio(Gspecies):
    print(Gspecies.number_of_nodes(),Gspecies.number_of_edges())
    print(nx.is_weakly_connected(Gspecies))
    DictFFL = subgraph3N(Gspecies)
    motifs = {v:0 for v in DictFFL.values()}
    N=len(DictFFL)
    print(N)
    for i in DictFFL.keys():
        motifs[DictFFL[i]]+=1
    print(sorted(motifs.items(), key=lambda item: item[1]))
    print(motifs['FFL']/N)

def graphCreator(species):
    N = len(species[0])
    Graph = nx.DiGraph()
    for i in range(N):
        Graph.add_edge(species[0][i],species[1][i])
    return Graph

def main():
    print("ratioFFL", FFLratio(graphCreator(homoSapiens)))


if __name__ == "__main__":
    main()