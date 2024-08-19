import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from multiprocessing import Pool, cpu_count

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
    
    raise

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
    nbNodes = Graph.number_of_nodes()
    nodes = list(Graph.nodes)
    for u in range(nbNodes-2):
        if u%50==0:
            print(u)
        for v in range(u+1,nbNodes-1):
            for w in range(v+1,nbNodes):
                if issubgraphconnected(Graph, nodes[u],nodes[v],nodes[w]):
                    Subgraph = nx.subgraph(Graph, [nodes[u],nodes[v],nodes[w]])
                    resDic[f"{u}-{v}-{w}"] = findGroupGraph(Subgraph)
    return resDic

def issubgraphconnected(Graph, nodeA, nodeB, nodeC):
    neigboursC = Graph
    i=0
    if nodeA in Graph[nodeC] or nodeC in Graph[nodeA]:
        i+=1
    if nodeA in Graph[nodeB] or nodeB in Graph[nodeA]:
        i+=1
    if nodeB in Graph[nodeC] or nodeC in Graph[nodeB]:
        i+=1
    return i>=2
    
def process_subgraph(args):
    Graph, u, v, w = args
    if issubgraphconnected(Graph, u, v, w):
        Subgraph = Graph.subgraph([u, v, w])
        return f"{u}-{v}-{w}", findGroupGraph(Subgraph)
    return None

def subgraph3N_parallel(Graph, num_processes=4):
    resDic = {}
    nodes = list(Graph.nodes())
    total_combinations = len(nodes) * (len(nodes) - 1) * (len(nodes) - 2) // 6
    
    for i, (u, v, w) in enumerate(combinations(nodes, 3)):
        if i % 1000000 == 0:
            print(f"Processed {i}/{total_combinations} combinations ({i/total_combinations*100:.2f}%)")
        
        if issubgraphconnected(Graph, u, v, w):
            Subgraph = Graph.subgraph([u, v, w])
            resDic[f"{u}-{v}-{w}"] = findGroupGraph(Subgraph)
    
    return resDic

def main():
    G = nx.scale_free_graph(1000)
    print(subgraph3N(G))
    nx.draw_circular(G,with_labels=True)
    plt.show()
    TGraph = [nx.DiGraph() for i in range(13)]
    TGraph[0].add_edges_from([(0,1),(1,2),(2,0)])
    TGraph[1].add_edges_from([(1,0),(1,2),(2,0)])
    TGraph[2].add_edges_from([(0,1),(1,2)])
    TGraph[3].add_edges_from([(1,0),(1,2)])
    TGraph[4].add_edges_from([(0,1),(2,1)])
    TGraph[5].add_edges_from([(0,1),(1,2),(2,0),(0,2)])
    TGraph[6].add_edges_from([(1,0),(1,2),(2,0),(0,2)])
    TGraph[7].add_edges_from([(0,1),(2,1),(2,0),(0,2)])
    TGraph[8].add_edges_from([(0,1),(1,2),(2,1)])
    TGraph[9].add_edges_from([(0,1),(1,2),(1,0)])
    TGraph[10].add_edges_from([(0,1),(1,2),(1,0),(2,1)])
    TGraph[11].add_edges_from([(0,1),(1,2),(1,0),(2,1),(2,0)])
    TGraph[12].add_edges_from([(0,1),(1,2),(1,0),(2,1),(2,0),(0,2)])
    for i in range(3):
        plt.subplot(3,5,i+2)
        nx.draw_circular(TGraph[i],arrowsize=30)
        plt.title(findGroupGraph(TGraph[i]))
    for i in range(5):
        plt.subplot(3,5,i+6)
        nx.draw_circular(TGraph[i+3],arrowsize=30)
        plt.title(findGroupGraph(TGraph[i+3])) 
    for i in range(5):
        plt.subplot(3,5,i+11)
        nx.draw_circular(TGraph[i+8],arrowsize=30)
        plt.title(findGroupGraph(TGraph[i+8]))
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.savefig("Melvin/Images/allGroups.png")


if __name__ == "__main__":
    main()
        
