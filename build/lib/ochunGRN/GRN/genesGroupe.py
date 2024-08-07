import networkx as nx
import matplotlib.pyplot as plt

def findGroupGraph(threeNGraph):
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
    resDic = {}
    nbNodes = Graph.number_of_nodes()
    for u in range(nbNodes-2):
        for v in range(u+1,nbNodes-1):
            for w in range(v+1,nbNodes):
                Subgraph = nx.subgraph(Graph, [u,v,w])
                UndirectG = Subgraph.to_undirected()
                if nx.is_connected(UndirectG):
                    resDic[f"{u}-{v}-{w}"] = findGroupGraph(Subgraph)
    return resDic


def main():
    G = nx.scale_free_graph(10)
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
        
