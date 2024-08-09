from GRNCreationUtils import meanClustering, BarabasiAlbertAlgorithm, adjacenteDiMatriceFromGraph, addColors
from genesGroupe import subgraph3N
import networkx as nx

def randomGrn(genesNb:int, autoRG:float, duoRG:float, resDict:dict=None):
    if resDict is None:
        resDict = {}
    Graph = BarabasiAlbertAlgorithm(genesNb, 2)
    Graph,M = adjacenteDiMatriceFromGraph(Graph, autoRG, duoRG)
    resDict["Graph"] = Graph
    resDict["genesNb"], resDict["autoRG"], resDict["duoRG"] = genesNb, autoRG, duoRG
    resDict["AdjMatrice"] = M
    resDict["meanClustering"] = meanClustering(Graph)
    resDict["subGraph"] = subgraph3N(Graph)
    
    return resDict

def GrnFromAdj(AdjMatrice, resDict:dict=None):
    if resDict is None:
        resDict = {}
    genesNb = len(AdjMatrice)
    Graph = nx.from_numpy_array(AdjMatrice, create_using=nx.DiGraph)
    autoRG = 0
    for i in range(genesNb):
        if AdjMatrice[i][i]!=0:
            autoRG+=1
    N = Graph.number_of_edges()
    N-=autoRG
    duoRG = 0
    for i in range(genesNb):
        for j in range(i+1,genesNb):
            if AdjMatrice[i][j]!=0 and AdjMatrice[j][i]!=0:
                N-=1
                duoRG+=1
    autoRG/=Graph.number_of_edges()
    duoRG/=N
    addColors(Graph,AdjMatrice)
    resDict["Graph"] = Graph
    resDict["genesNb"], resDict["autoRG"], resDict["duoRG"] = genesNb, autoRG, duoRG
    resDict["AdjMatrice"] = AdjMatrice
    resDict["meanClustering"] = meanClustering(Graph)
    resDict["subGraph"] = subgraph3N(Graph)
    return resDict
