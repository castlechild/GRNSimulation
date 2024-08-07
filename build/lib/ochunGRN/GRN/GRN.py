from GRNCreationUtils import meanClustering, BarabasiAlbertAlgorithm, adjacenteDiMatriceFromGraph
from genesGroupe import subgraph3N

def Grn(genesNb:int, autoRG:float, duoRG:float, resDict:dict=None):
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
