#!/usr/bin/env python3
from scipy.integrate import solve_ivp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from GRN.GRNCreation import meanClustering, BarabasiAlbertAlgorithm, adjacenteDiMatriceFromGraph
from ODESystems.coefficientFinder import getCoefficient
from GRN.genesGroupe import subgraph3N
from ODESystems.massAction import MaMatrice, massAction
from ODESystems.Hill import HillEquation


def simulation(ODEs:list, T:tuple, genesNb:int=None, autoRG:float=None, duoRG:float=None, Graph=None, Coeff:dict=None, plot:bool=False, saveName:str=None) -> dict :
    """
    Simulates the given ODE systems over the specified time interval.
    
    Parameters:
    - odes: list of ODE systems to simulate.
    - T: tuple representing the time interval.
    - genes_nb: number of genes.
    - auto_rg: auto-regulation parameter.
    - duo_rg: dual regulation parameter.
    - graph: graph structure.
    - coeff: coefficients dictionary.
    - plot: boolean indicating whether to plot the results.
    - save_name: filename to save the plot.
    
    Returns:
    - res_dict: dictionary containing the simulation results.
    """

    def otherODE(L):
        for ode in L:
            if ode not in ["massAction", "Hill"]:
                return True
        return False 
    
    if len(ODEs) == 0 and otherODE(ODEs): 
        raise ValueError("ODEs not valid")
    if Graph is None and ((genesNb is None) or (autoRG is None) or (duoRG is None)):
        raise ValueError("the generation did not have a graph entry or parameter to create one")
    if len(T) != 2 and T[0] >= T[1]:
        raise ValueError("the temporal pair is invalid")
    
    if Graph is None:
        Graph = BarabasiAlbertAlgorithm(genesNb, 2)
        Graph,M = adjacenteDiMatriceFromGraph(Graph, autoRG, duoRG)
    else:
        M = np.transpose(nx.to_numpy_array(Graph))
        genesNb = Graph.number_of_nodes()
    
    if Coeff is None:
        Coeff = getCoefficient(genesNb)
    
        
    resDict = {}
    t0, tf = T
    G0 = Coeff["mRNAAvg"]
    kk = [1]*genesNb
    resDict["Graph"] = Graph
    resDict["AdjMatrice"] = M
    resDict["Coefficients"] = Coeff
    resDict["genesNb"], resDict["autoRG"], resDict["duoRG"] = genesNb, autoRG, duoRG
    resDict["meanClustering"] = meanClustering(Graph)
    resDict["subGraph"] = subgraph3N(Graph)

    if "massAction" in ODEs:
        K = np.resize(Coeff["TranslationsRate"],(genesNb,genesNb))
        Ma = MaMatrice(M, K)
        
        equation = lambda t,G: massAction(t, G, Ma)     
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        resDict["massActionY"] = solution.y
        resDict["massActionX"] = solution.t
    
    if "Hill" in ODEs:
        K = Coeff["TranslationsRate"]
        Kdeg = np.sqrt(2)/Coeff["ProtsHalfTime"]
        equation = lambda t,G: HillEquation(t, G, M, K, G0, [0]*genesNb, Kdeg, 1)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        resDict["HillsY"] = solution.y
        resDict["HillsX"] = solution.t

    if plot:
        plt.figure()
        edges = Graph.edges()
        colors = [Graph[u][v]['color'] for u,v in edges]
        font = {'family':'serif','color':'darkred','size':10}
        plt.subplot(1 + int("Hill" in ODEs), 2, 1)
        nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=colors)
        plt.title("Graph RGG")
        
        if "Hill" in ODEs :
            acInColors = [Graph[u][v]['acInColor'] for u,v in edges]
            plt.subplot(2,2,3)
            nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=acInColors, connectionstyle="arc3,rad=0.05" )
            #pos = nx.circular_layout(Graph)
            #nx.draw_networkx_nodes(Graph, pos, label= range(genesNb))
            #nx.draw_networkx_edges(Graph, pos,connectionstyle="arc3,rad=0.1",edge_color= acInColors)

            plt.title("Graph RGG Activator/Inhibitor")

        if "massAction" in ODEs:
            plt.subplot(len(ODEs),2,2)
            for solGenes in range(genesNb):
                plt.plot(resDict["massActionX"], resDict["massActionY"][solGenes], label=solGenes)
            plt.xlabel("time (h)", fontdict=font)
            plt.ylabel("Genes concentrations", fontdict=font)
            plt.title("Mass Action law simulation")
            plt.legend()

        if "Hill" in ODEs:
            plt.subplot(len(ODEs),2,2*len(ODEs))
            for solGenes in range(genesNb):
                plt.plot(resDict["HillsX"], resDict["HillsY"][solGenes], label=solGenes)
            plt.xlabel("time (h)", fontdict=font)
            plt.ylabel("Genes concentrations", fontdict=font)
            plt.title("Hill law Simulation")
            plt.legend()
        manager = plt.get_current_fig_manager()
        manager.full_screen_toggle()
        if saveName is None:     
            plt.savefig("simulation.png")
        else :
            plt.savefig(saveName)

    return resDict
  

def main():
    NB_GENES = 7
    AUTO_RG = 0.1
    DUO_RG = 0.2
    simulation(["massAction","Hill"],(0,0.2),NB_GENES,AUTO_RG,DUO_RG,plot=True)


if __name__ == "__main__":
    main()