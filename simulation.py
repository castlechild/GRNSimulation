#!/usr/bin/env python3
from scipy.integrate import solve_ivp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from GRN.GRNCreationUtils import meanClustering, BarabasiAlbertAlgorithm, adjacenteDiMatriceFromGraph
from GRN.genesGroupe import subgraph3N
from ODESystems.coefficientFinder import getCoefficient
from ODESystems.massAction import MaMatrice, massAction
from ODESystems.Hill import HillEquation
from ODESystems.indirect import indirect


def simulation(ODEs:list, T:tuple, genesNb:int=None, autoRG:float=None, duoRG:float=None, nHill=1, Graph=None, M=None, Coeff:dict=None, plot:bool=False, saveName:str=None) -> dict :
    """
    Simulates the given ODE systems over the specified time interval.
    
    Args:
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
    - dictionary containing the simulation results.
    """

    def otherODE(L):
        for ode in L:
            if ode not in ["massAction", "Hill","indirect"]:
                return True
        return False 
    
    if len(ODEs) == 0 and otherODE(ODEs): 
        raise ValueError("ODEs not valid")
    if (Graph is None and M is None) and ((genesNb is None) or (autoRG is None) or (duoRG is None)):
        raise ValueError("the generation did not have a graph entry or parameter to create one")
    if len(T) != 2 and T[0] >= T[1]:
        raise ValueError("the temporal pair is invalid")
    
    if Graph is None:
        if M is None:    
            Graph = BarabasiAlbertAlgorithm(genesNb, 2)
            Graph,M = adjacenteDiMatriceFromGraph(Graph, autoRG, duoRG)
        else :
            genesNb = len(M)
            Graph = nx.DiGraph()
            for i in range(genesNb):
                for j in range(genesNb):
                    if M[i][j] != 0:
                        Graph.add_edge(i, j)
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
        Prout = [Coeff["TranslationsRate"][i]/Coeff["mRNAAvg"][i] for i in range(genesNb)]
        K = np.resize(Prout,(genesNb,genesNb))
        Ma = MaMatrice(M, K)
        
        equation = lambda t,G: massAction(t, G, Ma)     
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        resDict["massActionY"] = solution.y
        resDict["massActionX"] = solution.t
    
    if "Hill" in ODEs:
        K = Coeff["TranscriptionsRate"]
        Kdeg = Coeff["mRNAsDeg"]
        equation = lambda t,G: HillEquation(t, G, M, K, G0, [0]*genesNb, Kdeg, 5)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        resDict["HillsY"] = solution.y
        resDict["HillsX"] = solution.t

    if "indirect" in ODEs:
        k_P = Coeff["TranslationsRate"]
        Ka_P = G0
        k_mRNA = Coeff["TranscriptionsRate"]
        K_degP = Coeff["ProtsDeg"]
        K_degMRNA = Coeff["mRNAsDeg"]
        def equation(t, G):
            mRNA = G[:genesNb]
            P = G[genesNb:]
            return indirect(t, mRNA, P, M, k_mRNA, k_P, Ka_P, K_degP, K_degMRNA, 5)
        G0_indirect = np.concatenate((G0, Coeff["ProtAvg"]))
        solution = solve_ivp(equation, [t0, tf], G0_indirect, max_step= 0.5)
        resDict["indirectY"] = solution.y[:genesNb]
        resDict["indirectX"] = solution.t

    if plot:
        plt.figure()
        edges = Graph.edges()
        colors = [Graph[u][v]['color'] for u,v in edges]
        font = {'family':'serif','color':'darkred','size':8}
        plt.subplot(2, 2, 1)
        nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=colors)
        plt.title("Graph RGG")
        
        if "Hill" or "indirect" in ODEs :
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
            plt.subplot(len(ODEs),2,min(2*len(ODEs),4))
            for solGenes in range(genesNb):
                plt.plot(resDict["HillsX"], resDict["HillsY"][solGenes], label=solGenes)
            plt.xlabel("time (h)", fontdict=font)
            plt.ylabel("Genes concentrations", fontdict=font)
            plt.title("Hill law Simulation")
            plt.legend()

        if "indirect" in ODEs:
            plt.subplot(len(ODEs), 2, len(ODEs)*2)
            for solGenes in range(genesNb):
                plt.plot(resDict["indirectX"], resDict["indirectY"][solGenes], label=solGenes)
            plt.xlabel("time (h)", fontdict=font)
            plt.ylabel("mRNA concentrations", fontdict=font)
            plt.title("indirect (massAction & Hill laws) Simulation")
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
    simulation(["massAction","Hill","indirect"],(0,100),NB_GENES,AUTO_RG,DUO_RG,plot=True)


if __name__ == "__main__":
    main()