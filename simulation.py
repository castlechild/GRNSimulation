#!/usr/bin/env python3
from scipy.integrate import solve_ivp
from reseauGenes.RGGCreation import *
from systemeEDO.massAction import *
from systemeEDO.Hills import *


def simulation(ODEs, T, genesNb=None, autoRG=None, duoRG=None, Graph=None, plot=False, saveName=None):

    def otherODE(L):
        for ode in L:
            if ode not in ["massAction", "Hills"]:
                return True
        return False 
    if len(ODEs) == 0 and otherODE(ODEs): 
        raise ValueError("Les ODEs ne sont pas valides")
    if Graph is None and ((genesNb is None) or (autoRG is None) or (duoRG is None)):
        raise ValueError("la génération n'a pas eu d'entrée de graphe ni de paramètre pour en créé un")
    if len(T) != 2 and T[0] >= T[1]:
        raise ValueError("le couple temporel est invalide")
    
    if Graph is None:
        Graph = BarabasiAlbertAlgorithm(genesNb, 2)
        Graph,M = adjacenteDiMatriceFromGraph(Graph, autoRG, duoRG)
    else:
        M = np.transpose(nx.to_numpy_array(Graph))
        genesNb = Graph.number_of_nodes()
        
    res = {}
    t0, tf = T
    G0 = [2]*genesNb
    kk = [1]*genesNb
    res["Graph"] = Graph
    res["AdjMatrice"] = M
    res["genesNb"], res["autoRG"], res["duoRG"] = genesNb, autoRG, duoRG
    res["meanClustering"] = meanClustering(Graph)

    if "massAction" in ODEs:
        K = np.random.random((genesNb, genesNb))
        Ma = MaMatrice(M, K)
        
        equation = lambda t,G: masseAction(t, G, Ma)     
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.05)
        res["massActionY"] = solution.y
        res["massActionX"] = solution.t
    
    if "Hills" in ODEs:
        equation = lambda t,G: HillEquation(t, G, M, kk*2, kk, kk, kk, 1)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        res["HillsY"] = solution.y
        res["HillsX"] = solution.t

    if plot:
        plt.clf()
        edges = Graph.edges()
        colors = [Graph[u][v]['color'] for u,v in edges]
        font = {'family':'serif','color':'darkred','size':10}
        plt.subplot(2,2,1)
        nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=colors)
        plt.title("Graph RGG")
        
        acInColors = [Graph[u][v]['acInColor'] for u,v in edges]
        plt.subplot(2,2,3)
        nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=acInColors, connectionstyle="arc3,rad=0.05" )
        #pos = nx.circular_layout(Graph)
        #nx.draw_networkx_nodes(Graph, pos, label= range(genesNb))
        #nx.draw_networkx_edges(Graph, pos,connectionstyle="arc3,rad=0.1",edge_color= acInColors)

        plt.title("Graph RGG Activator/Inhibitor")

        if "massAction" in ODEs:
            print(len(ODEs))
            plt.subplot(len(ODEs),2,2)
            for solGenes in range(genesNb):
                plt.plot(res["massActionX"], res["massActionY"][solGenes], label=solGenes)
            plt.xlabel("time", fontdict=font)
            plt.ylabel("Genes concentrations", fontdict=font)
            plt.title("Mass Action law simulation")
            plt.legend()

        if "Hills" in ODEs:
            plt.subplot(len(ODEs),2,2*len(ODEs))
            for solGenes in range(genesNb):
                plt.plot(res["HillsX"], res["HillsY"][solGenes], label=solGenes)
            plt.xlabel("time", fontdict=font)
            plt.ylabel("Genes concentrations", fontdict=font)
            plt.title("Hills law Simulation")
            plt.legend()
        manager = plt.get_current_fig_manager()
        manager.full_screen_toggle()
        if saveName is None:     
            plt.savefig("simulation.png")
        else :
            plt.savefig(saveName)

    return res
  

def main():
    NB_GENES = 7
    AUTO_RG = 0.1
    DUO_RG = 0.2
    simulation(["massAction","Hills"],(0,20),NB_GENES,AUTO_RG,DUO_RG,plot=True)


if __name__ == "__main__":
    main()