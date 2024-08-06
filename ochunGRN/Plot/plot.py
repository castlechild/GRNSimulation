import networkx as nx
import matplotlib.pyplot as plt

def plotGraph(GenesDict, saveName=None):
    Graph = GenesDict["Graph"]
    ODEs = GenesDict["ODEs"]
    plt.figure()
    edges = Graph.edges()
    colors = [Graph[u][v]['color'] for u,v in edges]
    
    plt.subplot(2, 1, 1)
    nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=colors)
    plt.title("Graph RGG")
    
    if "Hill" or "indirect" in ODEs :
        acInColors = [Graph[u][v]['acInColor'] for u,v in edges]
        plt.subplot(2,1,2)
        nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=acInColors, connectionstyle="arc3,rad=0.05" )
        #pos = nx.circular_layout(Graph)
        #nx.draw_networkx_nodes(Graph, pos, label= range(genesNb))
        #nx.draw_networkx_edges(Graph, pos,connectionstyle="arc3,rad=0.1",edge_color= acInColors)
        plt.title("Graph RGG Activator/Inhibitor")
    if saveName is not None:     
        plt.savefig(saveName)

def plotSim(GenesDict, saveName=None):
    font = {'family':'serif','color':'darkred','size':8}
    ODEs = GenesDict["ODEs"]
    genesNb = GenesDict["genesNb"]
    if "massAction" in ODEs:
        plt.subplot(len(ODEs),1,1)
        for solGenes in range(genesNb):
            plt.plot(GenesDict["massActionX"], GenesDict["massActionY"][solGenes], label=solGenes)
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("Genes concentrations", fontdict=font)
        plt.title("Mass Action law simulation")
        

    if "Hill" in ODEs:
        plt.subplot(len(ODEs),1,min(len(ODEs),2))
        for solGenes in range(genesNb):
            plt.plot(GenesDict["HillsX"], GenesDict["HillsY"][solGenes], label=solGenes)
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("Genes concentrations", fontdict=font)
        plt.title("Hill law Simulation")

    if "indirect" in ODEs:
        plt.subplot(len(ODEs), 1, len(ODEs))
        for solGenes in range(genesNb):
            plt.plot(GenesDict["indirectX"], GenesDict["indirectY"][solGenes], label=solGenes)
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("mRNA concentrations", fontdict=font)
        plt.title("indirect (massAction & Hill laws) Simulation")
    plt.legend()
        

    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    if saveName is not None:     
        plt.savefig(saveName)