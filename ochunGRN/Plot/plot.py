import networkx as nx
import matplotlib.pyplot as plt

def plotGraph(GenesDict, actInhPlotBool=False, saveName=None):
    Graph = GenesDict["Graph"]
    plt.figure()
    edges = Graph.edges()
    colors = [Graph[u][v]['color'] for u,v in edges]
    
    plt.subplot(1+ actInhPlotBool, 1, 1)
    nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=colors)
    plt.title("GRN Graph")
    
    if actInhPlotBool:
        acInColors = [Graph[u][v]['acInColor'] for u,v in edges]
        plt.subplot(2,1,2)
        nx.draw_circular(Graph, with_labels=True, font_weight='bold',edge_color=acInColors, connectionstyle="arc3,rad=0.05" )
        plt.title("GRN Graph Activator/Inhibitor")
    if saveName is not None:     
        plt.savefig(saveName)

def plotSim(GenesDict, saveName=None):
    font = {'family':'serif','color':'darkred','size':8}
    ODEs = GenesDict["ODEs"]
    genesNb = GenesDict["genesNb"]
    for i in range(len(ODEs)):
        plt.subplot(len(ODEs),1,i+1)
        for solGenes in range(genesNb):
            plt.plot(GenesDict[f"{ODEs[i]}X"], GenesDict[f"{ODEs[i]}Y"][solGenes], label=solGenes)
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("mRNA concentrations", fontdict=font)
        plt.title(f"{ODEs[i]} law simulation")
    labels = []
    for solGenes in range(genesNb):
        labels.append(f"Gene {solGenes}")
    handles, _ = plt.gca().get_legend_handles_labels()
    plt.figlegend(handles, labels, loc='upper right')
    if saveName is not None:     
        plt.savefig(saveName)