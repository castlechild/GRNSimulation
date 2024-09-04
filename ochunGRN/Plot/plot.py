import networkx as nx
import matplotlib.pyplot as plt


def plotGraph(GenesDict: dict,
              actInhPlotBool: bool = False,
              saveName: str = None) -> None:
    """
    Plot the Gene Regulatory Network (GRN) graph from the provided
    gene dictionary.

    This function plots a graph using a circular layout. If the actInhPlotBool
    option is enabled, it also creates a separate plot to represent
    the activators and inhibitors in the network.

    Parameters:
    - GenesDict (dict): A dictionary containing information about the graph,
    under the key "Graph."
    - actInhPlotBool (bool, optional): If True, an additional graph showing
    activators and inhibitors is plotted.
    - saveName (str, optional): If provided, the graph will be saved with this
    file name.

    Returns:
    - None: The graph is displayed or saved depending on the parameters
    """
    Graph = GenesDict["Graph"]
    plt.figure()
    edges = Graph.edges()
    colors = [Graph[u][v]['color'] for u, v in edges]

    plt.subplot(1 + actInhPlotBool, 1, 1)
    nx.draw_circular(Graph, with_labels=True,
                     font_weight='bold', edge_color=colors)
    plt.title("GRN Graph")

    if actInhPlotBool:
        acInColors = [Graph[u][v]['acInColor'] for u, v in edges]
        plt.subplot(2, 1, 2)
        nx.draw_circular(Graph, with_labels=True,
                         font_weight='bold', edge_color=acInColors,
                         connectionstyle="arc3,rad=0.05")
        plt.title("GRN Graph Activator/Inhibitor")
    if saveName is not None:
        plt.savefig(saveName)


def plotSim(GenesDict: dict,
            ODEs: list = None,
            saveName: str = None) -> None:
    """
    Plot the simulation results of the ODEs for the gene regulatory network.

    This function plots the simulated mRNA concentrations for each type of ODE
    specified in GenesDict. The results are displayed on separate subplots
    for each ODE.

    Parameters:
    - GenesDict (dict): A dictionary containing the simulation results
    for each type of ODE.
    - saveName (str, optional): If provided,
    the graph will be saved with this file name.

    Returns:
    - None: The graph is displayed or saved depending on the parameters.
    """
    font = {'family': 'serif', 'color': 'darkred', 'size': 8}
    if ODEs is None:
        ODEs = GenesDict["ODEs"]
    if isinstance(ODEs, str):
        ODEs = [ODEs]
    genesNb = GenesDict["genesNb"]
    for i in range(len(ODEs)):
        plt.subplot(len(ODEs), 1, i+1)
        for solGenes in range(genesNb):
            plt.plot(GenesDict[f"{ODEs[i]}X"],
                     GenesDict[f"{ODEs[i]}Y"][solGenes], label=solGenes)
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


def plotProt(GenesDict, saveName=None):
    if "indirect" not in GenesDict["ODEs"]:
        raise
    genesNb = GenesDict["genesNb"]
    font = {'family': 'serif', 'color': 'darkred', 'size': 8}
    for solGenes in range(genesNb):
        plt.plot(GenesDict["indirectX"],
                 GenesDict["indirectProt"][solGenes], label=solGenes)
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("protein concentrations", fontdict=font)
        plt.title("protein concentrations from indirect law simulation")
    labels = []
    for solGenes in range(genesNb):
        labels.append(f"Gene {solGenes}")
    handles, _ = plt.gca().get_legend_handles_labels()
    plt.figlegend(handles, labels, loc='upper right')
    if saveName is not None:
        plt.savefig(saveName)


def main():
    pass


if __name__ == "__main__":
    main()
