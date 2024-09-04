import networkx as nx
import matplotlib.pyplot as plt


def plotGraph(GenesDict: dict,
              actInhPlotBool: bool = False,
              saveName: str = None) -> None:
    """
    Plot the Gene Regulatory Network (GRN) graph from the provided gene
    dictionary.

    This function plots a graph using a circular layout.
    If the `actInhPlotBool` option is enabled, it also creates a separate plot
    to represent
    the activators and inhibitors in the network.

        Parameters:
            - GenesDict (dict): A dictionary containing information about
            the graph, under the key "Graph."
            - actInhPlotBool (bool, optional): If True, an additional graph
            showing activators and inhibitors is plotted. Default is False.
            - saveName (str, optional): If provided, the graph will be saved
            with this file name.

        Returns:
            - None: The graph is displayed or saved depending on the
            parameters.
    """
    # Extract the graph from the gene dictionary
    Graph = GenesDict["Graph"]

    # Create a new figure for plotting
    plt.figure()

    # Get edges and their colors from the graph
    edges = Graph.edges()
    colors = [Graph[u][v]['color'] for u, v in edges]

    # Plot the main GRN graph using a circular layout
    plt.subplot(1 + actInhPlotBool, 1, 1)
    nx.draw_circular(Graph, with_labels=True,
                     font_weight='bold', edge_color=colors)
    plt.title("GRN Graph")

    # Plot activator/inhibitor graph if specified
    if actInhPlotBool:
        # Get activator/inhibitor edge colors
        acInColors = [Graph[u][v]['acInColor'] for u, v in edges]

        # Create a new subplot for activator/inhibitor graph
        plt.subplot(2, 1, 2)
        nx.draw_circular(Graph, with_labels=True,
                         font_weight='bold', edge_color=acInColors,
                         connectionstyle="arc3,rad=0.05")
        plt.title("GRN Graph Activator/Inhibitor")

    # Save the graph to a file if saveName is provided
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
            - ODEs (list, optional): A list of ODE types to plot.
            If not provided, it defaults to the ODEs in `GenesDict["ODEs"]`.
            - saveName (str, optional): If provided,
            the graph will be saved with this file name.

        Returns:
            - None: The graph is displayed or saved depending on
            the parameters.
    """
    # Set the font for plot labels
    font = {'family': 'serif', 'color': 'darkred', 'size': 8}

    # If ODEs are not provided, use all ODEs existing in GenesDict
    if ODEs is None:
        ODEs = GenesDict["ODEs"]

    # Ensure ODEs is a list for consistency
    if isinstance(ODEs, str):
        ODEs = [ODEs]

    # Get the number of genes from the gene dictionnary
    genesNb = GenesDict["genesNb"]

    # Loop through each ODE type and plot results
    for i in range(len(ODEs)):
        # Create a subplot for each ODE
        plt.subplot(len(ODEs), 1, i+1)

        # Plot the mRNA concentrations for each gene
        for solGenes in range(genesNb):
            plt.plot(GenesDict[f"{ODEs[i]}X"],
                     GenesDict[f"{ODEs[i]}Y"][solGenes], label=solGenes)

        # Set plot labels and titles
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("mRNA concentrations", fontdict=font)
        plt.title(f"{ODEs[i]} law simulation")

    # Create a legend with gene labels
    labels = []
    for solGenes in range(genesNb):
        labels.append(f"Gene {solGenes}")
    handles, _ = plt.gca().get_legend_handles_labels()
    plt.figlegend(handles, labels, loc='upper right')

    # Save the plot to a file if saveName is provided
    if saveName is not None:
        plt.savefig(saveName)


def plotProt(GenesDict, saveName=None):
    """
    Plot the protein concentrations over time from
    the "indirect" ODE simulation.

    This function plots the simulated protein concentrations using the results
    from the "indirect" ODE model in `GenesDict`.

        Parameters:
            - GenesDict (dict): A dictionary containing the simulation results,
            including the "indirect" ODE results.
            - saveName (str, optional): If provided, the graph will be saved
            with this file name.

        Raises:
            - KeyError: If "indirect" is not in `GenesDict["ODEs"]`.

        Returns:
            - None: The graph is displayed or saved depending
            on the parameters.
    """
    # Ensure that "indirect" ODE results are available
    if "indirect" not in GenesDict["ODEs"]:
        raise

    # Get the number of genes from the gene dictionnary
    genesNb = GenesDict["genesNb"]

    # Set the font for plot labels
    font = {'family': 'serif', 'color': 'darkred', 'size': 8}

    # Plot protein concentrations for each gene
    for solGenes in range(genesNb):
        plt.plot(GenesDict["indirectX"],
                 GenesDict["indirectProt"][solGenes], label=solGenes)

        # Set plot labels and title
        plt.xlabel("time (h)", fontdict=font)
        plt.ylabel("protein concentrations", fontdict=font)
        plt.title("protein concentrations from indirect law simulation")

    # Create a legend with gene labels
    labels = []
    for solGenes in range(genesNb):
        labels.append(f"Gene {solGenes}")
    handles, _ = plt.gca().get_legend_handles_labels()
    plt.figlegend(handles, labels, loc='upper right')

    # Save the plot to a file if saveName is provided
    if saveName is not None:
        plt.savefig(saveName)


def main():
    pass


if __name__ == "__main__":
    main()
