#!/usr/bin/env python3
import numpy as np
from noise import stochastiqueNoise


def HillEquation(t: float,
                 G: np.ndarray,
                 Adj: np.ndarray,
                 K: np.ndarray,
                 Ka: np.ndarray,
                 K_synt: np.ndarray,
                 K_deg: np.ndarray,
                 n: int,
                 noise_amplitude: int | float = 2) -> np.ndarray:
    """
    Computes the rate of change of gene expression levels based on Hill
    functions for activation and inhibition.
    The function calculates the derivative of gene expression levels using
    a system of Hill equations.
    This model incorporates both activation and inhibition effects based
    on the Hill function.

        Parameters:
            - t (float): The current time point (not used in the calculation).
            - G (numpy.ndarray): A 1D array of gene expression levels
            at time t.
            - Adj (numpy.ndarray): A 2D adjacency matrix representing
            interactions between genes. `Adj[i][j]` indicates the effect of
            gene j on gene i (1 for activation, -1 for inhibition).
            - K (numpy.ndarray): A 1D array of baseline rate constants
            for each gene.
            - Ka (numpy.ndarray): A 1D array of Hill constants (Ka)
            for the activation and inhibition functions.
            - K_synt (numpy.ndarray): A 1D array of synthesis rate constants
            for each gene.
            - K_deg (numpy.ndarray): A 1D array of degradation rate constants
            for each gene.
            - n (int): The Hill coefficient, which defines the steepness
            of the activation/inhibition response.

        Returns:
            numpy.ndarray: A 1D array where each element represents the rate
            of change of the corresponding gene's expression level.
    """
    del t
    dG = np.zeros(len(G))
    for i in range(len(dG)):
        production = K[i]
        for j in range(len(dG)):
            stateEdge = Adj[j][i]
            if stateEdge == 1:
                production *= Hill_activation(G[j], Ka[j], n)
            elif stateEdge == -1:
                production *= Hill_inhibition(G[j], Ka[j], n)
        production += K_synt[i]
        degradation = K_deg[i]*G[i]
        dG[i] = production - degradation
        noise_mRNA = stochastiqueNoise(production, degradation,
                                       noise_amplitude)
        dG[i] += noise_mRNA
    return dG


def Hill_inhibition(G_i: float, ka: float, n: int) -> float:
    """
    Computes the Hill inhibition function value for a given gene
    expression level.

        Args:
            - G_i (float): Expression level of the inhibiting gene.
            - ka (float): Hill constant for inhibition.
            - n (int): Hill coefficient.

        Returns:
            float: The value of the Hill inhibition function.
    """
    return ka**n / (ka**n + G_i**n)


def Hill_activation(G_a: float, ka: float, n: int) -> float:
    """
    Computes the Hill activation function value for a given gene
    expression level.

        Args:
            - G_a (float): Expression level of the activating gene.
            - ka (float): Hill constant for activation.
            - n (int): Hill coefficient.

        Returns:
            float: The value of the Hill activation function.
    """
    return 1 + (G_a**n / (ka**n + G_a**n))


def main():
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt
    import networkx as nx
    G = nx.DiGraph()
    G.add_edge(0, 1)
    plt.subplot(1, 2, 1)
    nx.draw(G, with_labels=True)
    A = np.array([[0, 0], [1, 0]])

    def equationAct(t, G): HillEquation(t, G, A, [1, 1], [1, 1], [0, 0],
                                        [1, 1], 1)

    def equationInh(t, G): HillEquation(t, G, -A, [1, 1], [1, 1], [0, 0],
                                        [1, 1], 1)
    solutionAct = solve_ivp(equationAct, [0, 50], [1, 1], max_step=0.5)
    solutionInh = solve_ivp(equationInh, [0, 50], [1, 1], max_step=0.5)
    plt.subplot(1, 2, 2)
    plt.plot(solutionAct.t, solutionAct.y[0], label="0 activator", color='b')
    plt.plot(solutionAct.t, solutionAct.y[1], label="1 activator", color='r')
    plt.plot(solutionInh.t, solutionInh.y[0], label="0 inhibitor",
             linestyle='dashed', color='b')
    plt.plot(solutionInh.t, solutionInh.y[1], label="1 inhibitor",
             linestyle='dashed', color='r')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
