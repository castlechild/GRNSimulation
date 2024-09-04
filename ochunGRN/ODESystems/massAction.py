#!/usr/bin/env python3
import numpy as np
from noise import stochastiqueNoise
# from noise import stochastiqueNoise


def massAction(t: float,
               G: np.ndarray,
               Adj: np.ndarray,
               k: np.ndarray,
               noise_amplitude: int | float = 2) -> np.ndarray:
    """
    Simulates the mass-action dynamics for gene regulatory networks.

        Args:
            - t (float): The current time point (not used in the calculation).
            - G (numpy.ndarray): Array of gene expression levels.
            - Adj (numpy.ndarray): Adjacency matrix indicating gene
            interactions.
            - k (numpy.ndarray): Rate constants for the interactions.
            - noise_amplitude (int or float, optional): Amplitude
            of stochastic noise. Default is 2.

        Returns:
            numpy.ndarray: An array containing the rate of change of gene
            expression levels.
    """
    del t
    genesNb = len(G)
    dG = np.zeros(genesNb)

    for i in range(genesNb):
        production, degradation = 0, 0
        for j in range(genesNb):
            stateEdgeprod = Adj[j][i]
            if stateEdgeprod != 0:
                production += G[j]*k[j]
            stateEdgedeg = Adj[i][j]
            if stateEdgedeg != 0:
                degradation += G[i]*k[i]
        dG[i] = production - degradation

        # Add stochastic noise
        noise = stochastiqueNoise(production, degradation, noise_amplitude)
        dG[i] += noise
    return dG

#####################################################################


def main():
    pass


if __name__ == "__main__":
    main()
