#!/usr/bin/env python3
import numpy as np
from noise import stochastiqueNoise

from Hill import Hill_activation, Hill_inhibition


def indirect(t: float,
             mRNA: np.ndarray,
             P: np.ndarray,
             Adj: np.ndarray,
             kTranscription: np.ndarray,
             kTranslation: np.ndarray,
             degP: np.ndarray,
             degMRNA: np.ndarray,
             Pavg: np.ndarray,
             mRNAavg: np.ndarray,
             n: int,
             noise_amplitude: int | float = 2) -> np.ndarray:
    """
    Simulates indirect gene regulation dynamics involving mRNA
    and protein levels.

        Args:
            - t (float): The current time point (not used in the calculation).
            - mRNA (numpy.ndarray): Array of mRNA levels for each gene.
            - P (numpy.ndarray): Array of protein levels for each gene.
            - Adj (numpy.ndarray): Adjacency matrix indicating
            gene interactions.
            - kTranscription (numpy.ndarray): Transcription rate constants.
            - kTranslation (numpy.ndarray): Translation rate constants.
            - degP (numpy.ndarray): Protein degradation rate constants.
            - degMRNA (numpy.ndarray): mRNA degradation rate constants.
            - Pavg (numpy.ndarray): Average protein levels.
            - mRNAavg (numpy.ndarray): Average mRNA levels.
            - n (int): Hill coefficient for activation/inhibition.
            - noise_amplitude (int or float, optional): Amplitude
            of stochastic noise. Default is 2.

        Returns:
            numpy.ndarray: An array containing the rate of change of mRNA
            and protein levels.
    """
    del t
    genesNb = len(mRNA)
    dmNRA = np.zeros(genesNb)
    dP = np.zeros(genesNb)
    for i in range(genesNb):

        # Protein production and degradation terms
        degradation_protein = degP[i]*P[i]
        production_protein = kTranslation[i]*mRNA[i]
        # mRNA production and degradation terms
        production_mRNA = kTranscription[i]
        for j in range(genesNb):
            stateEdge = Adj[j][i]
            if stateEdge == 1:
                #  degradation_protein += k_P[i]
                production_mRNA *= Hill_activation(P[j], mRNAavg[j], n)
            elif stateEdge == -1:
                #  degradation_protein += k_P[i]
                production_mRNA *= Hill_inhibition(P[j], mRNAavg[j], n)
        degradation_mRNA = degMRNA[i] * mRNA[i]

        dP[i] = production_protein - degradation_protein
        dmNRA[i] = production_mRNA - degradation_mRNA

        # Add stochastic noise
        noise_mRNA = stochastiqueNoise(production_mRNA,
                                       degradation_mRNA, noise_amplitude)
        noise_protein = stochastiqueNoise(production_protein,
                                          degradation_protein, noise_amplitude)

        dmNRA[i] += noise_mRNA
        dP[i] += noise_protein

    return np.concatenate((dmNRA, dP))


def main():
    pass


if __name__ == "__main__":
    main()
