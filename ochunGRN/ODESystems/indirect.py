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
    del t
    genesNb = len(mRNA)
    dmNRA = np.zeros(genesNb)
    dP = np.zeros(genesNb)
    for i in range(genesNb):
        # Terme de production et dégradation pour les protéines
        production_protein = kTranslation[i]*mRNA[i]
        degradation_protein = degP[i]*P[i]

        # Terme de production et dégradation pour les ARNm
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

        # Ajout du terme de bruit stochastique
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
