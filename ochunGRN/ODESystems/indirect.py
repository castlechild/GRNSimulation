#!/usr/bin/env python3
import numpy as np

from Hill import *

def indirect(t, mRNA, P, Adj, k_mRNA, k_P, Ka_P, K_degP, K_degMRNA, n, noise_amplitude=0.5):
    del t
    genesNb = len(mRNA)
    dmNRA = np.zeros(genesNb)
    dP = np.zeros(genesNb)
    for i in range(genesNb):
        # Terme de production et dégradation pour les protéines
        production_protein = k_mRNA[i]*mRNA[i]
        degradation_protein = K_degP[i]*P[i]
        dP[i] =  production_protein - degradation_protein

        # Terme de production et dégradation pour les ARNm
        production_mRNA = k_P[i]
        for j in range(genesNb):
            stateEdge = Adj[j][i]
            if stateEdge == 1:
                dP[i] -= k_mRNA[j]*P[i]
                production_mRNA*= Hill_activation(P[j],Ka_P[j],n)
            elif stateEdge == -1:
                dP[i] -= k_mRNA[j]*P[i]
                production_mRNA*= Hill_inhibition(P[j],Ka_P[j],n)
        degradation_mRNA = K_degMRNA[i] * mRNA[i]

        dmNRA[i] = production_mRNA - degradation_mRNA

        # Ajout du terme de bruit stochastique
        noise_mRNA = noise_amplitude * (np.sqrt(np.abs(production_mRNA)) * np.random.normal() +
                                        np.sqrt(np.abs(degradation_mRNA)) * np.random.normal())
        noise_protein = noise_amplitude *(np.sqrt(np.abs(production_protein)) * np.random.normal() +
                                           np.sqrt(np.abs(degradation_protein)) * np.random.normal())
        
        dmNRA[i] += noise_mRNA
        dP[i] += noise_protein
    
    return np.concatenate((dmNRA,dP))

def main():
    pass

if __name__ == "__main__":
    main()
