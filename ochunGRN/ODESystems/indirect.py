#!/usr/bin/env python3
import numpy as np

from Hill import *

def indirect(t, mRNA, P, Adj, k_mRNA, k_P, Ka_P, K_degP, K_degMRNA, n):
    del t
    genesNb = len(mRNA)
    dmNRA = np.zeros(genesNb)
    dP = np.zeros(genesNb)
    for i in range(genesNb):
        dP[i] = k_mRNA[i]*mRNA[i] - K_degP[i]*P[i]
        dmNRA[i]= k_P[i]
        for j in range(genesNb):
            stateEdge = Adj[i][j]
            if stateEdge == 1:
                dP[i] -= k_mRNA[j]*P[i]
                dmNRA[i]*= Hill_activation(P[j],Ka_P[j],n)
            elif stateEdge == -1:
                dP[i] -= k_mRNA[j]*P[i]
                dmNRA[i]*= Hill_inhibition(P[j],Ka_P[j],n)
        dmNRA[i] -= K_degMRNA[i]*mRNA[i]
    return np.concatenate((dmNRA,dP))

def main():
    pass

if __name__ == "__main__":
    main()
