#!/usr/bin/env python3
import numpy as np
from noise import stochastiqueNoise
# from noise import stochastiqueNoise


def MaMatrice(M, K):
    MaPlus = np.abs(M) * K
    MaMinus = np.dot(-np.ones(np.shape(M)), MaPlus)
    MaMinus = np.diag(np.diag(MaMinus))
    # print("minus:",MaMinus)
    # print("plus:",MaPlus)
    return MaPlus + MaMinus


def massAction(t, G, A):
    del t
    return np.dot(A, G)


def massAction2(t, G, Adj, k, noise_amplitude=2):
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
                degradation += G[i]*k[j]
        dG[i] = production - degradation
        noise = stochastiqueNoise(production, degradation, noise_amplitude)
        dG[i] += noise
    return dG

#####################################################################


def main():
    Ma = np.random.randint(2, size=(3, 3))
    K = np.random.random((3, 3))
    print(Ma)
    print(K)
    print(MaMatrice(Ma, K))


if __name__ == "__main__":
    main()
