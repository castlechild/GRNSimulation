#!/usr/bin/env python3
import numpy as np

# dG/dt = Hill_inhibition(G_i) + Hill_activation(G_a) + k_synthesis - k_degratation . G
# Hill_activation(G) = k . G^n / (K^n + G^n)
# Hill_inhibition(G) = k . K^n / (K^n + G^n)

def HillEquation(t, G, Adj, K, Ka, K_synt, K_deg, n):
    dG = np.zeros(len(G))
    for i in range(len(dG)):
        dG[i]+= K_synt[i]
        dG[i]-= K_deg[i]*G[i]
        for j in range(len(dG)):
            stateEdge = Adj[i][j]
            if stateEdge == 1:
                dG[i]+= Hill_activation(G[j],K[j],Ka[j],n)
            elif stateEdge == -1:
                dG[i]+= Hill_inhibition(G[j],K[j],Ka[j],n)
    return dG

def Hill_inhibition(G_i, k, ka, n):
    return k * ka**n / (ka**n + G_i**n)

def Hill_activation(G_a, k, ka, n):
    return k * G_a**n / (ka**n + G_a**n)


def main():
    pass

if __name__ == "__main__":
    main()