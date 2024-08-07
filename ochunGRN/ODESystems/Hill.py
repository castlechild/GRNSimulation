#!/usr/bin/env python3
import numpy as np

# dG/dt = k * \prod( Hill_inhibition(G_i) * Hill_activation(G_a) )  - k_degratation . G
# Hill_activation(G) = 1 + G^n / (K^n + G^n)
# Hill_inhibition(G) =  K^n / (K^n + G^n)

def HillEquation(t, G, Adj, K, Ka, K_synt, K_deg, n):
    del t
    dG = np.zeros(len(G))
    for i in range(len(dG)):
        dG[i]=K[i]
        for j in range(len(dG)):
            stateEdge = Adj[i][j]
            if stateEdge == 1:
                dG[i]*= Hill_activation(G[j],Ka[j],n)
            elif stateEdge == -1:
                dG[i]*= Hill_inhibition(G[j],Ka[j],n)
        dG[i]+= K_synt[i]
        dG[i]-= K_deg[i]*G[i]
    return dG

def Hill_inhibition(G_i, ka, n):
    return ka**n / (ka**n + G_i**n)

def Hill_activation(G_a, ka, n):
    return  1 + (G_a**n / (ka**n + G_a**n))


def main():
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt
    import networkx as nx
    G = nx.DiGraph()
    G.add_edge(0,1)
    plt.subplot(1,2,1)
    nx.draw(G,with_labels=True)
    A = np.array([[0,0],[1,0]])
    equationAct = lambda t,G : HillEquation(t, G, A, [1,1], [1,1], [0,0], [1,1], 1)
    equationInh = lambda t,G : HillEquation(t, G, -A, [1,1], [1,1], [0,0], [1,1], 1)
    solutionAct = solve_ivp(equationAct, [0, 50], [1,1], max_step=0.5)
    solutionInh = solve_ivp(equationInh, [0, 50], [1,1], max_step=0.5)
    plt.subplot(1,2,2)
    plt.plot(solutionAct.t, solutionAct.y[0],label="0 activator",color='b')
    plt.plot(solutionAct.t, solutionAct.y[1],label="1 activator",color='r')
    plt.plot(solutionInh.t, solutionInh.y[0],label="0 inhibitor",linestyle = 'dashed',color='b')
    plt.plot(solutionInh.t, solutionInh.y[1],label="1 inhibitor",linestyle = 'dashed',color='r')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()