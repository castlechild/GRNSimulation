import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from simulation import *

G = nx.DiGraph()
G.number_of_nodes()
#G.add_edge(1,2,color='r')
#G.add_edge(2,3,color='b')
#G.add_edge(3,4,color='b')
#G.add_edge(3,1,color='b')
#G.add_edge(4,4,color='g')
#G.add_edge(2,1,color='r')
#
#
#edges = G.edges()
#for u,v in edges:
#    G[u][v]['weight']=0.5
#    G[u][v]['colorweight']='b'
#for node in G :
#    print(node)
#print(G.degree(1))
#
#colors = [G[u][v]['color'] for u,v in edges]
#arrow = [G[u][v]['colorweight'] for u,v in edges]
#print(arrow)
##plt.subplot(1,1,1)
##nx.draw_circular(G, with_labels=True, font_weight='bold',edge_color=colors)
##plt.show()
#
#def twice(A):
#    A*=2
#
#M = np.ones((2,2))
#twice(M)
#print(M)
#
#M=nx.to_numpy_array(G)
#print(M)
#
#print([0]*2)

#def test(a,b=2,c=3):
#    print(a,b,c)
#
#test(1)
#test(1,4)
#test(1,2,5)
#test(c=1,b=2,a=0)
#a= None
#if a is not None:
#    print("hello")
#main()
def betterPrint(D:dict):
    for obj in D:
        print(obj,":")
        if type(D[obj]) == dict:
            for obj2 in D[obj]:
                print("     ",obj2,":",D[obj][obj2])
            print("")
        else :
            print(D[obj],"\n")

def sim():
    NB_GENES = 7
    AUTO_RG = 0.1
    DUO_RG = 0.2
    test = simulation(["massAction","Hill"], (0,0.50), 7, 0.1, 0.2)
    G = test["Graph"]
    Coeff = test["Coefficients"]
    betterPrint(test)
    simulation(["massAction","Hill"],(0,0.50),Graph=G,plot=True,saveName="1.png",Coeff=Coeff)
    simulation(["massAction","Hill"],(0,50),Graph=G,plot=True,saveName="1bis.png",Coeff=Coeff)
    simulation(["massAction","Hill"],(0,0.50),Graph=G,plot=True,saveName="2.png")
    simulation(["massAction"],(0,0.50),5,AUTO_RG,DUO_RG,plot=True,saveName="3.png")
    simulation(["Hill"],[0,50],5,AUTO_RG,DUO_RG,plot=True,saveName="4.png")

def powerPointImage():
    nbrCourbe = 100
    simPP = simulation(["massAction"],(0,0.50),nbrCourbe,.1,.2)
    for solGenes in range(nbrCourbe):
        plt.plot(simPP["massActionX"], simPP["massActionY"][solGenes], color='k', linewidth=0.5)
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.savefig("graphiquePowerPoint.png")


def main():
    powerPointImage()
    pass

if __name__ == "__main__":
    main()
        

