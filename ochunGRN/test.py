import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
#from simulation import *
import ochunGRN as oGRN

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
    test = Melvin.simulation(["massAction","Hill"], (0,0.50), 7, 0.1, 0.2)
    G = test["Graph"]
    Coeff = test["Coefficients"]
    betterPrint(test)
    Melvin.simulation(["massAction","Hill"],(0,0.50),Graph=G,plot=True,saveName="1.png",Coeff=Coeff)
    Melvin.simulation(["massAction","Hill"],(0,50),Graph=G,plot=True,saveName="1bis.png",Coeff=Coeff)
    Melvin.simulation(["massAction","Hill"],(0,0.50),Graph=G,plot=True,saveName="2.png")
    Melvin.simulation(["massAction"],(0,0.50),5,AUTO_RG,DUO_RG,plot=True,saveName="3.png")
    Melvin.simulation(["Hill"],[0,50],5,AUTO_RG,DUO_RG,plot=True,saveName="4.png")

def graphsim():
    coeff = {'ProtsDeg': [0.026364906084509605,
   0.045341890425556114,
   0.05210808999163947,
   0.0876759803083134,
   0.025527320620452983,
   0.014292203763245024,
   0.034509847788508914],
  'mRNAsDeg': [0.32812379637426803,
   0.04437444500700016,
   0.21993990083562912,
   0.17655600029626656,
   0.2767541217951262,
   0.1033026707357995,
   0.14029896452114038],
  'TranscriptionsRate': [2.87, 3.05, 0.61, 4.55, 4.1, 0.68, 1.8],
  'TranslationsRate': [311.21, 6.03, 5.25, 50.08, 213.21, 72.25, 25.22],
  'mRNAAvg': [15.91, 69.72, 5.46, 42.67, 26.46, 9.54, 20.11],
  'ProtAvg': [137923.26,
   9217.4,
   218.75,
   32435.06,
   158546.41,
   22593.89,
   12535.34]}
    Adj = np.array([[ 0.,  0.,  0.,  1.,  0.,  0., -1.],
        [ 0.,  0.,  1.,  0.,  0.,  1.,  0.],
        [-1.,  1.,  0.,  0.,  0.,  0.,  0.],
        [-1.,  0.,  0.,  0.,  1.,  0.,  0.],
        [ 0.,  0., -1.,  0.,  0.,  1.,  0.],
        [ 0.,  1.,  0.,  0.,  0.,  0., -1.],
        [-1.,  0.,  0.,  0.,  0.,  0.,  0.]])
    #MelDict=Melvin.simulation(["massAction","Hill","indirect"],(0,100),Coeff=coeff,M=Adj)
    MelDict=Melvin.simulation(["indirect"],(0,100),7,0.1,0.2,Coeff=coeff,plot=True)


def powerPointImage():
    nbrCourbe = 100
    simPP = Melvin.simulation(["massAction"],(0,0.50),nbrCourbe,.1,.2)
    for solGenes in range(nbrCourbe):
        plt.plot(simPP["massActionX"], simPP["massActionY"][solGenes], color='k', linewidth=0.5)
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.savefig("graphiquePowerPoint.png")


def main():
    graphsim()
    pass

if __name__ == "__main__":
    main()
        

