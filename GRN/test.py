from GRNCreation import *

G1 = nx.barabasi_albert_graph(10,2)
G2 = BarabasiAlbertAlgorithm(10,2)
#plt.subplot(1,2,1)
#nx.draw(G1)
#plt.subplot(1,2,2)
#nx.draw(G2)

nx.draw(G1,with_labels=True)
partition = nx.community.greedy_modularity_communities(G1)

print(partition)
plt.show()