
import networkx as nx
import random as rd
import numpy as np
import pandas as pd
import pkg_resources
import ochunGRN as oG


excel_path = pkg_resources.resource_filename("ochunGRN", "GRN/41598_2021_3625_MOESM5_ESM.xlsx") 
document = pd.read_excel(excel_path)

arabidopsisThaliana=(document["Supplementary Table S1: Networks. A spreadsheet file with filtered networks"].tolist()[2:],document["Unnamed: 1"].tolist()[2:], "Arabidopsis thaliana" )
drosophilaMelanogaster=(document["Unnamed: 2"].tolist()[2:],document["Unnamed: 3"].tolist()[2:],"Drosophila Melanogaster")
escherichniaColi=(document["Unnamed: 4"].tolist()[2:],document["Unnamed: 5"].tolist()[2:], "Escherichnia coli")
homoSapiens=(document["Unnamed: 6"].tolist()[2:],document["Unnamed: 7"].tolist()[2:], "Homo sapiens")
saccharomycesCerevisiae=(document["Unnamed: 8"].tolist()[2:],document["Unnamed: 9"].tolist()[2:], "Saccharomyces cerevisiae")
for species in [arabidopsisThaliana, drosophilaMelanogaster, escherichniaColi, homoSapiens, saccharomycesCerevisiae]:
    oG.nanDestructor(species)
    GSpecies=oG.graphCreator(species)
    print(f"\nspecies : {species[2]}") 
    print("nodes :",GSpecies.number_of_nodes())
    print("edges :",GSpecies.number_of_edges())
    print("Auto Rg :",oG.findAutoRegulationGenes(species))
    print("Double Rg:", oG.findDoubleRegulationGenes(GSpecies))
    L=list(nx.clustering(GSpecies).items())
    print("moyenne coefficient de clustering", np.mean( [L[i][1] for i in range(GSpecies.number_of_nodes())]))

print("ratioFFL", oG.FFLratio(oG.graphCreator(escherichniaColi)))
G = oG.graphCreator(drosophilaMelanogaster)
L=list(nx.clustering(G).items())
print([L[i] for i in range(G.number_of_nodes())])

options = {             
    'node_color': 'black',
    'node_size': 0.4,
    'width': 0.02,
}
#G = graphCreator(saccharomycesCerevisiae)
#subax1 = plt.subplot(111)
#nx.draw(G, **options)
#plt.show()