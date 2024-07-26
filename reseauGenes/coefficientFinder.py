#!/usr/bin/env python3

import networkx as nx
import random as rd
import numpy as np
import pandas as pd

document = pd.read_excel("reseauGenes/41586_2011_BFnature10098_MOESM304_ESM.xls")
Attribut = document.columns.tolist()
document = document.T
i = rd.randint(0, 5027)
gene = document[i].tolist()
    
for j in range(len(gene)):
    print(Attribut[j],":\n",gene[j],"\n")