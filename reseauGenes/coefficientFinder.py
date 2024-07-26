#!/usr/bin/env python3

import networkx as nx
import random as rd
import numpy as np
import pandas as pd

document = pd.read_excel("reseauGenes/41586_2011_BFnature10098_MOESM304_ESM.xls")
Attribut = document.columns.tolist()
document = document.T

def getCoefficient(GenesNb):
    resDict = {}
    resDict["ProtsHalfTime"] = []
    resDict["mRNAsHalfTime"] = []
    resDict["TranscriptionsRate"] = []
    resDict["TranslationsRate"] = []
    resDict["mRNAAvg"] = []
    for i in range(GenesNb):
        randomNb = rd.randint(0,5027)
        L=document[randomNb].tolist()
        while True in pd.isna([L[i] for i in [16,19,22,25,28]]):
            randomNb = rd.randint(0,5027)
            L=document[randomNb].tolist()
        resDict["ProtsHalfTime"].append(L[19])
        resDict["mRNAsHalfTime"].append(L[22])
        resDict["TranscriptionsRate"].append(L[25])
        resDict["TranslationsRate"].append(L[28])
        resDict["mRNAAvg"].append(L[16])
    return resDict

def main():
    #i = rd.randint(0, 5027)
    #gene = document[i].tolist()    
    #for j in range(len(gene)):
    #    print(Attribut[j],":\n",gene[j],"\n")
    print(getCoefficient(10))
if __name__ == "__main__" :
    main()