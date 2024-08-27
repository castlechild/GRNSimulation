#!/usr/bin/env python3
import numpy as np
import random as rd


def stochastiqueNoise(production, degradation, noise_amplitude):
    return noise_amplitude * (
        np.sqrt(np.abs(production)) * np.random.normal()
        + np.sqrt(np.abs(degradation)) * np.random.normal()
        )


def physicalNoise(concentrationsList, noise_amplitude):
    geneNB = len(concentrationsList)
    for gene in range(geneNB):
        valuesNB = len(concentrationsList[gene])
        for value in range(valuesNB):
            concentrationsList[gene][value] += (
                noise_amplitude * np.random.normal()
                )
            if concentrationsList[gene][value] < 0:
                concentrationsList[gene][value] = 0


def dropOut(concentrationsList, droupOutRate):
    geneNB = len(concentrationsList)
    for gene in range(geneNB):
        valuesNB = len(concentrationsList[gene])
        for value in range(valuesNB):
            if rd.random() < droupOutRate:
                concentrationsList[gene][value] = 0


def dataSelector(concentrationsList, timeList, dataNB):
    geneNB = len(concentrationsList)
    res = ([[] for _ in range(geneNB)], [])
    valuesNB = len(timeList)
    tf = timeList[-1]
    if valuesNB > dataNB:
        a = -1
        for value in range(valuesNB):
            ti = timeList[value]
            b = (ti+1) * dataNB / (tf+1)
            if int(b) != a:
                a += 1
                for gene in range(geneNB):
                    res[0][gene].append(concentrationsList[gene][value])
                res[1].append(timeList[value])
    else:
        res = (concentrationsList, timeList)
    return res


def main():
    pass


if __name__ == "__main__":
    main()
