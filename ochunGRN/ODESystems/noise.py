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


def dropOut(concentrationsList, droupOutRate):
    geneNB = len(concentrationsList)
    for gene in range(geneNB):
        valuesNB = len(concentrationsList[gene])
        for value in range(valuesNB):
            if rd.random() < droupOutRate:
                concentrationsList[gene][value] = 0


def main():
    pass


if __name__ == "__main__":
    main()
