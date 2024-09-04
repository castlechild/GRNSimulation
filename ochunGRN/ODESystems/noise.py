#!/usr/bin/env python3
import numpy as np
import random as rd


def stochastiqueNoise(production: float,
                      degradation: float,
                      noise_amplitude: int | float) -> float:
    """
    Computes stochastic noise for gene expression dynamics.

        Args:
            - production (float): The production term for a gene.
            - degradation (float): The degradation term for a gene.
            - noise_amplitude (int or float): Amplitude of stochastic noise.

        Returns:
            float: The computed stochastic noise.
    """
    return noise_amplitude * (
        np.sqrt(np.abs(production)) * np.random.normal()
        + np.sqrt(np.abs(degradation)) * np.random.normal()
        )


def physicalNoise(concentrationsList: np.ndarray,
                  noise_amplitude: int | float) -> None:
    """
    Adds physical noise to a list of gene concentrations.

        Args:
            - concentrationsList (numpy.ndarray): Array of gene concentrations.
            - noise_amplitude (int or float): Amplitude of noise to add.

        Returns:
            None: Modifies the input array in place.
    """
    geneNB = len(concentrationsList)
    for gene in range(geneNB):
        valuesNB = len(concentrationsList[gene])
        for value in range(valuesNB):
            concentrationsList[gene][value] += (
                noise_amplitude * np.random.normal()
                )
            if concentrationsList[gene][value] < 0:
                concentrationsList[gene][value] = 0


def dropOut(concentrationsList: np.ndarray, droupOutRate: float) -> None:
    """
    Applies dropout noise to a list of gene concentrations.

        Args:
            - concentrationsList (numpy.ndarray): Array of gene concentrations.
            - droupOutRate (float): Probability of setting a concentration
            to zero.

        Returns:
            None: Modifies the input array in place.
    """
    geneNB = len(concentrationsList)
    for gene in range(geneNB):
        valuesNB = len(concentrationsList[gene])
        for value in range(valuesNB):
            if rd.random() < droupOutRate:
                concentrationsList[gene][value] = 0


def dataSelector(concentrationsList: np.ndarray,
                 timeList: np.ndarray,
                 dataNB: int) -> np.ndarray:
    """
    Selects a subset of data points from gene concentration and time lists.

        Args:
            - concentrationsList (numpy.ndarray): Array of gene concentrations.
            - timeList (numpy.ndarray): Array of time points corresponding
            to concentrations.
            - dataNB (int): Number of data points to select.

        Returns:
            numpy.ndarray: A tuple of two arrays - selected concentrations
            and times.
    """
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
