#!/usr/bin/env python3
import numpy as np
from scipy.integrate import solve_ivp
from coefficientFinder import getCoefficient
from massAction import MaMatrice, massAction
from Hill import HillEquation
from indirect import indirect
from noise import physicalNoise, dropOut, dataSelector


def simulationODEs(GenesDict: dict,
                   ODEs: any,
                   T: tuple,
                   Coeff: dict = None,
                   normalisationBool: bool = False):
    """
    Simulate the temporal evolution of a gene regulatory network using
    ordinary differential equations (ODEs).

    This function simulates the dynamics of a gene network based on
    the specified types of ODEs. It supports mass action models,
    the Hill equation, and an indirect model (a hybrid of the first two),
    each applied to the provided genetic data.

    Parameters:
    - GenesDict (dict): A dictionary containing information about the genes, including the adjacency matrix and the number of genes.
    - ODEs (list, tuple, str): A list of strings indicating the types of ODEs to simulate (e.g., "massAction", "Hill", "indirect").
    - T (tuple): A tuple of two values indicating the start and end of the simulation (t0, tf).
    - Coeff (dict, optional): A dictionary containing the coefficients for the simulation. If no coefficients are provided, they will be generated automatically.
    # noqa: E501

    Raises:
    - ValueError: If the specified ODEs are invalid or if the time interval
    is incorrect.

    Returns:
    - None: The simulation results are stored in `GenesDict`.
"""

    def otherODE(L: list):
        """
        Check if the list of ODEs contains elements other than
        the supported models.

        Parameters:
        - L (list): The list of ODEs to check.

        Returns:
        - bool: True if unsupported ODEs are detected, otherwise False.
        """
        for ode in L:
            if ode not in ["massAction", "Hill", "indirect"]:
                return True
        return False
    if isinstance(ODEs, str):
        ODEs = [ODEs]
    if len(ODEs) == 0 and otherODE(ODEs):
        raise ValueError("ODEs not valid")
    if len(T) != 2 and T[0] >= T[1]:
        raise ValueError("the temporal pair is invalid")
    genesNb = GenesDict["genesNb"]
    if Coeff is None:
        Coeff = getCoefficient(genesNb)
    G0 = Coeff["mRNAAvg"]

    M = GenesDict["AdjMatrice"]
    t0, tf = T

    GenesDict["Coefficients"] = Coeff
    GenesDict["ODEs"] = ODEs
    if "massAction" in ODEs:
        RatioCoeff = [Coeff["TranscriptionsRate"][i]/Coeff["mRNAAvg"][i]
                      for i in range(genesNb)]
        K = np.resize(RatioCoeff, (genesNb, genesNb))
        Ma = MaMatrice(np.transpose(M), K)

        def equation(t, G): return massAction(t, G, Ma)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        dataY, dataX = dataSelector(solution.y, solution.t, 200)
        physicalNoise(dataY, 0.5)
        if normalisationBool:
            normalisation(dataY)
        dropOut(dataY, 0.02)
        GenesDict["massActionY"] = dataY
        GenesDict["massActionX"] = dataX

    if "Hill" in ODEs:
        K = Coeff["TranscriptionsRate"]
        Kdeg = Coeff["mRNAsDeg"]

        def equation(t, G):
            return HillEquation(t, G, M, K, G0, [0]*genesNb, Kdeg, 2)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        dataY, dataX = dataSelector(solution.y, solution.t, 200)
        physicalNoise(dataY, 0.5)
        if normalisationBool:
            normalisation(dataY)
        dropOut(dataY, 0.02)
        GenesDict["HillY"] = dataY
        GenesDict["HillX"] = dataX

    if "indirect" in ODEs:
        k_P = Coeff["TranslationsRate"]
        Ka_P = G0
        k_mRNA = Coeff["TranscriptionsRate"]
        K_degP = Coeff["ProtsDeg"]
        K_degMRNA = Coeff["mRNAsDeg"]

        def equation(t, G):
            mRNA = G[:genesNb]
            P = G[genesNb:]
            return indirect(t, mRNA, P, M, k_mRNA, k_P,
                            Ka_P, K_degP, K_degMRNA, 2)
        G0_indirect = np.concatenate((G0, Coeff["ProtAvg"]))
        solution = solve_ivp(equation, [t0, tf], G0_indirect, max_step=0.5)
        dataY, dataX = dataSelector(solution.y[:genesNb], solution.t, 200)
        physicalNoise(dataY, 20)
        if normalisationBool:
            normalisation(dataY)
        dropOut(dataY, 0.02)
        GenesDict["indirectY"] = dataY
        GenesDict["indirectX"] = dataX


def normalisation(YLists):
    """
    Normalize the lists of values by dividing each element by the maximum
    of its respective list.

    Parameters:
    - YLists (list of list): A list containing lists of values to normalize.

    Returns:
    - None: The normalization is performed in place on `YLists`.
    """
    listNb = len(YLists)
    listSize = len(YLists[0])
    for i in range(listNb):
        ratio = max(YLists[i])
        if ratio != 0:
            for j in range(listSize):
                YLists[i][j] /= ratio
