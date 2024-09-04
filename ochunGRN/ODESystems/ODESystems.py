#!/usr/bin/env python3
import numpy as np
from scipy.integrate import solve_ivp
from extensisq import BS5
from coefficientFinder import getCoefficient
from massAction import massAction
from Hill import HillEquation
from indirect import indirect
from noise import physicalNoise, dropOut, dataSelector


def simulationODEs(GenesDict: dict,
                   ODEs: list | tuple | str,
                   T: tuple,
                   Coeff: dict | None = None,
                   **kwargs) -> None:
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
    stochasticNoiseAmplitude = kwargs.get("stochasticNoiseAmplitude", 0.2)
    physicalNoiseAmplitude = kwargs.get("physicalNoiseAmplitude", 0.2)
    normalisationBool = kwargs.get("normalisationBool", False)
    dropOutRate = kwargs.get("dropOutRate", 0.02)
    timeMeasuresNbMax = kwargs.get("timeMeasuresNbMax", 200)

    def otherODE(L: list) -> bool:
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
        stochasticNoiseAmplitudeMA = kwargs.get("stochasticNoiseAmplitudeMA",
                                                stochasticNoiseAmplitude)
        physicalNoiseAmplitudeMA = kwargs.get("physicalNoiseAmplitudeMA",
                                              physicalNoiseAmplitude)
        normalisationBoolMA = kwargs.get("normalisationBoolMA",
                                         normalisationBool)
        dropOutRateMA = kwargs.get("dropOutRateMA", dropOutRate)
        timeMeasuresNbMaxMA = kwargs.get("timeMeasuresNbMaxMa",
                                         timeMeasuresNbMax)
        K = [Coeff["TranscriptionsRate"][i]/Coeff["mRNAAvg"][i]
             for i in range(genesNb)]

        def equation(t: float, G: np.ndarray) -> np.ndarray:
            return massAction(t, G, M, K, stochasticNoiseAmplitudeMA)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5, method=BS5)
        dataY, dataX = dataSelector(solution.y, solution.t,
                                    timeMeasuresNbMaxMA)
        physicalNoise(dataY, physicalNoiseAmplitudeMA)
        if normalisationBoolMA:
            normalisation(dataY)
        dropOut(dataY, dropOutRateMA)
        GenesDict["massActionY"] = dataY
        GenesDict["massActionX"] = dataX

    if "Hill" in ODEs:
        stochasticNoiseAmplitudeHI = kwargs.get("stochasticNoiseAmplitudeHI",
                                                stochasticNoiseAmplitude)
        physicalNoiseAmplitudeHI = kwargs.get("physicalNoiseAmplitudeHI",
                                              physicalNoiseAmplitude)
        normalisationBoolHI = kwargs.get("normalisationBoolHI",
                                         normalisationBool)
        dropOutRateHI = kwargs.get("dropOutRateHI", dropOutRate)
        timeMeasuresNbMaxHI = kwargs.get("timeMeasuresNbMaxHI",
                                         timeMeasuresNbMax)
        K = Coeff["TranscriptionsRate"]
        Kdeg = Coeff["mRNAsDeg"]

        def equation(t: float, G: np.ndarray) -> np.ndarray:
            return HillEquation(t, G, M, K, G0,
                                [0]*genesNb, Kdeg, 2,
                                stochasticNoiseAmplitudeHI)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5, method=BS5)
        dataY, dataX = dataSelector(solution.y, solution.t,
                                    timeMeasuresNbMaxHI)
        physicalNoise(dataY, physicalNoiseAmplitudeHI)
        if normalisationBoolHI:
            normalisation(dataY)
        dropOut(dataY, dropOutRateHI)
        GenesDict["HillY"] = dataY
        GenesDict["HillX"] = dataX

    if "indirect" in ODEs:
        stochasticNoiseAmplitudeIN = kwargs.get("stochasticNoiseAmplitudeIN",
                                                stochasticNoiseAmplitude)
        physicalNoiseAmplitudeIN = kwargs.get("physicalNoiseAmplitudeIN",
                                              physicalNoiseAmplitude)
        normalisationBoolIN = kwargs.get("normalisationBoolIN",
                                         normalisationBool)
        dropOutRateIN = kwargs.get("dropOutRateIN", dropOutRate)
        timeMeasuresNbMaxIN = kwargs.get("timeMeasuresNbMaxIN",
                                         timeMeasuresNbMax)
        kTranslation = Coeff["TranslationsRate"]
        mRNAavg = G0
        kTranscription = Coeff["TranscriptionsRate"]
        degP = Coeff["ProtsDeg"]
        degmRNA = Coeff["mRNAsDeg"]
        Pavg = Coeff["ProtAvg"]

        def equation(t: float, G: np.ndarray) -> np.ndarray:
            mRNA = G[:genesNb]
            P = G[genesNb:]
            return indirect(t, mRNA, P, M, kTranscription,
                            kTranslation, degP, degmRNA,
                            Pavg, mRNAavg, 2, stochasticNoiseAmplitudeIN)
        G0_indirect = np.concatenate((G0, Coeff["ProtAvg"]))
        solution = solve_ivp(equation, [t0, tf], G0_indirect,
                             max_step=0.5, method=BS5)
        dataY, dataX = dataSelector(solution.y[:genesNb], solution.t,
                                    timeMeasuresNbMaxIN)
        dataProt = dataSelector(solution.y[genesNb:], solution.t,
                                timeMeasuresNbMaxIN)[0]
        physicalNoise(dataY, physicalNoiseAmplitudeIN)
        if normalisationBoolIN:
            normalisation(dataY)
        dropOut(dataY, dropOutRateIN)
        GenesDict["indirectY"] = dataY
        GenesDict["indirectX"] = dataX
        GenesDict["indirectProt"] = dataProt


def normalisation(YLists: np.ndarray) -> None:
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


def main():
    pass


if __name__ == "__main__":
    main()
