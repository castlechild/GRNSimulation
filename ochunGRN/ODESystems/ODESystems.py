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
            - GenesDict (dict): A dictionary containing information about
            the genes, including the adjacency matrix (`AdjMatrice`) and
            the number of genes (`genesNb`).
            - ODEs (list, tuple, str): A list, tuple, or string indicating
            the types of ODEs to simulate
            (e.g., "massAction", "Hill", "indirect").
            - T (tuple): A tuple of two values indicating the start and end
            of the simulation (t0, tf).
            - Coeff (dict, optional): A dictionary containing the coefficients
            for the simulation. If no coefficients are provided, they will
            be generated automatically.
            - **kwargs: Optional keyword arguments for configuring
            noise parameters, normalization, dropout rate, and the number
            of data points to select.

        Keyword Args (common to all ODEs):
            - stochasticNoiseAmplitude (float, optional): Amplitude
            of stochastic noise for all ODEs. Default is 0.2.
            - physicalNoiseAmplitude (float, optional): Amplitude of
            physical noise for all ODEs. Default is 0.2.
            - normalisationBool (bool, optional): Flag indicating whether
            to normalize results. Default is False.
            - dropOutRate (float, optional): Rate at which to drop out data
            points. Default is 0.02.
            - timeMeasuresNbMax (int, optional): Maximum number of time points
            to measure. Default is 200.

        Keyword Args (specific to "massAction" ODE):
            - stochasticNoiseAmplitudeMA (float, optional): Specific
            stochastic noise amplitude for "massAction".
            - physicalNoiseAmplitudeMA (float, optional): Specific physical
            noise amplitude for "massAction".
            - normalisationBoolMA (bool, optional): Specific normalization
            flag for "massAction".
            - dropOutRateMA (float, optional): Specific dropout rate for
            "massAction".
            - timeMeasuresNbMaxMA (int, optional): Specific maximum time
            points for "massAction".

        Keyword Args (specific to "Hill" ODE):
            - stochasticNoiseAmplitudeHI (float, optional): Specific
            stochastic noise amplitude for "Hill".
            - physicalNoiseAmplitudeHI (float, optional): Specific physical
            noise amplitude for "Hill".
            - normalisationBoolHI (bool, optional): Specific normalization
            flag for "Hill".
            - dropOutRateHI (float, optional): Specific dropout rate for
            "Hill".
            - timeMeasuresNbMaxHI (int, optional): Specific maximum time
            points for "Hill".

        Keyword Args (specific to "indirect" ODE):
            - stochasticNoiseAmplitudeIN (float, optional): Specific
            stochastic noise amplitude for "indirect".
            - physicalNoiseAmplitudeIN (float, optional): Specific physical
            noise amplitude for "indirect".
            - normalisationBoolIN (bool, optional): Specific normalization
            flag for "indirect".
            - dropOutRateIN (float, optional): Specific dropout rate for
            "indirect".
            - timeMeasuresNbMaxIN (int, optional): Specific maximum time
            points for "indirect".

        Raises:
            - ValueError: If the specified ODEs are invalid or if
            the time interval is incorrect.

        Returns:
            - None: The simulation results are stored in `GenesDict`
            under keys corresponding to each ODE type
            (e.g., 'massActionY', 'HillY', 'indirectY', etc.).
    """
    # Default parameter values
    stochasticNoiseAmplitude = kwargs.get("stochasticNoiseAmplitude", 0.2)
    physicalNoiseAmplitude = kwargs.get("physicalNoiseAmplitude", 0.2)
    normalisationBool = kwargs.get("normalisationBool", False)
    dropOutRate = kwargs.get("dropOutRate", 0.02)
    timeMeasuresNbMax = kwargs.get("timeMeasuresNbMax", 200)

    # Helper function to check for unsupported ODE types
    def otherODE(L: list) -> bool:
        """
        Check if the list of ODEs contains elements other than
        the supported models.

            Parameters:
                - L (list): The list of ODEs to check.

            Returns:
                bool: True if unsupported ODEs are detected, otherwise False.
        """
        for ode in L:
            if ode not in ["massAction", "Hill", "indirect"]:
                return True
        return False

    # Validate ODE list and time interval
    if isinstance(ODEs, str):
        ODEs = [ODEs]
    if len(ODEs) == 0 and otherODE(ODEs):
        raise ValueError("ODEs not valid")
    if len(T) != 2 and T[0] >= T[1]:
        raise ValueError("the temporal pair is invalid")

    # Extract number of genes and coefficients
    genesNb = GenesDict["genesNb"]
    if Coeff is None:
        # Generate coefficients if not provided
        Coeff = getCoefficient(genesNb)
    G0 = Coeff["mRNAAvg"]  # Initial mRNA levels

    # Adjacency matrix for gene interactions
    M = GenesDict["AdjMatrice"]
    t0, tf = T

    # Store coefficients and ODEs in `GenesDict`
    GenesDict["Coefficients"] = Coeff
    GenesDict["ODEs"] = ODEs

    # Simulation for massAction law
    if "massAction" in ODEs:
        # Specific parameters for massAction law
        stochasticNoiseAmplitudeMA = kwargs.get("stochasticNoiseAmplitudeMA",
                                                stochasticNoiseAmplitude)
        physicalNoiseAmplitudeMA = kwargs.get("physicalNoiseAmplitudeMA",
                                              physicalNoiseAmplitude)
        normalisationBoolMA = kwargs.get("normalisationBoolMA",
                                         normalisationBool)
        dropOutRateMA = kwargs.get("dropOutRateMA", dropOutRate)
        timeMeasuresNbMaxMA = kwargs.get("timeMeasuresNbMaxMa",
                                         timeMeasuresNbMax)

        # Calculate specific constants
        K = [Coeff["TranscriptionsRate"][i]/Coeff["mRNAAvg"][i]
             for i in range(genesNb)]

        # Define the differential equation for massAction law
        def equation(t: float, G: np.ndarray) -> np.ndarray:
            return massAction(t, G, M, K, stochasticNoiseAmplitudeMA)

        # Solve ODEs using `solve_ivp`
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5, method=BS5)
        dataY, dataX = dataSelector(solution.y, solution.t,
                                    timeMeasuresNbMaxMA)

        # Apply physical noise, normalization
        physicalNoise(dataY, physicalNoiseAmplitudeMA)
        if normalisationBoolMA:
            normalization(dataY)
        dropOut(dataY, dropOutRateMA)

        # Store results in `GenesDict`
        GenesDict["massActionY"] = dataY
        GenesDict["massActionX"] = dataX

    # Simulation for Hill law
    if "Hill" in ODEs:
        # Specific parameters for "Hill"
        stochasticNoiseAmplitudeHI = kwargs.get("stochasticNoiseAmplitudeHI",
                                                stochasticNoiseAmplitude)
        physicalNoiseAmplitudeHI = kwargs.get("physicalNoiseAmplitudeHI",
                                              physicalNoiseAmplitude)
        normalisationBoolHI = kwargs.get("normalisationBoolHI",
                                         normalisationBool)
        dropOutRateHI = kwargs.get("dropOutRateHI", dropOutRate)
        timeMeasuresNbMaxHI = kwargs.get("timeMeasuresNbMaxHI",
                                         timeMeasuresNbMax)

        # Specific constants for "Hill"
        K = Coeff["TranscriptionsRate"]
        Kdeg = Coeff["mRNAsDeg"]

        # Define the differential equation for "Hill"
        def equation(t: float, G: np.ndarray) -> np.ndarray:
            return HillEquation(t, G, M, K, G0,
                                [0]*genesNb, Kdeg, 2,
                                stochasticNoiseAmplitudeHI)

        # Solve ODEs using `solve_ivp`
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5, method=BS5)
        dataY, dataX = dataSelector(solution.y, solution.t,
                                    timeMeasuresNbMaxHI)

        # Apply noise, normalization, and dropout
        physicalNoise(dataY, physicalNoiseAmplitudeHI)
        if normalisationBoolHI:
            normalization(dataY)
        dropOut(dataY, dropOutRateHI)

        # Store results in `GenesDict`
        GenesDict["HillY"] = dataY
        GenesDict["HillX"] = dataX

    # Simulation for "indirect" ODE
    if "indirect" in ODEs:
        # Specific parameters for "indirect"
        stochasticNoiseAmplitudeIN = kwargs.get("stochasticNoiseAmplitudeIN",
                                                stochasticNoiseAmplitude)
        physicalNoiseAmplitudeIN = kwargs.get("physicalNoiseAmplitudeIN",
                                              physicalNoiseAmplitude)
        normalisationBoolIN = kwargs.get("normalisationBoolIN",
                                         normalisationBool)
        dropOutRateIN = kwargs.get("dropOutRateIN", dropOutRate)
        timeMeasuresNbMaxIN = kwargs.get("timeMeasuresNbMaxIN",
                                         timeMeasuresNbMax)

        # Specific constants for "indirect"
        kTranslation = Coeff["TranslationsRate"]
        mRNAavg = G0
        kTranscription = Coeff["TranscriptionsRate"]
        degP = Coeff["ProtsDeg"]
        degmRNA = Coeff["mRNAsDeg"]
        Pavg = Coeff["ProtAvg"]

        # Define the differential equation for "indirect"
        def equation(t: float, G: np.ndarray) -> np.ndarray:
            mRNA = G[:genesNb]
            P = G[genesNb:]
            return indirect(t, mRNA, P, M, kTranscription,
                            kTranslation, degP, degmRNA,
                            Pavg, mRNAavg, 2, stochasticNoiseAmplitudeIN)
        G0_indirect = np.concatenate((G0, Coeff["ProtAvg"]))

        # Solve ODEs using `solve_ivp`
        solution = solve_ivp(equation, [t0, tf], G0_indirect,
                             max_step=0.5, method=BS5)
        dataY, dataX = dataSelector(solution.y[:genesNb], solution.t,
                                    timeMeasuresNbMaxIN)
        dataProt = dataSelector(solution.y[genesNb:], solution.t,
                                timeMeasuresNbMaxIN)[0]

        # Apply noise, normalization, and dropout
        physicalNoise(dataY, physicalNoiseAmplitudeIN)
        if normalisationBoolIN:
            normalization(dataY)
        dropOut(dataY, dropOutRateIN)

        # Store results in `GenesDict`
        GenesDict["indirectY"] = dataY
        GenesDict["indirectX"] = dataX
        GenesDict["indirectProt"] = dataProt


def normalization(YLists: np.ndarray) -> None:
    """
    Normalize the lists of values by dividing each element by the maximum
    of its respective list.

        Parameters:
            - YLists (list of list): A list containing lists of values
            to normalize.
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
