import numpy as np
from scipy.integrate import solve_ivp

from coefficientFinder import getCoefficient
from massAction import MaMatrice, massAction
from Hill import HillEquation
from indirect import indirect

def simulationODEs(GenesDict:dict, ODEs:list, T:tuple, Coeff:dict=None):
    def otherODE(L):
        for ode in L:
            if ode not in ["massAction", "Hill","indirect"]:
                return True
        return False 
    
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
        RatioCoeff = [Coeff["TranscriptionsRate"][i]/Coeff["mRNAAvg"][i] for i in range(genesNb)]
        K = np.resize(RatioCoeff,(genesNb,genesNb))
        Ma = MaMatrice(M, K)
        
        equation = lambda t,G: massAction(t, G, Ma)     
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        GenesDict["massActionY"] = solution.y
        GenesDict["massActionX"] = solution.t
    
    if "Hill" in ODEs:
        K = Coeff["TranscriptionsRate"]
        Kdeg = Coeff["mRNAsDeg"]
        equation = lambda t,G: HillEquation(t, G, M, K, G0, [0]*genesNb, Kdeg, 5)
        solution = solve_ivp(equation, [t0, tf], G0, max_step=0.5)
        GenesDict["HillY"] = solution.y
        GenesDict["HillX"] = solution.t

    if "indirect" in ODEs:
        k_P = Coeff["TranslationsRate"]
        Ka_P = G0
        k_mRNA = Coeff["TranscriptionsRate"]
        K_degP = Coeff["ProtsDeg"]
        K_degMRNA = Coeff["mRNAsDeg"]
        def equation(t, G):
            mRNA = G[:genesNb]
            P = G[genesNb:]
            return indirect(t, mRNA, P, M, k_mRNA, k_P, Ka_P, K_degP, K_degMRNA, 5)
        G0_indirect = np.concatenate((G0, Coeff["ProtAvg"]))
        solution = solve_ivp(equation, [t0, tf], G0_indirect, max_step= 0.5)
        GenesDict["indirectY"] = solution.y[:genesNb]
        GenesDict["indirectX"] = solution.t