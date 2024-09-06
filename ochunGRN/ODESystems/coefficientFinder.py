#!/usr/bin/env python3

import random as rd
import numpy as np
import pandas as pd
import pkg_resources

excel_path = pkg_resources.resource_filename(
    "ochunGRN", "ODESystems/41586_2011_BFnature10098_MOESM304_ESM.xls")
document = pd.read_excel(excel_path)
Attribut = document.columns.tolist()
document = document.T


def getCoefficient(GenesNb: int) -> dict:
    """
    Extracts random coefficients for a specified number of genes.

        Parameters:
            - GenesNb (int): The number of genes to extract coefficients for.

        Returns:
            - dict: A dictionary containing the coefficients.
    """
    resDict = {
        "ProtsDeg": [],
        "mRNAsDeg": [],
        "TranscriptionsRate": [],
        "TranslationsRate": [],
        "mRNAAvg": [],
        "ProtAvg": []
    }

    for i in range(GenesNb):
        randomNb = rd.randint(0, 5027)
        L = document[randomNb].tolist()

        # Ensure that no selected values are NaN for required indices
        while True in pd.isna([L[i] for i in [13, 16, 19, 22, 25, 28]]):
            randomNb = rd.randint(0, 5027)
            L = document[randomNb].tolist()

        # Extract coefficients and calculate necessary values
        resDict["ProtsDeg"].append(np.log(2)/L[19])
        resDict["mRNAsDeg"].append(np.log(2)/L[22])
        resDict["TranscriptionsRate"].append(L[25])
        resDict["TranslationsRate"].append(L[28])
        resDict["mRNAAvg"].append(L[16])
        resDict["ProtAvg"].append(L[13])

    return resDict


def main():
    # i = rd.randint(0, 5027)
    # gene = document[i].tolist()
    # for j in range(len(gene)):
    #     print(Attribut[j],":\n",gene[j],"\n")
    print(getCoefficient(10))


if __name__ == "__main__":
    main()
