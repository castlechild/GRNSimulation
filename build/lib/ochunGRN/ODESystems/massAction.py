#!/usr/bin/env python3
import numpy as np


def MaMatrice(M, K):
    MaPlus = np.abs(M) * K
    MaMinus = np.dot(-np.ones(np.shape(M)), MaPlus)
    MaMinus = np.diag(np.diag(MaMinus))
    #print("minus:",MaMinus)
    #print("plus:",MaPlus)
    return MaPlus + MaMinus

def massAction(t, G, A):
    del t
    return np.dot(A,G)

#####################################################################
def main():
    Ma = np.random.randint(2,size=(3,3))
    K = np.random.random((3,3))
    print(Ma)
    print(K)
    print(MaMatrice(Ma,K))


if __name__ == "__main__":
    main()