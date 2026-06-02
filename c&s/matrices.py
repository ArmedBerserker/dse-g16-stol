import numpy as np


def C1SymGen(cessna):
    cV = cessna.c / cessna.V0
    C1Sym = np.zeros((4, 4))
    C1Sym[0, 0] = -2 * cessna.muc / cessna.V0 * cV
    C1Sym[1, 1] = (cessna.CZadot - 2 * cessna.muc) * cV
    C1Sym[2, 2] = -cV
    C1Sym[3, 3] = -2 * cessna.muc * cessna.KY2 * cV ** 2
    C1Sym[3, 1] = cessna.Cmadot * cV
    return C1Sym


def C2SymGen(cessna):
    cV = cessna.c / cessna.V0
    C2Sym = np.zeros((4, 4))
    C2Sym[0, 0] = cessna.CXu / cessna.V0
    C2Sym[0, 1] = cessna.CXa
    C2Sym[0, 2] = cessna.CZ0
    C2Sym[0, 3] = cessna.CXq * cV
    C2Sym[1, 0] = cessna.CZu / cessna.V0
    C2Sym[1, 1] = cessna.CZa
    C2Sym[1, 2] = -cessna.CX0
    C2Sym[1, 3] = (cessna.CZq + 2 * cessna.muc) * cV
    C2Sym[2, 3] = cV
    C2Sym[3, 0] = cessna.Cmu / cessna.V0
    C2Sym[3, 1] = cessna.Cma
    C2Sym[3, 3] = cessna.Cmq * cV
    return C2Sym


def C3SymGen(cessna):
    C3Sym = np.zeros((4, 1))
    C3Sym[0] = cessna.CXde
    C3Sym[1] = cessna.CZde
    C3Sym[3] = cessna.Cmde
    return C3Sym


def ASymGen(cessna):
    C1Sym = C1SymGen(cessna)
    C2Sym = C2SymGen(cessna)
    A = - np.matmul(np.linalg.inv(C1Sym), C2Sym)
    return A


def BSymGen(cessna):
    C1Sym = C1SymGen(cessna)
    C3Sym = C3SymGen(cessna)
    B = -np.matmul(np.linalg.inv(C1Sym), C3Sym)
    return B


def CSymGen():
    C = np.eye(4)
    return C


def DSymGen():
    D = np.zeros((4, 1))
    return D


def C1AsymGen(cessna):
    bV = cessna.b / cessna.V0
    C1Asym = np.zeros((4, 4))
    C1Asym[0, 0] = (cessna.CYbdot - 2 * cessna.mub) * bV / cessna.V0
    C1Asym[1, 1] = -0.5 * bV
    C1Asym[2, 2] = -4 * cessna.mub * cessna.KX2 * bV ** 2 * 0.5
    C1Asym[2, 3] = 4 * cessna.mub * cessna.KXZ * bV ** 2 * 0.5
    C1Asym[3, 0] = cessna.Cnbdot * bV
    C1Asym[3, 2] = 4 * cessna.mub * cessna.KXZ * bV ** 2 * 0.5
    C1Asym[3, 3] = -4 * cessna.mub * cessna.KZ2 * bV ** 2 * 0.5
    return C1Asym


def C2AsymGen(cessna):
    bV = cessna.b / cessna.V0
    C2Asym = np.zeros((4, 4))
    C2Asym[0, 0] = cessna.CYb / cessna.V0
    C2Asym[0, 1] = cessna.CL
    C2Asym[0, 2] = cessna.CYp * bV * 0.5
    C2Asym[0, 3] = (cessna.CYr - 4 * cessna.mub) * bV * 0.5
    C2Asym[1, 2] = bV * 0.5
    C2Asym[2, 0] = cessna.Clb / cessna.V0
    C2Asym[2, 2] = cessna.Clp * bV * 0.5
    C2Asym[2, 3] = cessna.Clr * bV * 0.5
    C2Asym[3, 0] = cessna.Cnb / cessna.V0
    C2Asym[3, 2] = cessna.Cnp * bV * 0.5
    C2Asym[3, 3] = cessna.Cnr * bV * 0.5
    return C2Asym


def C3AsymGen(cessna):
    C3Asym = np.zeros((4, 2))
    C3Asym[0, 0] = cessna.CYda
    C3Asym[0, 1] = cessna.CYdr
    C3Asym[2, 0] = cessna.Clda
    C3Asym[2, 1] = cessna.Cldr
    C3Asym[3, 0] = cessna.Cnda
    C3Asym[3, 1] = cessna.Cndr
    return C3Asym


def AAsymGen(cessna):
    C1Asym = C1AsymGen(cessna)
    C2Asym = C2AsymGen(cessna)
    A = -np.matmul(np.linalg.inv(C1Asym), C2Asym)
    return A


def BAsymGen(cessna):
    C1Asym = C1AsymGen(cessna)
    C3Asym = C3AsymGen(cessna)
    B = -np.matmul(np.linalg.inv(C1Asym), C3Asym)
    return B


def CAsymGen():
    return np.eye(4)


def DAsymGen():
    return np.zeros((4, 2))
