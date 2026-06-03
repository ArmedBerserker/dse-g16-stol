import numpy as np
import control as ctr
#from derivatives_parameters import *
from test_params import *

def C1SymGen(muc, dc,CZadot, KY2, Cmadot):
    C1Sym = np.zeros((4, 4))
    C1Sym[0, 0] = -2 * muc*dc
    C1Sym[1, 1] = (CZadot - 2 * muc) * dc
    C1Sym[2, 2] = -dc
    C1Sym[3, 3] = -2 * muc * KY2* dc 
    C1Sym[3, 1] = Cmadot * dc
    return C1Sym


def C2SymGen(CXu, CXa,CZ0,CXq,CZu,CZa,CX0,CZq,muc,Cmu,Cma,Cmq):
    C2Sym = np.zeros((4, 4))
    C2Sym[0, 0] = CXu
    C2Sym[0, 1] = CXa
    C2Sym[0, 2] = CZ0
    C2Sym[0, 3] = CXq
    C2Sym[1, 0] = CZu
    C2Sym[1, 1] = CZa
    C2Sym[1, 2] = -CX0
    C2Sym[1, 3] = CZq + 2 * muc
    C2Sym[2, 3] = 1
    C2Sym[3, 0] = Cmu
    C2Sym[3, 1] = Cma
    C2Sym[3, 3] = Cmq
    return C2Sym


def C3SymGen(CXde,CZde,Cmde):
    C3Sym = np.zeros((4, 1))
    C3Sym[0] = -CXde
    C3Sym[1] = -CZde
    C3Sym[3] = -Cmde
    return C3Sym


A_sym = -np.linalg.solve(C1SymGen(muc, dc,CZadot, KY2, Cmadot), C2SymGen(CXu, CXa,CZ0,CXq,CZu,CZa,CX0,CZq, muc, Cmu,Cma,Cmq)) 
B_sym = -np.linalg.solve(C1SymGen(muc, dc,CZadot, KY2, Cmadot),C3SymGen(CXde,CZde,Cmde))
C_sym = np.eye(4)
D_sym = np.zeros((4,1))


def C1AsymGen(CYbdot,mub,db,KX2,KXZ,Cnbdot,KZ2):
    C1Asym = np.zeros((4, 4))
    C1Asym[0, 0] = (CYbdot - 2 *mub) * db
    C1Asym[1, 1] = -0.5 * db
    C1Asym[2, 2] = -4 * mub * KX2 * db
    C1Asym[2, 3] = 4 * mub * KXZ * db
    C1Asym[3, 0] = Cnbdot * db
    C1Asym[3, 2] = 4 * mub * KXZ *db
    C1Asym[3, 3] = -4 * mub * KZ2 * db
    return C1Asym


def C2AsymGen(CYb,CL,CYp, CYr,mub,Clb,Clp,Clr,Cnb,Cnp, Cnr):
    C2Asym = np.zeros((4, 4))
    C2Asym[0, 0] = CYb
    C2Asym[0, 1] = CL
    C2Asym[0, 2] = CYp
    C2Asym[0, 3] = CYr - 4 * mub
    C2Asym[1, 2] = 1
    C2Asym[2, 0] = Clb
    C2Asym[2, 2] = Clp
    C2Asym[2, 3] = Clr
    C2Asym[3, 0] = Cnb
    C2Asym[3, 2] = Cnp
    C2Asym[3, 3] = Cnr
    return C2Asym


def C3AsymGen(CYda, CYdr,Clda, Cldr, Cnda, Cndr):
    C3Asym = np.zeros((4, 2))
    C3Asym[0, 0] = -CYda
    C3Asym[0, 1] = -CYdr
    C3Asym[2, 0] = -Clda
    C3Asym[2, 1] = -Cldr
    C3Asym[3, 0] = -Cnda
    C3Asym[3, 1] = -Cndr
    return C3Asym


A_Asym = -np.linalg.solve(C1AsymGen(CYbdot,mub,db,KX2,KXZ,Cnbdot,KZ2), C2AsymGen(CYb,CL,CYp, CYr,mub,Clb,Clp,Clr,Cnb,Cnp, Cnr)) 
B_Asym = -np.linalg.solve(C1AsymGen(CYbdot,mub,db,KX2,KXZ,Cnbdot,KZ2),C3AsymGen(CYda, CYdr,Clda, Cldr, Cnda, Cndr))
C_Asym = np.eye(4)
D_Asym = np.zeros((4,2))

sys_sym = ctr.ss(A_sym, B_sym, C_sym, D_sym)
sys_asym = ctr.ss(A_Asym, B_Asym, C_Asym, D_Asym)

# print(A_sym)
# print(A_Asym)
# print(B_Asym)
# print(B_sym)
# def ASymGen(cessna):
#     C1Sym = C1SymGen(cessna)
#     C2Sym = C2SymGen(cessna)
#     A = - np.matmul(np.linalg.inv(C1Sym), C2Sym)
#     return A


# def BSymGen(cessna):
#     C1Sym = C1SymGen(cessna)
#     C3Sym = C3SymGen(cessna)
#     B = -np.matmul(np.linalg.inv(C1Sym), C3Sym)
#     return B


# def CSymGen():
#     C = np.eye(4)
#     return C


# def DSymGen():
#     D = np.zeros((4, 1))
#     return D

# def AAsymGen(cessna):
#     C1Asym = C1AsymGen(cessna)
#     C2Asym = C2AsymGen(cessna)
#     A = -np.matmul(np.linalg.inv(C1Asym), C2Asym)
#     return A


# def BAsymGen(cessna):
#     C1Asym = C1AsymGen(cessna)
#     C3Asym = C3AsymGen(cessna)
#     B = -np.matmul(np.linalg.inv(C1Asym), C3Asym)
#     return B


# def CAsymGen():
#     return np.eye(4)


# def DAsymGen():
#     return np.zeros((4, 2))
