
import numpy as np
import numpy.linalg as lg
import time as time

#Parameters
N_given = 100
N = N_given - 2 #Bo mamy juz dwa pierwsze rozwiazania!!!
H = np.float64(1e-4)-2
gamma = 1/2
wypisz = 1
iteration = 0

#Help functions
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

#Dot product for Tridiagonal matricies
def TriDot(A,b):
    out = np.zeros(b.size)
    under = A[0]
    diag = A[1]
    up = A[2]
    out[0] = diag[0] * b[0] + b[1]*up[0]
    out[b.size-1] = diag[b.size-1]*b[b.size-1] + under[b.size-2]*b[b.size-2]
    for i in range(1, b.size-1):
        out[i] = under[i-1]*b[i-1] + diag[i]*b[i] + up[i]*b[i+1]
    return out

#RichardSonForTridagonal
def Richardson(matrixA,b,xn):
    xn1 = xn + gamma*(b - TriDot(matrixA,xn))
    return xn1
    

#Init Symetryczna-Dodatnio Okreslona
a = np.ones(N-1,np.float64)
b = np.array([H for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.zeros(N,np.float64)
d[0] = -1

#Setup for matrix
matrixA = [a,b,c]

#Gamma setup
help = lg.eig(TriDiag(a,b,c)).eigenvalues
gamma = 2/(min(help) + max(help))

#Iteracje
xn=np.zeros(N)
while(True):
    xn1 = Richardson(matrixA,d,xn)
    if(lg.norm(xn1-xn) < 1e-10):
        if(wypisz):
            print(TriDiag(a,b,c))
            print(xn1)
            print(lg.solve(TriDiag(a,b,c),d))
        print(f"Liczba iteracji: {iteration}")
        break
    else:
        iteration += 1
        xn = xn1

