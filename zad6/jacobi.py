import numpy as np
import numpy.linalg as lg
import time as time

#Parameters
N_given = 100
N = N_given - 2 #Bo mamy juz dwa pierwsze rozwiazania!!!
H = np.float64(1e-4)-2
error = 1e-10
wypisz = 1
iteration = 0

#Make Tridiagonal Matrix
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

#Jacobi
def Jacobi(D,R,b,xn):
    xn1 = TriDot(D,(b-TriDot(R,xn)))
    return xn1

#Init Symetryczna-Dodatnio Okreslona
a = np.ones(N-1,np.float64)
b = np.array([H for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.zeros(N,np.float64)
d[0] = -1

#Setup for Rest-Matrix and D_Inverse-Matrix
R = [a,np.zeros(N),c]
D_Inverse = [np.zeros(N-1), 1/b, np.zeros(N-1)]

#Iteracja
xnhelp = np.ones(N)
xn=np.zeros(N)
while(True):
    xn1 = Jacobi(D_Inverse,R,d,xn)
    if(lg.norm(xn1-xnhelp)< 1e-10):
        if(wypisz):
            print(TriDiag(a,b,c))
            print(xn1)
            print(lg.solve(TriDiag(a,b,c),d))
        print(f"Liczba iteracji: {iteration}")
        break
    else:
        iteration += 1
        xnhelp = xn
        xn = xn1

