import numpy as np
import numpy.linalg as lg
import time as time

#Parameters
N_given = 1000
N = N_given - 1 #Bo mamy juz dwa pierwsze rozwiazania!!!
H = np.float64(1e-4)+2
error = 1e-10

#Help function
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

#Thomas Solver
def TDMA(A,d):
    a = A[0]
    b = A[1]
    c = A[2]
    n = len(d)
    beta= np.zeros(n,float)
    gamma= np.zeros(n, float)
    x = np.zeros(n,float)
    
    beta[0] = -c[0]/b[0]
    gamma[0] = d[0]/b[0]

    for i in range(1,n-1):
        beta[i] = -(c[i]/(b[i] + a[i-1]*beta[i-1]))
    for i in range(1,n):
        gamma[i] = (d[i] - a[i-1]*gamma[i-1])/(b[i] + a[i-1]*beta[i-1])
    x[n-1] = gamma[n-1]
    for i in range(n-2,-1,-1):
        x[i] = gamma[i] + beta[i]*x[i+1]
    return x

#Gaus-Seidle
def GausSeidel(L,U,b):
    iteration = 0
    xn=np.zeros(N)
    while(True):
        newb = b - TriDot(U,xn)
        xn1 = TDMA(L,newb)
        if(lg.norm(xn1-xn) < 1e-10):
            print(f"Norma dla Gaus: {lg.norm(xn1)}")
            print(f"Liczba iteracji Gaus: {iteration}")
            return xn1
        else:
            iteration += 1
            xn = xn1

#Init Symetryczna-Dodatnio Okreslona
a = np.array([-1 for i in range(0,N-1)])
b = np.array([H for i in range(0,N)],np.float64)
b[N-1] = 2
c = np.array([-1 for i in range(0,N-1)])
d = np.zeros(N,np.float64)
d[0] = 1
d[N-1] = 1

#LowerWithDiagonal-Matrix and UpperNoDiagonal-Matrix
U = [np.zeros(N-1,np.float64), np.zeros(N,np.float64), c]
L = [a,b,np.zeros(N-1,np.float64)]

#Iteracje
print(f"Norma dla numpy: {lg.norm(lg.solve(TriDiag(a,b,c),d))}")
GausSeidel(L,U,d)
