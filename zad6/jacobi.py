import numpy as np
import numpy.linalg as lg
import time as time

#Parameters
N_given = 1000
N = N_given - 1 #Bo mamy juz dwa pierwsze rozwiazania!!!
H = np.float64(1e-4)+2
error = 1e-10

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
def Jacobi(D,R,b):
    iteration = 0
    xn=np.zeros(N)
    while(True):
        xn1 = TriDot(D,(b-TriDot(R,xn))) # Tutaj jest funckja cala
        if(lg.norm(xn1-xn)< 1e-10):
            print(f"Norma dla Jacobi: {lg.norm(xn1)}")
            print(f"Liczba iteracji Jacobi: {iteration}")
            return xn1
        else:
            iteration += 1
            xn = xn1

#Init
a = np.array([-1 for i in range(0,N-1)])
b = np.array([H for i in range(0,N)],np.float64)
b[N-1] = 2
c = np.array([-1 for i in range(0,N-1)])
d = np.zeros(N,np.float64)
d[0] = 1
d[N-1] = 1

#Setup for Rest-Matrix and D_Inverse-Matrix
R = [a,np.zeros(N),c]
D_Inverse = [np.zeros(N-1), 1/b, np.zeros(N-1)]

#Iteracja
print(f"Norma dla numpy: {lg.norm(lg.solve(TriDiag(a,b,c),d))}")
Jacobi(D_Inverse,R,d)

