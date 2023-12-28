
import numpy as np
import time as time

#Parameters
N = 8
error = 1e-10
powh2 = np.float64(0.0001)
gamma = 1/2
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
def Richardson(matrixA,d,xn):
    xn1 = xn + gamma*(d - TriDot(matrixA,xn))
    return xn1
    

#Init Symetryczna-Dodatnio Okreslona
a = np.ones(N-1,np.float64)
b = np.array([powh2+2 for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.zeros(N,np.float64)
d[0] = 1
matrixA = [a,b,c]

#Porownanie
print(np.linalg.solve(TriDiag(a,b,c),d))

#Iteracje
xn=np.zeros(N)
while(True):
    xn1 = Richardson(matrixA,d,xn)
    if(np.abs(max(xn1) - max(xn)) < 1e-10):
        print(iteration)
        print(xn1)
        break
    else:
        iteration += 1
        xn = xn1

