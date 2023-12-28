import numpy as np
import time as time

#Parameters
N = 8
powh2 = np.float64(0.0001)
iteration = 0
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

#Jacobi
def Jacobi(D,R,b,xn):
    xn1 = TriDot(D,(b-TriDot(R,xn)))
    return xn1

#Init Symetryczna-Dodatnio Okreslona
a = np.ones(N-1,np.float64)
b = np.array([(powh2+2) for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.zeros(N,np.float64)
d[0] = 1

#Setup for Rest-Matrix and D_Inverse-Matrix
R = [a,np.zeros(N),c]
D_Inverse = [np.zeros(N-1), 1/b, np.zeros(N-1)]

#Porownainie
print(np.linalg.solve(TriDiag(a,b,c),d))

#Iteracja
xn=np.zeros(N)
while(True):
    xn1 = Jacobi(D_Inverse,R,d,xn)
    if(iteration == 300):
        print(iteration)
        print(xn1)
        break
    else:
        iteration += 1
        xn = xn1

