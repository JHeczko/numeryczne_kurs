import numpy.linalg as lg
import numpy as np
import time as time

#Parameters
N = 6
H = np.float64(0.0001)+2
iteration = 0
gamma = 1/3
error = 1e-10

#Help function
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

#OverRelaxation No Optimalization
# def OverRelaxationNoOptimalization(A,b,xn):
#     xn1 = np.zeros(N)
#     for i in range(1,N+1):
#         sum1 = 0
#         sum2 = 0
#         for j in range(1,i):
#             sum1 += A[i-1][j-1]*xn1[j-1]
#         for j in range(i+1,N+1):
#             sum2 += A[i-1][j-1]*xn[j-1]
#         xn1[i-1] = (1-gamma)*xn[i-1]+(gamma/A[i-1][i-1])*(b[i-1] - sum1 - sum2)
#     return xn1

#OverRelaxation
def OverRelaxation(A,b,xn):
    xn1 = np.zeros(N)
    for i in range(1,N+1):
        sum1 = 0
        sum2 = 0
        if(i != 1):
            sum1+=A[i-1][i-2]*xn1[i-2]
        if(i != N):
            sum2 += A[i-1][i]*xn[i]
        xn1[i-1] = (1-gamma)*xn[i-1]+(gamma/A[i-1][i-1])*(b[i-1] - sum1 - sum2)
    return xn1

#Init Symetryczna-Dodatnio Okreslona
a = np.ones(N-1,np.float64)
b = np.array([H for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.zeros(N,np.float64)
d[0] = 1

#Gamma Setup
mi = max(lg.eig(np.diag(np.ones(N)) - np.diag(b).dot(TriDiag(a,b,c))).eigenvalues)
gamma = 1 + ((mi)/(1 + np.sqrt(1 - (mi*mi))))*((mi)/(1 + np.sqrt(1-(mi*mi))))

#Setup for A-Matrix
A = TriDiag(a,b,c)

#Porownanie
print(np.linalg.solve(A,d))

#Iteracja
xn=np.zeros(N)
while(True):
    xn1 = OverRelaxation(A,d,xn)
    if(np.abs(max(xn1) - max(xn)) < 1e-10):
        print(iteration)
        print(xn1)
        break
    else:
        iteration += 1
        xn = xn1