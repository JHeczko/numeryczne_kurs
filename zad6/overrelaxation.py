import numpy.linalg as lg
import numpy as np
import time as time

#Parameters
N_given = 1000
N = N_given - 1 #Bo mamy juz dwa pierwsze rozwiazania!!!
H = np.float64(1e-4)+2
iteration = 0
gamma = 1/3
error = 1e-10

#Help function
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

# #OverRelaxation With No memory optimalization
# def OverRelaxation(A,b,xn):
#     xn1 = np.zeros(N)
#     for i in range(1,N+1):
#         sum1 = 0
#         sum2 = 0
#         if(i != 1):
#             sum1+=A[i-1][i-2]*xn1[i-2]
#         if(i != N):
#             sum2 += A[i-1][i]*xn[i]
#         xn1[i-1] = (1-gamma)*xn[i-1]+(gamma/A[i-1][i-1])*(b[i-1] - sum1 - sum2)
#     return xn1

#OverRelaxation 
def OverRelaxationHelp(A,d,xn):
    xn1 = np.zeros(N)
    a = A[0]
    b = A[1]
    c = A[2]
    for i in range(1,N+1):
        sum1 = 0
        sum2 = 0
        if(i != 1):
            sum1+=c[i-2]*xn1[i-2]
        if(i != N):
            sum2 += a[i-1]*xn[i]
        xn1[i-1] = (1-gamma)*xn[i-1]+(gamma/b[i-1])*(d[i-1] - sum1 - sum2)
    return xn1
def OverRelaxation(A,d):
    xn=np.zeros(N)
    iteration = 0
    while(True):
        xn1 = OverRelaxationHelp(A,d,xn)
        if(lg.norm(xn1-xn)< 1e-10):
            print(f"Norma dla numpy: {lg.norm(lg.solve(TriDiag(a,b,c),d))}")
            print(f"Norma dla OverRelaxation: {lg.norm(xn1)}")
            print(f"Liczba iteracji OverRelaxation: {iteration}")
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

#Gamma Setup
mi = max(lg.eig(np.diag(np.ones(N)) - np.diag(b).dot(TriDiag(a,b,c))).eigenvalues)
gamma = 1 + ((mi)/(1 + np.sqrt(1 - (mi*mi))))*((mi)/(1 + np.sqrt(1-(mi*mi))))

#Setup for A-Matrix
A = [a,b,c]

#Iteracja
print(gamma)
OverRelaxation(A,d)