from math import gamma
from re import X
from matplotlib.pylab import beta
import numpy as np
import numpy.linalg as lg
import time
import matplotlib.pyplot as plt

N = 1000
h = 0.01

def plot(arr_y,arr_y2,title = "Graficzne rozwiazanie"):
    arr_x = np.array([i*h for i in range(0,N)])
    plt.title(title)
    plt.grid(True)
    plt.plot(arr_x,arr_y)
    plt.plot(arr_x,arr_y2)
    plt.show()

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def TDMA(a,b,c,d):
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

#setup for thomas
a = np.ones(N-1, np.float64)
b = np.array([(h*h)-2 for i in range(0,N)], np.float64)
b[0] = 1
b[N-1] = -2
c = np.ones(N-1, np.float64)
c[0] = 0
d = np.zeros(N)
d[0] = 1
d[N-1] = -1

#setup for LU
d_matrix = np.zeros(N)
d_matrix[0] = 1
matrixA = tridiag(a,b,c)
matrixA[N-1][0] = 1

start_thom = time.time()
y_thomas = TDMA(a,b,c,d)
end_thom = time.time()

start_lu = time.time()
y_lu = lg.solve(matrixA, d_matrix)
end_lu = time.time()

print(f"Solved matrixA with numpy no efficent with time {end_lu-start_lu}:\n")
print(f"Solved matrixA with thomas algorithm in time{end_thom-start_thom}:\n")
print(f"Porownanie wzgledne czasu LU do czasu thomasa:\n{np.float64((end_lu-start_lu)/(end_thom-start_thom))}")
plot(y_thomas,y_lu,"Thomas")
