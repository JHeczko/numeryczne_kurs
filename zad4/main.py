import numpy as np
import numpy.linalg as lg
import time
import matplotlib.pyplot as plt

N = 1000
h = 0.01

def plot(arr_y,title = "Graficzne rozwiazanie"):
    arr_x = np.array([i*h for i in range(0,N)])
    plt.title(title)
    plt.grid(True)
    plt.plot(arr_x,arr_y)
    plt.show()

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def TDMA(a,b,c,d):
    n = len(d)
    w= np.zeros(n-1,float)
    g= np.zeros(n, float)
    p = np.zeros(n,float)
    
    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]

    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p

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

print(matrixA)

start_lu = time.time()
y_lu = lg.solve(matrixA, d_matrix)
end_lu = time.time()

print(f"Solved matrixA with numpy no efficent with time {end_lu-start_lu}:\n")
print(f"Solved matrixA with thomas algorithm in time{end_thom-start_thom}:\n")
print(f"Porownanie wzgledne czasu LU do czasu thomasa:\n{np.float64((end_lu-start_lu)/(end_thom-start_thom))}")
plot(y_thomas,"Thomas")
plot(y_lu, "LU")