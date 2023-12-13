import time
import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt


N_INPUT = 1000
N = N_INPUT - 2
h=np.float64(0.01)

#PLOTTING
def plot(arr_y,title = "Graficzne rozwiazanie"):
    arr_x = np.array([i*h for i in range(0,N_INPUT)])
    plt.title(title)
    plt.grid(True)
    plt.plot(arr_x,arr_y)
    plt.show()

#MAKING TRIDIAGONAL MATRIX
def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

#THOMAS SOLVER MODIFIED
def TDMA(a,b,c,d):
    n = len(d)
    w= np.zeros(n-1,np.float64)
    g= np.zeros(n, np.float64)
    p = np.zeros(n+2,np.float64)
    p[0] = 1
    p[n+1] = 0
    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]
    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n] = g[n-1]
    for i in range(n-1,0,-1):
        p[i] = g[i-1] - w[i-1]*p[i+1]
    return p

#setup for thomas
a = np.ones(N-1, np.float64)
b = np.array([(h*h)-2 for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.zeros(N,np.float64)
d[0] = -1

#setup for matrix
a1 = np.ones(N_INPUT-1,np.float64)
a1[N_INPUT-2] = 0
b1 = np.array([(h*h)-2 for i in range(0,N_INPUT)],np.float64)
b1[N_INPUT-1] = 1
b1[0] = 1
c1 = np.ones(N_INPUT-1,np.float64)
c1[0] = 0
d1 = np.zeros(N_INPUT)
d1[0] = 1

matrixA = tridiag(a1,b1,c1)

#measure time for thomas
start = time.time()
y_thomas = TDMA(a,b,c,d)
end = time.time()
value = end - start

#measure time for LU
start1 = time.time()
y_lu = lg.solve(matrixA,d1)
end1 = time.time()
value1 = end1 - start1

#PRINTING AND PLOTTING
print(f"Solved matrixA with numpy no efficent with time: {value1}:\n")
print(f"Solved matrixA with thomas algorithm in time: {value}:\n")
print(f"Porownanie wzgledne czasu LU do czasu thomasa:\n{np.float64(value1/value)}")

plot(y_thomas)
plot(y_lu)