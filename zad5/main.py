import numpy as np
import matplotlib.pyplot as plt
import time

N = 1000
h = 0.01

#PLOT FUNCTION
def plot(arr_y,title = "Graficzne rozwiazanie"):
    arr_x = np.array([i*h for i in range(0,N)])
    plt.title(title)
    plt.grid(True)
    plt.plot(arr_x,arr_y)
    plt.show()

#CREATE TRIDIAGONAL MATRIX
def tridiag(matrixA_down_diag, matrixA_diag, matrixA_up_diag, k1=-1, k2=0, k3=1):
    return np.diag(matrixA_down_diag, k1) + np.diag(matrixA_diag, k2) + np.diag(matrixA_up_diag, k3)

#CLASSIC THOMAS SOLVER
def TDMA(a,b,c,d):
    n = len(d)
    w= np.zeros(n-1,np.float64)
    g= np.zeros(n, np.float64)
    p = np.zeros(n,np.float64)
    
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

#MORISON SOLVER WITH THOMAS
def morison_Solve(down_diag, diag, up_diag, b, u, v):
    z = TDMA(down_diag,diag,up_diag,b)
    y = TDMA(down_diag,diag,up_diag, u)
    x = z - (np.float64(y*np.dot(v,z))/np.float64(1+(np.dot(v,y))))
    return x

#setup for morison
matrixA_down_diag = np.ones(N-1,np.float64)
matrixA_down_diag[N-2] = 0

matrixA_diag = np.array([(h*h)-2 for i in range(0,N)],np.float64)
matrixA_diag[0] = 1
matrixA_diag[N-1] = 1

matrixA_up_diag = np.ones(N-1, np.float64)
matrixA_up_diag[0] = 0

d = np.zeros(N, np.float64)
d[0] = 1


#Morison setup
u = np.zeros(N,np.float64)
u[N-1] = 1

v = np.zeros(N,np.float64)
v[0] = -3
v[1] = 4
v[2] = 1
v[N-1] = -1


#Setup for LU
matrixA = tridiag(matrixA_down_diag, matrixA_diag, matrixA_up_diag)
matrixA += np.outer(u,v)

#Measuring time for morison solve with thomas
start = time.time()
y_morison = morison_Solve(matrixA_down_diag, matrixA_diag, matrixA_up_diag, d, u, v)
end = time.time()

#Measuing time for Lu solve
start1 = time.time()
y_lu = np.linalg.solve(matrixA,d)
end1 = time.time()

#Solving and plotting
print(f"Czas rozwiazanie thomasem i wzorem shermana morisona: {end-start}\nCzas rozwazania zswyklym rozkladem LU przy uzyciu numpy: {end1-start1}\nWzglednie porównanie LU do morisona: {(end1-start1)/(end-start)}")
plot(y_morison)
plot(y_lu)
