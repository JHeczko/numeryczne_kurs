import numpy as np
import matplotlib.pyplot as plt
import time

N = 1000
h = 0.01

#PLOT FUNCTION
def plot(arr_y,arr_y1,title = "Graficzne rozwiazanie"):
    arr_x = np.array([i*h for i in range(0,N)])
    plt.title(title)
    plt.grid(True)
    plt.plot(arr_x,arr_y1)
    plt.plot(arr_x,arr_y)
    plt.legend(["LU","Jawnie"])
    plt.show()

#CREATE TRIDIAGONAL MATRIX
def tridiag(matrixA_down_diag, matrixA_diag, matrixA_up_diag, k1=-1, k2=0, k3=1):
    return np.diag(matrixA_down_diag, k1) + np.diag(matrixA_diag, k2) + np.diag(matrixA_up_diag, k3)

def jawnieSolve(d):
    H = np.float64(h*h-2)
    for i in range(3,N):
        d[i] = -d[i-2]-(d[i-1]*H)
    return d

#setup for matrixA for LU
matrixA_down_diag = np.ones(N-1,np.float64)
matrixA_down_diag[N-2] = 0

matrixA_diag = np.array([(h*h)-2 for i in range(0,N)],np.float64)
matrixA_diag[0] = 1
matrixA_diag[N-1] = 0

matrixA_up_diag = np.ones(N-1, np.float64)
matrixA_up_diag[0] = 0
matrixA = tridiag(matrixA_down_diag, matrixA_diag, matrixA_up_diag)
matrixA[N-1][0] = -3
matrixA[N-1][1] = 4
matrixA[N-1][2] = 1
d1 = np.zeros(N, np.float64)
d1[0] = 1

#setup for jawnie
d = np.zeros(N, np.float64)
d[0] = 1
d[1] = np.float64(0.66667778)
d[2] = np.float64(0.33328889)

#Measuring time for jawnie solve
start = time.time()
y_jawnie = jawnieSolve(d)
end = time.time()

#Measuing time for Lu solve
start1 = time.time()
y_lu = np.linalg.solve(matrixA,d1)
end1 = time.time()

#Solving and plotting
print(f"Czas rozwiazania jawnie: {end-start}\n\nCzas rozwazania zswyklym rozkladem LU przy uzyciu numpy: {end1-start1}\n\nWzglednie porównanie LU do morisona: {(end1-start1)/(end-start)}")
plot(y_jawnie,y_lu,"Jawnie, LU")
