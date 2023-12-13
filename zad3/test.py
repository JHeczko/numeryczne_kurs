import time
import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def Thomas_Solver(a,b,c,d):
    n = len(d)
    for i in range(0,n-1):
        tmp = np.float64(b[i])
        tmp1 = np.float64(a[i])
        b[i]/=tmp
        c[i]/=tmp
        d[i]/=tmp
        d[i+1]/=tmp1
        if(i != n-2):
            c[i+1]/=tmp1
        b[i+1]/=tmp1
        a[i]/=tmp1
        a[i]-=b[i]
        b[i+1]-=c[i]
        d[i+1]-=d[i]
    d[n-1]/=d[n-1]
    b[n-1]/=b[n-1]
    for i in range(n-2,-1,-1):
        d[i]-= (d[i+1]*c[i])
        c[i]-= (b[i+1]*c[i])
    return d

# for i in 2:N
#     tmp = A[i, i-1] / A[i-1,i-1]
#     A[i,i] -= tmp * A[i-1,i]
#     b[i] -= tmp * b[i-1]
# end

# y = zeros(N+1)
# y[1] = 1
# y[end] = b[N] / A[N,N]

# for i in N-1:-1:1
#     y[i+1] = (b[i] - A[i,i+1] * y[i+2]) / A[i,i]
# end
N= 5
h=0.01

a = np.ones(N-1, np.float64)
b = np.array([2 for i in range(0,N)],np.float64)
c = np.ones(N-1,np.float64)
d = np.ones(N,np.float64)
d[0] = -1
print(tridiag(a,b,c))
print(np.linalg.solve(tridiag(a,b,c),d))
print(Thomas_Solver(a,b,c,d))