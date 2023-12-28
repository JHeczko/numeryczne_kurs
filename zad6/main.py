import numpy as np
import time as time

#Parameters
N = 8
error = 1e-10
powh2 = np.float64(0.0001)
gamma = 1/10
iteration = 0

#Help function
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

#RichardSon
def Richardson(matrixA,d,xn):
    xn1 = xn + gamma*(d - xn.dot(matrixA))
    return xn1
    

#Init
a = np.ones(N-1,np.float64)
a[N-2] = 0 #0
b = np.array([powh2+2 for i in range(0,N)],np.float64)
b[0] = 1
b[N-1] = 1
c = np.ones(N-1,np.float64)
c[0] = 0 # 0
d = np.zeros(N,np.float64)
d[0] = 1


print(np.linalg.solve(TriDiag(a,b,c),d))

xn=np.ones(N)
while(True):
    xn1 = Richardson(TriDiag(a,b,c),d,xn)
    if(xn[4]-xn1[4] < 1e-10):
        print(xn1)
        print(iteration)
        break
    else:
        iteration += 1
        xn = xn1

