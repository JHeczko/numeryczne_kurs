import numpy as np
import numpy.linalg as lg

#Parametry
N_given = 500
N = N_given - 2
error = np.float64(1e-10)
H = np.float64(1e-4)-2

#Make TriDiagonal matrix for numpy
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

#Is positive define?
def isPositiveDefined(x):
    return np.all(np.linalg.eigvals(x) > 0)

#Tridiagonal dot product
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

#Checking if matrix is positive defined
def isPositiveDefined(x):
    return np.all(np.linalg.eigvals(x) > 0)

#Macierz A musi byc dodatnio okreslona
def symetrycznieRozwiaz(matrixA, d, metoda):
    newA = np.matmul(matrixA, np.transpose(matrixA))
    print(f"Is positive define: {isPositiveDefined(newA)}")
    y = metoda(newA,d)
    x = np.transpose(matrixA).dot(y)
    return x

#Metoda Gradientow optymalizacja
def gradientSprzerzony(A,b):
    xn = np.zeros(N)
    rn = b - A.dot(xn)
    pn = rn
    xn1 = np.zeros(N) 
    rn1 = np.zeros(N) 
    pn1 = np.zeros(N) 
    global iteracje2
    iteracje2 = 0
    while(lg.norm(rn) > error):
        alpha = rn.dot(rn)/pn.dot(A.dot(pn))
        rn1 = rn - alpha*A.dot(pn)
        beta = rn1.dot(rn1)/rn.dot(rn)
        pn1 = rn1 + beta*pn
        xn1 = xn + alpha*pn
        xn = xn1
        pn = pn1
        rn = rn1
        iteracje2 += 1
    return xn

#Funckja robioca macierz
def giveMacierz():
    a = np.array([1 for i in range(0,N-1)])
    b = np.array([H for i in range(0,N)],np.float64)
    c = np.array([1 for i in range(0,N-1)])
    matrixA = TriDiag(a,b,c)
    for i in range(N-1,0,-1):
        matrixA[i-1]*=np.abs(matrixA[i][i])
        matrixA[i-1]+=matrixA[i]
    return matrixA*-1

#Macierz innit
matrixA = giveMacierz()
d = np.zeros(N,np.float64)
d[0] = -1*matrixA[1][1]

#Rozwiazanie
print(f"Norm for numpy{lg.norm(lg.solve(matrixA,d))}")
print(f"Norm for conquered gradients{lg.norm(symetrycznieRozwiaz(matrixA,d,gradientSprzerzony))}")
print(f"Ilosc iteracji: {iteracje2}")
