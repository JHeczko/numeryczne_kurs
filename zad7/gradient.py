import numpy as np
import numpy.linalg as lg



#Parametry
N_given = 1000
N = N_given - 1
error = np.float64(1e-10)
H = np.float64(1e-4)+2

#Make TriDiagonal matrix for numpy
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

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

#Preconditioner solver with specific iterative method
def preconditioningSolve(matrixA, matrixD1, matrixD2, d,method):
    matrixAMod = np.matmul(np.matmul(matrixD1,matrixA),matrixD2)
    solvedA = method(matrixAMod,matrixD1.dot(d))
    solvedATrue = matrixD2.dot(solvedA)
    return solvedATrue
def gradientSprzerzony2(A,b):
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

#Metoda Gradientow optymalizacja
def gradientSprzerzony(A,b):
    xn = np.zeros(N)
    rn = b - TriDot(A,xn)
    pn = rn
    xn1 = np.zeros(N) 
    rn1 = np.zeros(N) 
    pn1 = np.zeros(N) 
    global iteracje
    iteracje = 0
    while(lg.norm(rn) > error):
        alpha = rn.dot(rn)/pn.dot(TriDot(A,pn))
        rn1 = rn - alpha*TriDot(A,pn)
        beta = rn1.dot(rn1)/rn.dot(rn)
        pn1 = rn1 + beta*pn
        xn1 = xn + alpha*pn
        xn = xn1
        pn = pn1
        rn = rn1
        iteracje += 1
    return xn

#Macierz innit
a = np.array([-1 for i in range(0,N-1)])
b = np.array([H for i in range(0,N)],np.float64)
b[N-1] = 2
c = np.array([-1 for i in range(0,N-1)])
d = np.zeros(N,np.float64)
d[0] = 1
d[N-1] = 1
matrixA = [a,b,c]

#Polepszanie wspolczynika uwarunkowania naszej macierzy, poprzez preconditioner solver
D1 = np.array([(i+1) for i in range(0,N)])
D2 = np.array([(i+1) for i in range(0,N)])
matrixD1 = TriDiag(np.zeros(N-1),D1,np.zeros(N-1))
matrixD2 = TriDiag(np.zeros(N-1),D2,np.zeros(N-1))

#Rozwiazanie
print(f"MatrixA is positive definied: {isPositiveDefined(TriDiag(a,b,c))}")
print(f"Norm of vector using numpy: {lg.norm(lg.solve(TriDiag(a,b,c),d))}")
print(f"Norm of vector using conqurence gradients: {lg.norm(gradientSprzerzony(matrixA,d))}")
print(f"Ilosc iteracji: {iteracje}")
print(f"Norm of vector using conqurence gradients and preconditioning: {lg.norm(preconditioningSolve(TriDiag(a,b,c), matrixD1, matrixD2, d, gradientSprzerzony2))}")
print(f"Ilosc iteracji: {iteracje2}")
