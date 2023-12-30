import numpy as np
import numpy.linalg as lg

#Parametry
N_given = 100
N = N_given - 2
error = np.float64(1e-10)
H = np.float64(1e-4)-2

#Make TriDiagonal matrix for numpy
def TriDiag(a,b,c):
    return np.diag(a,-1) + np.diag(b,0) + np.diag(c,1)

#Checking if matrix is positive defined
def isPositiveDefined(x):
    return np.all(np.linalg.eigvals(x) > 0)

#Preconditioner solver with specific iterative method
def preconditioningSolve(matrixA, matrixD1, matrixD2, d,method):
    matrixAMod = np.matmul(np.matmul(matrixD1,matrixA),matrixD2)
    solvedA = method(matrixAMod,matrixD1.dot(d))
    solvedATrue = matrixD2.dot(solvedA)
    return solvedATrue

#Metoda Gradientow
def gradientSprzerzony(A,b):
    xn = np.zeros(N)
    rn = b - A.dot(xn)
    pn = rn
    xn1 = np.zeros(N) 
    rn1 = np.zeros(N) 
    pn1 = np.zeros(N) 
    global iteracje
    iteracje = 0
    while(lg.norm(rn) > error):
        alpha = rn.dot(rn)/pn.dot(A.dot(pn))
        rn1 = rn - alpha*A.dot(pn)
        beta = rn1.dot(rn1)/rn.dot(rn)
        pn1 = rn1 + beta*pn
        xn1 = xn + alpha*pn
        xn = xn1
        pn = pn1
        rn = rn1
        iteracje += 1
    return xn

#Macierz innit
a = np.array([1 for i in range(0,N-1)])
b = np.array([H for i in range(0,N)])
c = np.array([1 for i in range(0,N-1)])
d = np.zeros(N)
d[0] = 1
matrixA = TriDiag(a,b,c)

#Polepszanie wspolczynika uwarunkowania naszej macierzy
D1 = np.array([(i+1) for i in range(0,N)]) #D1 = D2
D2 = np.array([(i+1) for i in range(0,N)])
matrixD1 = TriDiag(np.zeros(N-1),D1,np.zeros(N-1))
matrixD2 = TriDiag(np.zeros(N-1),D2,np.zeros(N-1))

#Rozwiazanie
print(f"MatrixA is positive definied: {isPositiveDefined(matrixA)}")
print(f"Norm of vector using numpy: {lg.norm(lg.solve(matrixA,d))}")
print(f"Norm of vector using conqurence gradients: {lg.norm(preconditioningSolve(matrixA, matrixD1, matrixD2, d, gradientSprzerzony))}")
print(f"Ilosc iteracji: {iteracje}")
