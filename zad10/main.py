import numpy as np
import numpy.linalg as lg
import time

error = 1e-8
global iteration
iteration = 0
przedzialy = [[2,3],[3,5],[5,6]]

def sprawdzenie(x1, x2):
    if(np.abs(x1-x2) <= error):
        return True
    else:
        return False


def funkcja(x):
    return -x*x*x + 12*x*x - 46*x + 56

def bisekcja(f):
    miejsca_zerowe = [0,0,0]
    iteration = 0
    for i in range(0,3):
        x1 = przedzialy[i][0]
        x2 = przedzialy[i][1]
        while(True):
            x3 = (x1 + x2) / 2
            if(f(x3)*f(x2) < 0):
                x1 = x2
                x2 = x3
            if(f(x3)*f(x1) < 0):
                x1 = x1
                x2 = x3
            if(f(x3) == 0):
                miejsca_zerowe[i] = x3
                break
            if(sprawdzenie(x1,x2)):
                miejsca_zerowe[i] = x2
                break
            iteration += 1
    print(iteration)
    return miejsca_zerowe

def falsi(f):
    iteration = 0
    miejsca_zerowe = [0,0,0]
    for i in range(0,3):
        x1 = przedzialy[i][0]
        x2 = przedzialy[i][1]
        while(True):
            x3 = (f(x1)*x2 - f(x2)*x1) / (f(x1) - f(x2))
            if(f(x3)*f(x2) < 0):
                x1 = x3
                x2 = x2
            if(f(x3)*f(x1) < 0):
                x1 = x1
                x2 = x3
            if(f(x3) == 0):
                print("bo zero")
                miejsca_zerowe[i] = x3
                break
            if(sprawdzenie(x2,x1)):
                miejsca_zerowe[i] = x2
                break
            iteration += 1
    print(iteration)
    return miejsca_zerowe

def siecznych(f):
    iteration = 0
    miejsca_zerowe = [0,0,0]
    for i in range(0,3):
        x1 = przedzialy[i][0]
        x2 = przedzialy[i][1]
        while(True):
            x3 = (f(x1)*x2 - f(x2)*x1) / (f(x1) - f(x2))
            x1 = x2
            x2 = x3
            if(sprawdzenie(x1,x2)):
                miejsca_zerowe[i] = x2
                break
            iteration += 1
    print(iteration)
    return miejsca_zerowe

print("Bisekcja:",bisekcja(funkcja))
print("Falsi: ",falsi(funkcja))
print("Siecznych: ",siecznych(funkcja))