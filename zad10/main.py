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
        x3 = 0
        while(True):
            x3 = (x1 + x2) / 2
            if(sprawdzenie(x3,x2)):
                miejsca_zerowe[i] = x3
                break
            if(f(x3)*f(x2) <= 0):
                x1 = x2
                x2 = x3
            if(f(x3)*f(x1) <= 0):
                x1 = x1
                x2 = x3
            iteration += 1
    print("Bisekcja iteracje: ",iteration)
    return miejsca_zerowe

def falsi(f):
    iteration = 0
    miejsca_zerowe = [0,0,0]
    # x1_poprzednie = 0
    # x2_poprzednie = 0
    x3_poprzednie = 0
    for i in range(0,3):
        x1 = przedzialy[i][0]
        x2 = przedzialy[i][1]
        x3 = 30
        while(True):
            # x1_poprzednie = x1
            # x2_poprzednie = x2
            x3_poprzednie = x3
            x3 = (f(x1)*x2 - f(x2)*x1) / (f(x1) - f(x2))
            if(sprawdzenie(x3,x3_poprzednie)):
                miejsca_zerowe[i] = x3
                break
            if(f(x3)*f(x2) < 0):
                x1 = x3
                x2 = x2
            if(f(x3)*f(x1) < 0):
                x1 = x1
                x2 = x3
            if(f(x3) == 0):
                miejsca_zerowe[i] = x3
                break
            # if(sprawdzenie(x2,x2_poprzednie) & sprawdzenie(x1, x1_poprzednie)):
            #     miejsca_zerowe[i] = x2
            #     break
            iteration += 1
    print("Falsi iteracje: ",iteration)
    return miejsca_zerowe

def siecznych(f):
    iteration = 0
    miejsca_zerowe = [0,0,0]
    for i in range(0,3):
        x1 = przedzialy[i][0]
        x2 = przedzialy[i][1]
        while(True):
            x3 = (f(x1)*x2 - f(x2)*x1) / (f(x1) - f(x2))
            if(sprawdzenie(x2,x3)):
                miejsca_zerowe[i] = x3
                break
            if(iteration == 1000):
                print("Nie udalo sie znalezc miejsca zerowego")
                break
            x1 = x2
            x2 = x3
            iteration += 1
    print("Sieczne iteracje: ",iteration)
    return miejsca_zerowe


if __name__ == "__main__":
    print("Bisekcja:",bisekcja(funkcja))
    print("Falsi: ",falsi(funkcja))
    print("Siecznych: ",siecznych(funkcja))