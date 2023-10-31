import math
import random
import numpy as np
from numpy import linalg as LA

zakresLosowania = np.float64(2/math.sqrt(5))

def kreski():
    print("-"*50)

#Funckja generujaca losowy wektor o bardzo zblizonej normie 10^(-6)
def generateErrorVector(vector):
    flag = True
    while(flag):
        for i in range(1,6):
            random_error = random.uniform(-zakresLosowania, zakresLosowania)
            random_error *= pow(10,-6)
            vector = np.append(vector, random_error)
        if((pow(10, -6)/LA.norm(vector)) >= 0.99999 and (pow(10, -6)/LA.norm(vector)) <= 1.00001):
            flag = False
        else:
            for j in range(1,6):
                vector = np.delete(vector, 0)
    return vector

#Ustawienia parametrow z zadania
macierzA = np.array([[ 2.554219275 , 0.871733993 , 0.052575899 , 0.240740262 , 0.316022841 ] ,
                     [ 0.871733993 , 0.553460938 , -0.070921727 , 0.255463951 , 0.707334556 ] ,
                     [ 0.052575899 , -0.070921727 , 3.409888776 , 0.293510439 , 0.847758171 ],
                     [ 0.240740262 , 0.255463951 , 0.293510439 , 1.108336850 , -0.206925123 ],
                     [0.316022841 , 0.707334556 , 0.847758171 , -0.206925123 , 2.374094162 ]], 
                     np.float64)

macierzB = np.array([
                    [ 2.645152285 , 0.544589368 , 0.009976745 , 0.327869824 , 0.424193304 ],
                    [ 0.544589368 , 1.730410927 , 0.082334875 , -0.057997220 , 0.318175706 ],
                    [ 0.009976745 , 0.082334875 , 3.429845092 , 0.252693077 , 0.797083832],
                    [ 0.327869824 , -0.057997220 , 0.252693077 , 1.191822050 , -0.103279098 ],
                    [ 0.424193304 , 0.318175706 , 0.797083832 , -0.103279098 , 2.502769647 ]], 
                    np.float64)

wektorB = np.array([-0.642912346 , -1.408195475 , 4.595622394 , -5.073473196 , 2.178020609], np.float64)

#Generowanie wektora zaburzeń
wektorError = np.array([], np.float64)
wektorError = generateErrorVector(wektorError)

answearsMatrixA = LA.solve(macierzA, wektorB)
answearsMatrixB = LA.solve(macierzB, wektorB)

wektorB_zaburzony = np.add(wektorB, wektorError)

answearsMatrixA_zaburzone = LA.solve(macierzA, wektorB_zaburzony)
answearsMatrixB_zaburzone = LA.solve(macierzB, wektorB_zaburzony)

differenceMacierzA = np.absolute(answearsMatrixA - answearsMatrixA_zaburzone)
differenceMacierzB = np.absolute(answearsMatrixB - answearsMatrixB_zaburzone)

kreski()
print(f"\nMacierz A: {answearsMatrixA} \n\nZaburzona Macierz A: {answearsMatrixA_zaburzone}\n")
kreski()
print(f"\nMacierz B: {answearsMatrixB} \n\nZaburzona Macierz B: {answearsMatrixB_zaburzone}\n")
kreski()
print(f"\nTak wyglada roznica odpowiedzi odpowiednio:\n\nMacierzy A: {differenceMacierzA}\nMacierzy B: {differenceMacierzB}\n")
kreski()
print(f"\nKappa macierzy A: {LA.cond(macierzA)} \n\nKappa Macierzy B: {LA.cond(macierzB)}\n")
kreski()




