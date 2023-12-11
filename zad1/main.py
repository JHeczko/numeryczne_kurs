import numpy as np
def func(x):
    value1 = np.longdouble(0.0)
    for n in range(0,23):
        cos_val = np.cos(n*(x*x*x*x))
        value1 += np.longdouble( (np.exp(3*n) - np.exp(3*n)*cos_val*cos_val + cos_val) / np.exp(4*n))
    return value1
#tutaj koncze liczenie ilosci operacji, nie wiem jak dostac dokladna ilosc operacji, ktora wykonuje numpy, wiec zostawiam to jako optymalizacja, ktora uwzgeldnia kazda inna optymalizacje(bo dosc szybko liczy tę sume)
def func2(x):
    n = np.arange(0,23)
    wyrazenie = np.longdouble( (np.exp(3*n) - np.exp(3*n)*np.cos(n*(x**4))**2 + np.cos(n*(x**4))) / np.exp(4*n))
    return np.sum(wyrazenie)

print(func(1))
print(func2(1))