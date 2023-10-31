import numpy as np
def func(x):
    value1 = np.longdouble(0.0)
    for n in range(0,25):
        cos_val = np.cos(n*(x**4))
        value1 += np.longdouble( (np.exp(3*n) - np.exp(3*n)*cos_val**2 + cos_val) / np.exp(4*n))
    return value1
print(func(1))