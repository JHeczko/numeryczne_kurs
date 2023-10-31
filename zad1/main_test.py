from math import exp
import numpy as np
import matplotlib.pyplot as plt

def plot(arr_x ,arr_y):
    fig, ax = plt.subplots(figsize=(10, 6))
    plt.grid(True)
    plt.yscale("linear")
    plt.xscale("linear")
    ax.plot(arr_x,arr_y)
    ax.set(xlabel='x', ylabel='y',
        title="Sum plot")
    fig.savefig("plot_sum.png")
    plt.show()

arr_x1 = []
arr_y1 =[]
arr_x2 = []
arr_y2 = []
file = open("./dane.txt","w")
file.truncate
file2 = open("./dane2.txt","w")
file2.truncate

def func2(x):
    value1 = np.longdouble(0.0)
    for n in range(0,25):
        cos_val = np.cos(n*(x**4))
        value1 += np.longdouble( (exp(3*n) - exp(3*n)*cos_val**2 + cos_val) / exp(4*n))
        arr_x2.append(n)
        arr_y2.append(value1)
        file2.write(str(value1))
        file2.write("\n")
def func3(x):
    value1 = np.longdouble(0.0)
    for n in range(0,25):
        cos_val = np.cos(n*(x**4))
        value1 += np.longdouble( (np.exp(3*n) - np.exp(3*n)*cos_val**2 + cos_val) / np.exp(4*n))
        arr_x1.append(n)
        arr_y1.append(value1)
        file.write(str(value1))
        file.write("\n")
def func1(x):
    value1 = np.longdouble(0.0)
    for n in range(0,50):
        value1 += np.longdouble( np.power(np.sin(n*pow(x,4)),2)*np.exp(-n) + np.cos(n*pow(x,4))*np.exp(-4*n) )
        arr_x1.append(n)
        arr_y1.append(value1)
        file.write(str(value1))
        file.write("\n")

# plot(arr_x ,func1(2), "cos", "cos.png")
# z(n,x) = (e^(3n) - e^(3n)*cos(n*(x^4))^2 + cos(n*(x^4))) / e^(4*n) , x>=0 , n>=0
# sum (e^(3*n)-e^(3*n)*(cos(n))^2 + cos(n)) / e^(4*n), n=0 to 25
# 1.4007813613216868712 latex sumf(1)

func2(1)
