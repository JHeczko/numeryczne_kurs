import matplotlib.pyplot as plt
import numpy as np
import threading


#Parametry, mozemy je zmieniac w zaleznosci od tego od jakiej do jakiej potegi chcemy podnosic dana liczba, albo zmienic punkt w obliczanej pochodnej, albo zmienic podstawe wykladnicza h
type_val = "double" #typ ktory bedzie obslugiwany w calym programie "double" albo "float"
dolny_range = -55
gorny_range = -1

if(type_val == "double"):
    x_zmienna = np.float64(np.pi/2)
else:
    x_zmienna = np.float32(np.pi/2)

h_zmienna = 10 


def plot(arr_iloraz1, title,path):
    plt.title(title)
    plt.loglog(arr_iloraz1, color='r', label='A') 
    plt.grid(True)
    

    plt.xlabel("log^h") 
    plt.ylabel("|log^e-error|") 
    
    plt.legend() 
    plt.savefig(path)

    plt.show() 

def plot_multiple(arr_iloraz1, arr_iloraz2, arr_iloraz3, title,path):
    plot1 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
    plot2 = plt.subplot2grid((3, 3), (1, 0), colspan=3)
    plot3 = plt.subplot2grid((3, 3), (2, 0), colspan=3)

    plot1.grid(True)
    plot2.grid(True)
    plot3.grid(True)
    plot1.xaxis.set_tick_params(labelbottom=False)
    plot2.xaxis.set_tick_params(labelbottom=False)

    plt.title(title)
    plot1.set_title("Iloraz 1/A")
    plot2.set_title("Iloraz 2/B")
    plot3.set_title("Iloraz 3/C")

    plt.xlabel("log^h") 
    plt.ylabel("Error") 

    plot1.loglog(arr_iloraz1[0],arr_iloraz1[1], color="r")
    plot2.loglog(arr_iloraz2[0],arr_iloraz2[1], color="g")
    plot3.loglog(arr_iloraz3[0],arr_iloraz3[1], color="b")

    plt.savefig(path)
    plt.show()

def iloraz1(f, x, h):
    if(type_val == "double"):
        iloraz = np.float64((f(x + h) - f(x))/h)
    else:
        iloraz = np.float32((f(x + h) - f(x))/h)
    return iloraz

def iloraz2(f, x, h):
    if(type_val == "double"):
        iloraz = np.float64((f(x + h) - f(x - h)) / (2*h))
    else:
        iloraz = np.float32((f(x + h) - f(x - h)) / (2*h))
   # print(f"Value of second: {iloraz}")
    return iloraz

def iloraz3(f,x,h):
    if(type_val == "double"):
        iloraz = np.float64((-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12*h))
    else:
        iloraz = np.float32((-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12*h))
   # print(f"Value of third: {iloraz}")
    return iloraz

def funckja(x):
    if(type_val == "double"):
        value = np.float64(np.sin(x))
    else:
        value = np.float32(np.sin(x))
    #print(f"Value of first: {value}")
    return value

def pochodna_funckji(x):
    if(type_val == "double"):
        value = np.float64(np.cos(x))
    else:
        value = np.float32(np.cos(x))
    print(value)
    return value

def blad_przyblizenia(iloraz):
    arr = []
    arr_x = []
    arr_blad = []
    for e in np.arange(-16,-1,0.01):
        if(e == 0):
            continue
        blad = 2*pow(10,-17)*0.999999999999/pow(h_zmienna,e) + (1/2)*pow(h_zmienna,e)
        arr_x.append(blad)
    #print(arr_blad)
    arr.append(arr_x)
    return arr_x

iloraz1arr = blad_przyblizenia(iloraz1)

print(iloraz1arr)

plot(iloraz1arr, f"Ilorazy dla {type_val}", f"iloraz_{type_val}.png")

# plot_multiple(iloraz1arr, iloraz2arr, iloraz3arr, f"Ilorazy dla {type_val}", f"iloraz_multiple_{type_val}.png")
# plot(iloraz1arr, iloraz2arr, iloraz3arr, f"Ilorazy dla {type_val}", f"iloraz_{type_val}.png")