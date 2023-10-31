import matplotlib.pyplot as plt
import numpy as np

dolny_range = -55
gorny_range = -1
x_zmienna = np.float32(0.2)
h_zmienna = 2

def plot(arr_y,title,path):
    arr_x = []
    buf = -55.
    for x in range(0,len(arr_y)+1):
        if(x == 0):
            continue
        arr_x.append(2**buf)
        buf+=0.5
    fig, ax = plt.subplots(figsize=(10, 6))
    plt.grid(True)
    ax.loglog(arr_y)
    plt.yscale("log", base=10)
    # plt.xscale("log", base=10)
    ax.set(xlabel='h', ylabel=' e = log(e)',
        title=title)
    fig.savefig(path)
    plt.show()

def plot2(arr_y1, arr_y2):
    
    # Plotting both the curves simultaneously 
    plt.loglog(arr_y1, color='r', label='A') 
    plt.loglog(arr_y2, color='g', label='B') 
    plt.grid(True)
    
    # Naming the x-axis, y-axis and the whole graph 
    plt.xlabel("Angle") 
    plt.ylabel("Magnitude") 
    plt.title("Sine and Cosine functions") 
    
    # Adding legend, which helps us recognize the curve according to it's color 
    plt.legend() 
    
    # To load the display window 
    plt.show() 
    

file_data = []
file_data2 = []
file = open("./dane_iloraz1.txt", "r")
file2 = open("./dane_iloraz2.txt", "r")
for i in file:
    buf = np.float32(i[:-1])
    file_data.append(buf)
for i in file2:
    buf = np.float32(i[:-1])
    file_data2.append(buf)
plot2(file_data, file_data2)
