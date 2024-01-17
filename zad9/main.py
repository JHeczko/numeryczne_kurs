import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt

Macierz_Rozmiar = 4
N = Macierz_Rozmiar
L = 100
h = (2*L)/(Macierz_Rozmiar-1)

def plot1(vec,title):
    plt.grid(True)

    plt.title(title)

    plt.plot(vec[0])
    plt.plot(vec[1])
    plt.plot(vec[2])
    plt.plot(vec[3])

    plt.savefig("Wykres_"+title+".png")
    plt.show()

def plot(vec):
    plot1 = plt.subplot2grid((4, 4), (0, 0), colspan=3)
    plot2 = plt.subplot2grid((4, 4), (1, 0), colspan=3)
    plot3 = plt.subplot2grid((4, 4), (2, 0), colspan=3)
    plot4 = plt.subplot2grid((4, 4), (3, 0), colspan=3)

    plot1.yscale = "log"
    plot2.yscale = "log"
    plot3.yscale = "log"
    plot4.yscale = "log"

    plot1.grid(True)
    plot2.grid(True)
    plot3.grid(True)
    plot4.grid(True)
    plot1.xaxis.set_tick_params(labelbottom=False)
    plot2.xaxis.set_tick_params(labelbottom=False)
    plot3.xaxis.set_tick_params(labelbottom=False)

    plt.title("Wykres") 

    # plot1.loglog(vec[0], "o")
    # plot2.loglog(vec[1], "o")
    # plot3.loglog(vec[2], "o")
    # plot4.loglog(vec[3], "o")

    plot1.loglog(vec[0])
    plot2.loglog(vec[1])
    plot3.loglog(vec[2])
    plot4.loglog(vec[3])

    plt.savefig("wykres.png")
    plt.show()


#Make tridiagonal
def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

#CLASSIC THOMAS SOLVER
def TDMA(a,b,c,d):
    n = len(d)
    beta= np.zeros(n,float)
    gamma= np.zeros(n, float)
    x = np.zeros(n,float)
    
    beta[0] = -c[0]/b[0]
    gamma[0] = d[0]/b[0]

    for i in range(1,n-1):
        beta[i] = -(c[i]/(b[i] + a[i-1]*beta[i-1]))
    for i in range(1,n):
        gamma[i] = (d[i] - a[i-1]*gamma[i-1])/(b[i] + a[i-1]*beta[i-1])
    x[n-1] = gamma[n-1]
    for i in range(n-2,-1,-1):
        x[i] = gamma[i] + beta[i]*x[i+1]
    return x

def potegowa_odwrotna(a,b,c):
    e = [np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)]
    vk = np.zeros(N)
    wk = np.zeros(N)
    v_poprzednie = np.zeros(N)

    wartoscWlasnaA = 0
    wartoscWlasnaA_poprzednia = 0

    tabWartosciWlasnychA = np.zeros(4)

    iteration = 0
    for i in range(0,4):
        wartoscWlasnaA = 0
        wartoscWlasnaA_poprzednia = 0
        v_poprzednie = np.random.rand(N)
        v_poprzednie/=lg.norm(v_poprzednie)
        vk = np.zeros(N)
        wk = np.zeros(N)

        while(True):
            wartoscWlasnaA_poprzednia = wartoscWlasnaA

            wk = TDMA(a,b,c, v_poprzednie)

            vk = wk/lg.norm(wk)

            wartoscWlasnaA = (lg.norm(wk)/lg.norm(v_poprzednie))

            help = np.zeros(N)
            for j in range(0, i):
                iloczyn = e[j].dot(vk)
                help += e[j]*iloczyn
            vk -= help
            v_poprzednie = vk
            # print("-------------")
            # print(f"WLASNE WEKT: {e[0]}")
            # print(f"LAMBDA: {tabWartosciWlasnychA[0]}")

            if(np.abs(wartoscWlasnaA - wartoscWlasnaA_poprzednia) < 1e-14):
                tabWartosciWlasnychA[i] = 1/wartoscWlasnaA
                e[i] = v_poprzednie
                break

            tabWartosciWlasnychA[i] = wartoscWlasnaA
            e[i] = v_poprzednie

            iteration += 1
    print(iteration)
    return e,tabWartosciWlasnychA

def potegowa_odwrotna2(a,b,c):
    e = [np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)]
    vk = np.zeros(N)
    wk = np.zeros(N)

    v_poprzednie = np.zeros(N)
    v_poprzednie_poprzednie = np.zeros(N)

    wartoscWlasnaA = 0
    wartoscWlasnaA_poprzednia = 0
    wartoscWlasnaA_poprzednia_poprzednia = 0

    tabWartosciWlasnychA = np.zeros(4)

    iteration = 0

    for i in range(0,4):
        wartoscWlasnaA = 0
        wartoscWlasnaA_poprzednia = 3
        wartoscWlasnaA_poprzednia_poprzednia = 7
        v_poprzednie = np.random.rand(N)
        v_poprzednie/=lg.norm(v_poprzednie)
        v_poprzednie_poprzednie = np.zeros(N)
        vk = np.zeros(N)
        wk = np.zeros(N)

        while(True):
            wartoscWlasnaA_poprzednia_poprzednia = wartoscWlasnaA_poprzednia
            wartoscWlasnaA_poprzednia = wartoscWlasnaA



            wk = TDMA(a,b,c, v_poprzednie)

            vk = wk/lg.norm(wk)

            wartoscWlasnaA = (lg.norm(wk)/lg.norm(v_poprzednie))

            help = np.zeros(N)
            for j in range(0, i):
                iloczyn = e[j].dot(vk)
                help += e[j]*iloczyn
            vk -= help
            v_poprzednie = vk

            # print("-------------")
            # print(f"WLASNE WEKT: {e[0]}")
            # print(f"LAMBDA: {tabWartosciWlasnychA[0]}")

            q = np.abs(wartoscWlasnaA_poprzednia - wartoscWlasnaA) / np.abs(wartoscWlasnaA_poprzednia_poprzednia - wartoscWlasnaA_poprzednia)
            if(q*np.abs(wartoscWlasnaA - wartoscWlasnaA_poprzednia)/(q-1) < 1e-14):
                tabWartosciWlasnychA[i] = 1/wartoscWlasnaA
                e[i] = v_poprzednie
                break

            tabWartosciWlasnychA[i] = wartoscWlasnaA
            e[i] = v_poprzednie

            iteration += 1
    print(iteration)
    return e,tabWartosciWlasnychA

#Init gdzie A[0] = downDiag; A[1] = diag; A[2] = upDiag 
def x_zmienna(n): return -L+(n*h)

a = np.array([-1/(h*h) for i in range(0,Macierz_Rozmiar-1)])
b = np.array([(2/(h*h)) + 4 - (6/(np.cos(x_zmienna(i+1))**2)) for i in range(0,Macierz_Rozmiar)])
c = np.array([-1/(h*h) for i in range(0,Macierz_Rozmiar-1)])

values_true, vectors_true = lg.eig(tridiag(a,b,c))
vectors, values = potegowa_odwrotna(a,b,c)

print(vectors_true)
print(vectors)

plot1(vectors,"cosh")
plot1(vectors_true,"cosh_true")

