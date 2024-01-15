import numpy as np
import numpy.linalg as lg
import time

N=3

def getColumn(A,column):
    a = [0,0,0]
    for i in range(0,3):
        a[i] = A[column][i]
    return np.array(a)

def rqi(macierzA):
    e = [np.zeros(N),np.zeros(N),np.zeros(N)]
    values = [0,0,0]
    for i in range(0,N):
        n = macierzA.shape[0]
        lambda_val = 0
        lambda_val_previous = 0
        #Inicjalizacja wektora losowego o dlugości 1
        x_poprzednie = np.random.rand(n)
        x_poprzednie /= np.linalg.norm(x_poprzednie)
        while(True):
            # Obliczanie wartości własnej przy użyciu ilorazu Rayleigh
            lambda_val = (x_poprzednie @ macierzA @ x_poprzednie) / (x_poprzednie @ x_poprzednie)
            
            # Rozwiązanie układu równań z (macierzA - lambda*I)x = x_poprzednie, aby uzyskać nowy wektor własny
            x = np.linalg.solve(macierzA - lambda_val * np.eye(N), x_poprzednie)
            
            # Normalizacja wektora własnego
            x /= np.linalg.norm(x)

            help = np.zeros(N)
            for j in range(0, i):
                iloczyn = e[j].dot(x)
                help += e[j]*iloczyn
            x -= help

            # Warunek stopu: jeśli norma różnicy wektorów własnych jest mniejsza niż tolerancja
            if (np.abs(lambda_val_previous - lambda_val) < 1e-8):
                values[i] = lambda_val
                break
            
            lambda_val_previous = lambda_val
            e[i] = x
            x_poprzednie = x
    
    return e, values

def potegowa(macierzA):
    znak = 1
    e = [np.zeros(N),np.zeros(N),np.zeros(N)]
    vk = np.zeros(N)
    wk = np.zeros(N)
    v_poprzednie = np.zeros(N)

    wartoscWlasnaA = 0
    wartoscWlasnaA_poprzednia = 0

    tabWartosciWlasnychA = np.zeros(N)

    for i in range(0,N):
        v_poprzednie = np.zeros(N)
        v_poprzednie[0] = 1
        wartoscWlasnaA = 0
        iteration = 0
        while(True):
            wartoscWlasnaA_poprzednia = wartoscWlasnaA

            wk = macierzA.dot(v_poprzednie)

            vk = wk/lg.norm(wk)

            wartoscWlasnaA = lg.norm(wk)/lg.norm(v_poprzednie)

            help = np.zeros(N)
            for j in range(0, i):
                iloczyn = e[j].dot(vk)
                help += e[j]*iloczyn
            if(v_poprzednie[0]/vk[0] < 0):
                znak = -1
            else:
                znak = 1
            vk -= help
            v_poprzednie = vk

            # print("-------------")
            # print(f"WLASNE WEKT: {e}")
            # print(f"LAMBDA: {tabWartosciWlasnychA}")
            # time.sleep(1)

            if(np.abs(wartoscWlasnaA - wartoscWlasnaA_poprzednia) < 1e-8):
                tabWartosciWlasnychA[i] = wartoscWlasnaA*znak
                e[i] = v_poprzednie
                break

            tabWartosciWlasnychA[i] = wartoscWlasnaA
            e[i] = v_poprzednie

            iteration += 1

    return e,tabWartosciWlasnychA

def QR(macierzA):
    iteration = 0
    Bn1 = macierzA
    Bn = np.array([[3,1,1],[1,7,1],[1,9,1]])
    C = np.eye(3)
    while(lg.norm(np.diag(Bn1,0) - np.diag(Bn,0)) > 1e-8):
        qr = lg.qr(Bn1)
        Q = qr.Q
        C = np.matmul(C,Q)
        Bn = Bn1
        Bn1 = np.matmul(qr.R,Q)
        iteration += 1
    return C,np.diag(Bn1,0)

#Init macierz A
macierzA = np.array([[1,2,3],[2,4,5],[3,5,-1]])
print(f"Wartość własna przy pomocy numpy: {lg.eigvals(macierzA)}")
print(f"Wartość własna przy pomocy potegowej metody: {potegowa(macierzA)[1]}")
print(f"Wartość własna przy pomocy QR: {QR(macierzA)[1]}")
print(f"Wartość własna przy pomocy RQI: {rqi(macierzA)[1]}")