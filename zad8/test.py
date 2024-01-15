from optparse import Values
import numpy as np
import numpy.linalg as lg

N = 3

def potegowa(macierzA):
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
        while(True):
            wartoscWlasnaA_poprzednia = wartoscWlasnaA

            wk = macierzA.dot(v_poprzednie)

            vk = wk/lg.norm(wk)

            wartoscWlasnaA = lg.norm(wk)/lg.norm(v_poprzednie)

            help = np.zeros(N)
            for j in range(0, i):
                iloczyn = e[j].dot(vk)
                help += e[j]*iloczyn
            vk -= help
            v_poprzednie = vk

            # print("-------------")
            # print(f"WLASNE WEKT: {e}")
            # print(f"LAMBDA: {tabWartosciWlasnychA}")
            # time.sleep(1)

            # if(np.abs(wartoscWlasnaA - wartoscWlasnaA_poprzednia) < 1e-8):
            #     break
            time.sleep(0.5)

            tabWartosciWlasnychA[i] = wartoscWlasnaA
            e[i] = v_poprzednie
    return e,tabWartosciWlasnychA

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
            if np.abs(lambda_val_previous - lambda_val) < 1e-8:
                values[i] = lambda_val
                break
            
            lambda_val_previous = lambda_val
            e[i] = x
            x_poprzednie = x
    
    return e, values

# Przykład użycia
# Zakładamy macierz symetryczną
macierzA = np.array([[1,2,3],[2,4,5],[3,5,-1]])

e,values = rqi(macierzA)
print("Wartości własna:", e)
print("Wektor własny:", values)
print(lg.eig(macierzA))