#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

//PARAMETRY PROGRAMU, TO MOZNA ZMIENIAC
#define dol_index -16 
#define gora_index -1
#define x_zmienna 0.2 //PUNKT x W KTOREJ LICZYMY POCHODNO
#define h_zmienna 10 //BAZA POTEGI DO KTOREJ PODNOSIMI h
using typ = double; //precyzja calego programu
//KONIEC PARAMETRÓW

using namespace std;

void write_file(vector<typ> arr, const char* nazwa){
    ofstream File(nazwa);
    for(int i = 0; i < arr.size(); ++i){
        File << arr[i] << '\n';
    }
    File.close();
}

void to_string(vector<typ> arr){
    cout << "[ ";
    for(int i = 0; i < arr.size(); ++i){
        cout << arr[i] << " , ";
    }
    cout << " ]";
}

typ funckja(typ x){
    typ value = sin(pow(x,2));
    return value;
}

typ pochodna_funckji(typ x){
    typ value = 2*x*cos(pow(x,2));
    return value;
}

typ iloraz1(typ (*f)(typ x),typ x, typ h){
    typ value = (f(x + h) - f(x));
    cout << "Pierwszy iloraz: " << value << endl;
    typ value2 = value / h;
    return value2;
}

typ iloraz2(typ (*f)(typ x),typ x, typ h){
    typ value = (f(x + h) - f(x - h));
    typ value2 = value / (2*h);
    return value2;
}

vector<typ> blad_przyblizenia(typ (*iloraz)(typ (*f)(typ x),typ x, typ h)){
    vector<typ> arr_blad;
    for(typ e = dol_index; e < gora_index; e += 0.01){
        if(e == 0){
            continue;
        }
        //cout << iloraz1(funckja, x_zmienna, pow(h_zmienna, e)) << endl;
        //cout << pochodna_funckji(x_zmienna) << endl;
        typ blad = abs(pochodna_funckji(x_zmienna) - iloraz(funckja, x_zmienna, pow(h_zmienna, e)));
        arr_blad.push_back(blad);
    }
    return arr_blad;
}

int main(){
    vector<typ> arr_iloraz1 = blad_przyblizenia(iloraz1);
    vector<typ> arr_iloraz2 = blad_przyblizenia(iloraz2);
    to_string(arr_iloraz1);
    to_string(arr_iloraz2);
    write_file(arr_iloraz1, "dane_iloraz1.txt");
    write_file(arr_iloraz2, "dane_iloraz2.txt");
    
}

// def blad_przyblizenia(iloraz):
//     arr_blad = []
//     for e in range(dolny_range, gorny_range):
//         if(e == 0):
//             continue
//         blad = np.abs(pochodna_funckji(x_zmienna) - iloraz(funckja, x_zmienna, pow(h_zmienna,e)))
//         arr_blad.append(blad)
//     #print(arr_blad)
//     return arr_blad

// def iloraz1(f, x, h):
//     iloraz = np.float32((f(x + h) - f(x))/h)
//     return iloraz

// def iloraz2(f, x, h):
//     iloraz = np.float32((f(x + h) - f(x - h)) / (2*h))
//     return iloraz

// def funckja(x):
//     value = np.float32(np.sin(pow(x,2)))
//     return value

// def pochodna_funckji(x):
//     value = np.float32(2*x*np.cos(pow(x,2)))
//     return value

// def blad_przyblizenia(iloraz):
//     arr_blad = []
//     for e in range(dolny_range, gorny_range):
//         if(e == 0):
//             continue
//         blad = np.abs(pochodna_funckji(x_zmienna) - iloraz(funckja, x_zmienna, pow(h_zmienna,e)))
//         arr_blad.append(blad)
//     #print(arr_blad)
//     return arr_blad