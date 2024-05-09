from math import *
import numpy as np
import sys

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: 
        + Parametry planet: 
        """
      
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424518 
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "krasowskiego":
            self.a = 6378245.0
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ee = (2 * self.flat - self.flat ** 2)
        
    def Np(self, f):
        """ Funkcja oblicza promień krzywiznowy 
        Parameters 
        ----------
        f : FLOAT
                szerokosć geograficzna [radiany]
        
        Returns
        -------
        N : FLOAT
                promien krzywiznowy [metry]
        """
        
        N = self.a/np.sqrt(1 - self.ee * np.sin(f)**2)
        return(N)
        
  
    def xyz2plh(self, X_list, Y_list, Z_list):
        """" Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        """
        
        f = []
        l = []
        h = []

        for X, Y, Z in zip(X_list, Y_list, Z_list):
            p = np.sqrt(X**2 + Y**2)
            f = np.arctan(Z / (p * (1 - self.ee)))
         
        while True: 
             N = self.Np(f)
             h = (p / np.cos(f)) - N
             fs = f
             f = np.arctan(Z / (p * (1 - self.ee * (N / (N + h)))))
             
             if np.abs(fs - f) < (0.000001 / 206265):
                 break
                 
             l = np.arctan2(Y, X)
             f.append(f)
             l.append(l)
             h.append(h)

        return f,l,h 
  
        
   
        
    def flhXYZ(self, f, l, h):
        """ Funkcja zajmuje się transformacją współrzędnych geodezyjnych (fi, lambda, h) 
        na współrzędne ortokartezjańskie (X,Y,Z)
        Parameters 
        ----------
        f, l, h : FLOAT
                    współrzędne geodezyjne
                f - [radiany] - szerokość geodezyjna
                l - [radiany] - długość geodezyjna
                h - [metry] - wysokość elipsoidalna,
                
        Returns
        -------
        X, Y, Z - współrzędne w układzie ortokartezjańskim, 
                    : FLOAT
        """
        X_list = []
        Y_list = []
        Z_list = []

        for f, l, h in zip(f, l, h):
            N = self.Np(f)
            X = (N + h) * np.cos(f) * np.cos(l)
            Y = (N + h) * np.cos(f) * np.sin(l)
            Z = (N * (1 - self.ee) + h) * np.sin(f)
            X_list.append(X)
            Y_list.append(Y)
            Z_list.append(Z)

        return X_list, Y_list, Z_list
        
    def macierzR(self,f,l):
        """ Funkcja tworzy macierz R potrzebną do transformacji współrzędnych ortokartezjańskich  (X, Y, Z) do 
        współrzędnych układu NEU.
        Parameters
        ----------
        f, l : FLOAT
        parametry to szerokość i długość geodezyjna [radiany]
        
        Returns
        -------
        R : FLOAT
                R jest macierzą dwuwymiarową,
        """
        
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                           [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                           [np.cos(f),          0,     np.sin(f)]])
        return(R)
        
    
    
    def XYZ2neu(self, X_list, Y_list, Z_list, x0, y0, z0):
        """ Funcja transformuje współrzędne ortokartezjańskie (X,Y,Z) do układu NEU (North, East,Up ), wykorzystując 
        do tego celu macierz R
        Parameters 
        ----------
        X, Y, Z : FLOAT
                    współrzędne ortokartezjańskie punktów
        x0, y0, z0 : FLOAT
                        współrzędne ostokartezjańskie punktu odniesienia
        
        Returns
        -------
        dx - wektor współrzędnych przekształconych do układu NEU
        """
    
        wynikiNEU= []
        for X, Y, Z in zip(X_list, Y_list, Z_list):
         dX = np.array([X - x0, Y - y0, Z - z0])
         R = macierzR(f, l) 
         dx = np.dot(np.transpose(R), dX)
         wynikiNEU.append(dx.tolist())
        return wynikiNEU

    
    def sigma(self, f):
        A0 = 1 - self.ee/4 - 3*self.ee**2/64 - 5*self.ee**3/256
        A2 = (3/8) * (self.ee + self.ee**2/4 + 15*self.ee**3/128)
        A4 = (15/256) * (self.ee**2 + 3*self.ee**3/4)
        A6 = (35*self.ee**3)/3072
        si = a * (A0 * f - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
        return(sigma)
    
    
    def fl22000(self, f, l):
            m0 = 0.999923
            wsp2000 = []
            
            for f, l in zip(f,l):
                l0 = 0 
                strefa = 0
                if l >= 13.5 and l < 16.5:
                    l0 = 15 * np.pi / 180
                    nr = 5
                elif l >= 16.5 and l < 19.5:
                    l0 = 18 * np.pi / 180
                    nr = 6
                elif l >= 19.5 and l < 22.5:
                    l0 = 21 * np.pi / 180
                    nr = 7
                elif l >= 22.5 and l < 25.5:
                    l0 = 24 * np.pi / 180
                    nr = 8
                else:
                    l0 = None
                    nr = None
                    continue
                             
                e2prim = (self.a**2 - self.b**2) / self.b**2
                dl = l - l0
                t = np.tan(f)
                n = np.sqrt(e2prim * (np.cos(f))**2)
                N = self.Np(fi)
                Sigma = self.sigma(f)
                            
                XGK = Sigma + ((dl**2)/2) * N * np.sin(f)*np.cos(f) * ( 1+ ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(n**2) + 4*(n**4)     )  + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(n**2) - 330*(n**2)*(t**2))  )
                YGK = (dl*N* np.cos(f)) * (1+(((dl)**2/6)*(np.cos(f))**2) *(1-(t**2)+(n**2))+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(n**2)-58*(n**2)*(t**2)) )
                             
                X2000 = xgk * m0
                Y2000 = ygk * m0 + strefa*1000000 + 500000
                wsp.append([X2000, Y2000])
                    
            return(wsp2000) 

    def fl21992(self, f, l):
        l0 = (19 * np.pi)/180
        m0 = 0.9993
        wsp1992 = []
        for f,l in zip(f,l):
            e2prim = (self.a**2 - self.b**2) / self.b**2   
            dl = l - l0
            t = np.tan(f)
            n = np.sqrt(e2prim * (np.cos(f))**2)
            N = self.Np(f)
            Sigma = self.sigma(f)
                 
            XGK = Sigma + ((dl**2)/2)*N*np.sin(f)*np.cos(f) * ( 1+ ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(n**2) + 4*(n**4) ) + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(n**2) - 330*(n**2)*(t**2))  )
            YGK = (dl*N* np.cos(f)) * (1+(((dl)**2/6)*(np.cos(f))**2) *(1-(t**2)+(n**2))+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(n**2)-58*(n**2)*(t**2)) )
            
            X1992 = XGK * m0 - 5300000
            Y1992 = YGK * m0 + 500000
            wsp.append([X1992, Y1992]) 
            
        return(wsp1992)
  
 
   
    
  
        
        
if __name__ == "__main__":
    
    # utworzenie obiektu
    model = sys.argv
    try:
        model = sys.argv[3]
    except IndexError:
      print("model")
 
    geo = Transformacje(model = model)
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)        

def wczytaj_plik_wspolrzednych(nazwa_pliku):
    X = []
    Y = []
    Z = []
    try:
        with open(nazwa_pliku, 'r') as plik:
            linie = plik.readlines()
            for linia in linie:
                if not linia.startswith('#') and not linia.startswith('X'):  
                    kolumny = linia.strip().split(',')  
                    if len(kolumny) >= 3:
                        try:
                            x = float(kolumny[0])  
                            y = float(kolumny[1])  
                            z = float(kolumny[2])  
                            X.append(x)
                            Y.append(y)
                            Z.append(z)
                        except ValueError:
                            print("Nieprawidłowe dane w linii:", linia)
    except FileNotFoundError:
        print("Plik", nazwa_pliku, "nie został znaleziony.")
    return X, Y, Z


nazwa_pliku = 'wsp_inp.txt'  
X, Y, Z = wczytaj_plik_wspolrzednych(nazwa_pliku)




def zapisz_do_pliku_txt(wyniki, nazwa_pliku):
    try:
        with open(nazwa_pliku, 'w') as plik:
            for wiersz in wyniki:
                plik.write(str(wiersz) + '\n')
        print("Wyniki zostały zapisane do pliku", nazwa_pliku)
    except IOError:
        print("Błąd podczas zapisywania do pliku", nazwa_pliku)

wyniki = [f, l, h, X, Y, Z, dx, x2000, y2000, x92, y92]
nazwa_pliku = "WYNIKI.txt"
zapisz_do_pliku_txt(wyniki, nazwa_pliku)
    # if '--plh2xyz' in sys.argv:
    #     x, y,z = flhxyz(f, l, h)
    
    

model = {'WGS84': [6378137.000, 6356752.31424518], 'GRS80': [6378137.000, 6356752.31414036], 'KRASOWSKI': [6378245.000, 6356863.019]}
transformacje = {'XYZ2flh': 'XYZ2flh', 'flh2XYZ': 'flh2XYZ','XYZ2neu': 'XYZ2neu', 'fl22000': 'fl22000', 'fl21992': 'fl21992'}

try:
    while True:
        if len(sys.argv) <= 1:
            args_model = input('Podaj nazwe elipsoidy: ')
            args_dane = input('Wklej sciezke do pliku txt z danymi: ')
            args_transformacja = input('Podaj nazwę transformacji, którą chcesz wykonać: ')
        else:
            args_model = sys.argv[1]
            args_dane = sys.argv[2]
            args_transformacja = sys.argv[3]

        obiekt = Transformacje(model[args_model.upper()])
        wyniki = obiekt.wczytywanie(args_dane, Transformacje[args_transformacja.upper()])
        
        print('Plik z wynikami został utworzony.')
        
        wybor = input('Jezeli chcesz wykonac kolejna transformacje wpisz TAK, jeśli chcesz zakonczyc KONIEC: ').upper()
        if wybor != 'TAK':
            break

except FileNotFoundError:
    print('Podany plik nie istnieje.')
except KeyError:
    print('Zła podana elipsoida lub transformacja.')
except IndexError:
    print('Zły format danych w pliku.')
except ValueError:
    print('Zły format danych w pliku.')
finally:
    print('Koniec programu')
        
        
        
    
    
        
        
        
        
        