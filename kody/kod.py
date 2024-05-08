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
        N = self.a/np.sqrt(1 - self.ee * np.sin(f)**2)
        return(N)
        
    def xyz2plh(self, X, Y, Z):
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
            output [STR] - optional, defoulf 
                dec_degree - decimal degree
                dms - degree, minutes, sec"""
            p = np.sqrt(X**2 + Y**2)
            f = np.arctan(Z/(p * (1 -self.ee)))
            while True: 
                N = self.Np(f)
                h = (p/np.cos(f)) - N
                fs = f
                f = np.arctan(Z/(p *(1 - (self.ee * (N/(N + h )))))) 
                if np.abs(fs-f) < (0.000001/206265):
                    break
            l = np.arctan2(Y, X)
            return(f,l,h)
  
        
    def flhXYZ(self, f, l, h):
            """ Funkcja zmienia współrzędne geodezyjne na ortokartezjańskie"""
            N = self.Np(self, f)
            X = (N + h) * np.cos(f) * np.cos(l)
            Y = (N + h) * np.cos(f) * np.sin(l)
            Z = (N * (1-self.ee) + h) * np.sin(f)
            return(X, Y, Z)
        
    def macierzR(self,f,l):
             R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                           [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                           [np.cos(f),          0,     np.sin(f)]])
             return(R)
        
    def XYZ2neu(self, X, Y, Z, x0, y0, z0):
            f, l, h = self.xyz2bhl
            dX = array[X-x0,
                       Y-y0,
                       Z-z0]
            R = self.macierzR(f, l)
            dx= R.T @ dX
            return(dx)
    
    def ustal_parametry(self, l):
        if l >= 13.5 and l < 16.5:
            l02000 = 15 * np.pi / 180
            nr = 5
        elif l >= 16.5 and l < 19.5:
            l02000 = 18 * np.pi / 180
            nr = 6
        elif l >= 19.5 and l < 22.5:
            l02000 = 21 * np.pi / 180
            nr = 7
        elif l >= 22.5 and l < 25.5:
            l02000 = 24 * np.pi / 180
            nr = 8
        else:
            l02000 = None
            nr = None
        return l02000, nr

    
    def sigma2000(self, f):
        A0 = 1 - self.ee/4 - 3*self.ee**2/64 - 5*self.ee**3/256
        A2 = (3/8) * (self.ee + self.ee**2/4 + 15*self.ee**3/128)
        A4 = (15/256) * (self.ee**2 + 3*self.ee**3/4)
        A6 = (35*self.ee**3)/3072
        si = a * (A0 * f - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
        return(si2000)

    def wsploknagausskruger(self,si,l, l02000 ,f):
        b2 = self.a**2 * (1-self.ee)
        e22 = (self.a**2 - b2) / b2
        dl = l - l02000
        t = tan(f)
        eta2 = e22 * (cos(f))**2
        N = a / (sqrt(1-e2 * sin(f)**2))
        xgk = si + (dl**2/2) * N * np.sin(f) * np.cos(f) * ((1 + (dl**2/12)*(np.cos(f))**2 * (5 -t**2 +9*eta2 + 4*eta2**2) + (dl**4/360) * np.cos(f)**4 * (61 - 58 * t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2)))
        ygk = dl * N * np.cos(f) * (1 + (dl**2/6) * np.cos(f)**2 * (1 - t**2 + eta2) + (dl**4/120) * cos(f)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*eta2*t**2))
        return(xgk2000, ygk2000)
  
    def uklad2000(xgk2000,ygk2000,nr):
       m2000 = 0.999923
       x2000 = xgk * m2000
       y2000 = ygk * m2000 + nr * 1000000 + 500000
       return(x2000,y2000)
   
 
   
    
    def wsploknagausskruger(self,si,l,f):
        l01992 = 19*np.pi/180
        b2 = self.a**2 * (1-self.ee)
        e22 = (self.a**2 - b2) / b2
        dl = l - l092
        t = tan(f)
        eta2 = e22 * (cos(f))**2
        N = self.a / (sqrt(1-e2 * sin(f)**2))
        xgk = si + (dl**2/2) * N * np.sin(f) * np.cos(f) * ((1 + (dl**2/12)*(np.cos(f))**2 * (5 -t**2 +9*eta2 + 4*eta2**2) + (dl**4/360) * np.cos(f)**4 * (61 - 58 * t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2)))
        ygk = dl * N * np.cos(f) * (1 + (dl**2/6) * np.cos(f)**2 * (1 - t**2 + eta2) + (dl**4/120) * cos(f)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*eta2*t**2))
        return(xgk, ygk)
    
    def uklad1992(xgk,ygk):
        m1992 = 0.9993
        x1992 = xgk * m1992 - 5300000
        y1992 = ygk * m1992 + 500000
        return(x1992,y1992)
        
        
        
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
        
        
        
    
    
        
        
        
        
        