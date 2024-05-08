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
        
    def Np(self, f): 
            N = self.a/np.sqrt(1 - self.ee * np.sin(f)**2)
            return(N)
        
    def flhXYZ(self, f, l, h):
            """ Funkcja zmienia współrzędne geodezyjne na ortokartezjańskie"""
            N = Np(self, f)
            X = (N + h) * np.cos(f) * np.cos(l)
            Y = (N + h) * np.cos(f) * np.sin(l)
            Z = (N * (1-self.ee) + h) * np.sin(f)
            return(X, Y, Z)
        
    def macierzR(f,l):
             R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                           [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                           [np.cos(f),          0,     np.sin(f)]])
             return(R)
        
    def XYZ2neu(self, x, y, z, x0, y0, z0):
            f, l, h = self.xyz2bhl
            dX = array[x-x0,
                       y-y0,
                       z-z0]
            R = self.macierzR(f, l)
            dx= R.T @ dX
            return(dx)
    
    
    
        
        
        
if __name__ == "__main__":
    # utworzenie obiektu
    print(sys.argv)
    model = sys.argv[3]
    print(model)
    geo = Transformacje(model = model)
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)        
 
            
    if '--plh2xyz' in sys.argv:
        x, y,z = flhxyz(f, lh)
        
        
        
    
    
        
        
        
        
        