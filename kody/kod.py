from math import *
import numpy as np

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
        
        
        def xyz2bhl(self, X, Y, Z,):
            """" Funkcja zmienia współrządne ortokartezjańskie na współrzędne geodezyjne"""
            p = np.sqrt(X**2 + Y**2)
            f = np.arctan(Z/(p * (1 -self.ee)))
            while True: 
                N = Np(self,f)
                h = (p/np.cos(f)) - N
                fs = f
                f = np.arctan(Z/(p *(1 - (self.ee * (N/(N + h )))))) 
                if np.abs(fs-f) < (0.000001/206265):
                    break
            l = np.arctan2(Y, X)
            return(f,l,h)
        
        def Np(self, f): 
            N = a/np.sqrt(1 - ee * np.sin(f)**2)
            return(N)
            
    
        
        
        
        
        