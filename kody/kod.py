from math import *
import numpy as np
from argparse import *
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
        
  
    def xyz2flh(self,X, Y, Z): 
        """" Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokość elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
              współrzędne w układzie orto-kartezjańskim podawane w metrach

        Returns
        -------
        f : FLOAT
            [radiany] - szerokość geodezyjna
        l : FLOAT
            [radiany] - długośc geodezyjna.
        h : FLOAT
            [metry] - wysokość elipsoidalna
        """
        l = np.arctan2(Y, X)
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z / (p * (1 - self.ee)))
        while True :
            N = self.a/np.sqrt(1 - self.ee * np.sin(f)**2)
            h = p / np.cos(f) - N
            fs = f
            f = np.arctan(Z / (p * (1 - (self.ee * (N / (N + h))))))
            if np.abs(f - fs) < (0.000001/206265):
                break
        return(f, l, h)   
  

    
    def flh2xyz(self, f, l, h):
        """ Funkcja zajmuje się transformacją współrzędnych geodezyjnych (fi, lambda, h) 
        na współrzędne ortokartezjańskie (X,Y,Z)
        Parameters 
        ----------
        f, l, h : FLOAT
                    współrzędne geodezyjne:
                f - [radiany] - szerokość geodezyjna
                l - [radiany] - długość geodezyjna
                h - [metry] - wysokość elipsoidalna,
        
        Returns
        -------
        X, Y, Z : FLOAT
                    - współrzędne w układzie ortokartezjańskim, których wynik zwracany jest w metrach
        """
        N = self.a/np.sqrt(1 - self.ee * np.sin(f)**2)
        X = (N + h) * np.cos(f) * np.cos(l)
        Y = (N + h) * np.cos(f) * np.sin(l)
        Z = (N * (1-self.ee) + h) * np.sin(f)
        return(X, Y, Z)
   
    
    def xyz2neu(self, X, Y, Z, x_0, y_0, z_0):
        """ 
        Funkcja przeprowadza transformację ze współrzędnych ortokartezjańskich (X, Y, Z) 
        do topocentrycznych NEU (North, East, Up) z wykorzystaniem potrzebnej do tego macierzy R.
        Parameters
        ----------
        X, Y, Z : FLOAT
                    współrzędne ortokartezjańskie, te wartosci wprowadzamy w metrach
        x_0, y_0, z_0 : FLOAT
                    współrzędne ortokartezjańskie punktu odniesienia
        
        Returns
        -------
        n : FLOAT
                wartość współrzędnej North w metrach
        e : FLOAT
                wartość współrzędnej East w metrach
        up : FLOAT
                wartość współrzędnej Up w metrach

        """
        
        f, l, _ = [radians(coord) for coord in self.xyz2plh(X, Y, Z)]
        
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                           [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                           [np.cos(f),          0,     np.sin(f)]])
        
        xyz_t = np.array([[X-x_0],
                       [Y-y_0],
                       [Z-z_0]])
        dx = R.T @ xyz_t
        n = dx[0]
        e = dx[1]
        up = dx[2]
        return n,e,up

    
    def sigma(self, f):
        """ 
        Funkcja okresla wartoć sigma dzięki podanej wartosci f, funkcja jest kluczowa w transformacji współrzędnych geodezyjnych (fi, lambda, h)
        do układu PL2000, którego współrzędnycmi są (X2000, Y2000) podane w metrach.
        Parameters 
        ----------
        f : FLOAT
                szerokosć geodezyjna podawana w radianach
                
        Returns 
        -------
        sigma : FLOAT
                wartosć używana w transformacji ze współrzędnych geocentrycznych do układu PL2000, wartosc zwracana w metrach
        """
        A0 = 1 - self.ee/4 - 3*self.ee**2/64 - 5*self.ee**3/256
        A2 = (3/8) * (self.ee + self.ee**2/4 + 15*self.ee**3/128)
        A4 = (15/256) * (self.ee**2 + 3*self.ee**3/4)
        A6 = (35*self.ee**3)/3072
        sigma = self.a * (A0 * f - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
        return(sigma)
    
    
    
    
    def fl22000(self, f, l):
     """  
     Funkcja ma przekształca współrzędne geodezyjne (fi, lambda) na współrzędne w układzie PL2000 (X2000, Y2000), funkcja wykorzystuje
     wartosć sigma do transformacji współrzędnych
     Parameters
     ----------
     f : FLOAT
             szerokosć geodezyjna, wartosć podawana w radianach
     l : FLOAT
             dlugosc geodezyjna, wartosc podawana w radianach
             
     Returns
     -------
     X2000 : FLOAT
             współrzędna otrzymana w układzie PL2000, jednostką zwróconej wartosci są metry
     Y2000 : FLOAT
             współrzędna otrzymana w układzie PL2000, jednostką zwróconej wartosci są metry
     
     """
     m0 = 0.999923
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
         
     e2prim = (self.a**2 - self.b**2) / self.b**2
     dl = l - l0
     t = np.tan(f)
     n = np.sqrt(e2prim * (np.cos(f))**2)
     N = self.Np(f)
     Sigma = self.sigma(f)
    
     XGK = Sigma + ((dl**2)/2) * N * np.sin(f)*np.cos(f) * (1 + ((dl**2)/12)*(np.cos(f))**2 * (5 - (t**2) + 9*(n**2) + 4*(n**4)) + ((dl**4)/360)*(np.cos(f)**4) * (61 - 58*(t**2) + (t**4) + 270*(n**2) - 330*(n**2)*(t**2)))
     YGK = (dl*N* np.cos(f)) * (1 + (((dl)**2/6)*(np.cos(f))**2) * (1 - (t**2) + (n**2)) + ((dl**4)/120)*(np.cos(f)**4) * (5 - 18*(t**2) + (t**4) + 14*(n**2) - 58*(n**2)*(t**2)))
                             
     X2000 = XGK * m0
     Y2000 = YGK * m0 + strefa*1000000 + 500000
     return(X2000, Y2000)

 

    def fl21992(self, f, l):
        """ 
        Funcja przy pomocy wartosci sigma przekształca współrzędne geodezyjne (fi, lambda) do układu PL1992 (X1992, Y1992)
        Parameters
        ----------
        f : FLOAT
                szerokoć geodezyjna wyrażona w radianach
        l : FLOAT 
                długosć geodezyjna wyrażona w radianach
                
        Returns
        -------
        X1992 : FLOAT
                współrzędna w układzie PL1992, wyrażona w metrach
        Y1992 : FLOAT
                współrzędna w układzie PL1992, wyrażona w metrach
        
        """
        l0 = (19 * np.pi)/180
        m0 = 0.9993
        
        
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
            
            
        return(X1992, Y1992)
       
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    print(sys.argv)
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    f, l, h = geo.xyz2flh(X, Y, Z)
    print(f, l, h)
    X_new, Y_new, Z_new = geo.flh2xyz(f, l, h)
    print(X_new, Y_new, Z_new)

input_file_path = sys.argv[-1]
    
if '--xyz2flh' in sys.argv and '--flh2xyz' in sys.argv:
    print('mozezz podac tylko jedna flage')
elif'--xyz2flh' in sys.argv:
    with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[4:]
            # coords_neup = []   
            coords_flh = []
            for line in lines:
                line = line.strip()
                x_str, y_str, z_str = line.split(',')
                x, y,z = (float(x_str), float(y_str), float(z_str))
                f, l, h = geo.xyz2flh(x,y,z)
                h = '{:.3f}'.format(h)
                coords_flh.append([f,l,h])

                 
    with open('results_xyz2flh.txt', 'w') as f:
            f.write('f, l, h \n')
            for coords in coords_flh:    
                coords_flh_line = ','.join([str(coord) for coord in coords])
                f.write(coords_flh_line + '\n')

if '--flh2xyz' in sys.argv and '--xyz2flh' in sys.argv:
    print('mozezz podac tylko jedna flage')
elif'--flh2xyz' in sys.argv:
    with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            # coords_neup = []   
            coords_xyz = []
            for line in lines:
                line = line.strip()
                f_str, l_str, h_str = line.split(',')
                f,l,h = (float(f_str), float(l_str))
                x,y,z = geo.flh2xyz(f, l, h)
                h = '{:.3f}'.format(h)
                coords_xyz.append([x,y,z])

                 
    with open('results_flh2xyz.txt', 'w') as f:
            f.write('x, y, z \n')
            for coords in coords_flh:    
                coords_flh_line = ','.join([str(coord) for coord in coords])
                f.write(coords_xyz_line + '\n')


if '--fl22000' in sys.argv:
    print('mozezz podac tylko jedna flage')
elif'--flh22000' in sys.argv:
    with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            # coords_neup = []   
            coords_x2000y2000 = []
            for line in lines:
                line = line.strip()
                f_str, l_str = line.split(',')
                f,l = (float(f_str), float(l_str))
                x2000,y2000 = geo.flh22000(f, l)
                h = '{:.3f}'.format(h)
                coords_xyz.append([x2000,y2000])

                 
    with open('results_flh22000.txt', 'w') as f:
            f.write('x2000, y2000 \n')
            for coords in coords_flh:    
                coords_fl_line = ','.join([str(coord) for coord in coords])
                f.write(coords_x2000y2000_line + '\n')


if '--fl21992' in sys.argv:
    print('mozezz podac tylko jedna flage')
elif'--flh21992' in sys.argv:
    with open(input_file_path, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            # coords_neup = []   
            coords_x1992y1992 = []
            for line in lines:
                line = line.strip()
                f_str, l_str = line.split(',')
                f,l = (float(f_str), float(l_str))
                x1992,y1992 = geo.flh21992(f, l)
                h = '{:.3f}'.format(h)
                coords_xyz.append([x1992,y1992])

                 
    with open('results_flh21992.txt', 'w') as f:
            f.write('x1992, y1992 \n')
            for coords in coords_flh:    
                coords_fl_line = ','.join([str(coord) for coord in coords])
                f.write(coords_x1992y1992_line + '\n') 

elif'--xyz2neu' in sys.argv:
    with open(input_file_path, 'r') as f:
        lines = lines[4:]
        coords_neu = []
        for line in lines:
            line = line.strip()
            x_str, y_str, z_str = line.split(',')
            x, y,z = (float(x_str), float(y_str), float(z_str))
            x_0,y_0,z_0 = [float(coord) for coord in sys.argv[-4:-1]]
            n,e,up = geo.xyz2neup(x,y,z,x_0,y_0,z_0)
            coords_neup.append([n,e,up])
            
    with open('results_xyz2neu.txt', 'w') as f:
            f.write('n[m],e[m], u[m] \n')
            for coords in coords_neu:    
                coords_neu_line = ','.join([str(coord) for coord in coords])
                f.write(coords_neu_line + '\n') 

        

        
    
    
        
        
        
        
        