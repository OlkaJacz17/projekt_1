# projekt_1
Transformacja współrzędnych geodezyjnych

Utworzony przez nas program służy do transformacji współrzędnych
między układami. 

Obsługuje trzy elipsoidy : wgs84,GRS80 oraz elipsoidę Krasowskiego


__Transformacja współrzędnych ortokartezjańskich na geodezyjne (szerokość, długość i wysokość), to znaczy:

XYZ ==> BLH

Program przyjmuje współrzędne ortokartezjańskie i przy użyciu algorytmu Hirvonena przekształca je na współrzędne geodezyjne, gdzie:

B - szerokość geodezyjna, program zwraca tą wartość w radianach,

L - długość geodezyjna, również zwracana w radianach 

H - wysokość, odległość od elipsoidy tą wartość otrzymujemy w metrach.




__Transformacja odwrotna, przekształca współrzędne geodezyjne na ortokartezjańskie:

BLH ==> XYZ

Do programu wprowadzamy zmienne B, L podawane w radianach oraz H w metrach, w wyniku otrzymujemy współrzędne X, Y, Z w metrach.


__Transformacja ze współrzędnych ortokartezjańskich do topocentrycznych NEU (North, East, Up), w wyniku tej transformacji otrzymujemy tablicę z wartościami NEU, które są podane w metrach. 
XYZ ==> NEU

__Trnsformacja współrzędnych geodezyjnych na współrzędne w układzie PL2000, które program zwróci nam w metrach.
BL ==> PL2000

__Transformacja analogiczna do powyższej, wprowadzając współrzędne geodezyjne (szerokość oraz długość) program zwróci współrzędne w układzie PL1992 podane w metrach.

BL ==> PL1992
