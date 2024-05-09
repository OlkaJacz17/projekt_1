# projekt_1
Transformacja współrzędnych geodezyjnych

Utworzony przez nas program służy do transformacji współrzędnych
między układami. 

Obsługuje trzy elipsoidy : wgs84,GRS80 oraz elipsoidę Krasowskiego


____Transformacja współrzędnych ortokartezjańskich na geodezyjne (szerokość, długość i wysokość), to znaczy:

__XYZ ==> BLH__

Program przyjmuje współrzędne ortokartezjańskie i przy użyciu algorytmu Hirvonena przekształca je na współrzędne geodezyjne, gdzie:

B - szerokość geodezyjna, program zwraca tą wartość w radianach,

L - długość geodezyjna, również zwracana w radianach 

H - wysokość, odległość od elipsoidy tą wartość otrzymujemy w metrach.




____Transformacja odwrotna, przekształca współrzędne geodezyjne na ortokartezjańskie:

__BLH ==> XYZ__

Do programu wprowadzamy zmienne B, L podawane w radianach oraz H w metrach, w wyniku otrzymujemy współrzędne X, Y, Z w metrach.


____Transformacja ze współrzędnych ortokartezjańskich do topocentrycznych NEU (North, East, Up), w wyniku tej transformacji otrzymujemy tablicę z wartościami NEU, które są podane w metrach. 

__XYZ ==> NEU__

_____Transformacja współrzędnych geodezyjnych na współrzędne w układzie PL2000, które program zwróci nam w metrach.Program obsługuje tutaj wszytskie trzy elipsoidy.

__BL ==> PL2000__

_____Transformacja analogiczna do powyższej, wprowadzając współrzędne geodezyjne (szerokość oraz długość) program zwróci współrzędne w układzie PL1992 podane w metrach. Program obsługuje tutaj wszystkie trzy elipsoidy.

__BL ==> PL1992__

______Wymagania programu__
Do poprawnego działania programu należy skorzystać z pythona w wersji 3.11 a także zainstalowaną bibliotekę numpy oraz sys. Program został napisany dla systemu operacyjnego 
