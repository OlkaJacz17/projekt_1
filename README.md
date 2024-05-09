# projekt_1
Transformacja współrzędnych geodezyjnych

Utworzony przez nas program służy do transformacji współrzędnych
między układami. 

Obsługuje trzy elipsoidy : wgs84,GRS80 oraz elipsoidę Krasowskiego


Transformacja współrzędnych ortokartezjańskich na geodezyjne (szerokość, długość i wysokość), to znaczy:

XYZ ==> BLH

Program przyjmuje współrzędne ortokartezjańskie i przy użyciu algorytmu Hirvonena przekształca je na współrzędne geodezyjne, gdzie:

B - szerokość geodezyjna, program zwraca tą wartość w radianach,

L - długość geodezyjna, również zwracana w radianach 

H - wysokość, odległość od elipsoidy tą wartość otrzymujemy w metrach.


Transformacja odwrotna, przekształca współrzędne geodezyjne na ortokartezjańskie:

BLH ==> XYZ

Do programu wprowadzamy :


XYZ ==> NEU



