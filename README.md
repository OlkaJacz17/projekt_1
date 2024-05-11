# projekt_1
Transformacja współrzędnych geodezyjnych

Utworzony przez nas program służy do transformacji współrzędnych
między układami. 

Obsługuje trzy elipsoidy : __wgs84__,  __GRS80__ oraz __elipsoidę Krasowskiego__

# __TRANSFORMACJE__


__XYZ ==> BLH__

_Transformacja współrzędnych ortokartezjańskich na geodezyjne_ (szerokość, długość i wysokość), to znaczy:

Program przyjmuje współrzędne ortokartezjańskie i przy użyciu algorytmu Hirvonena przekształca je na współrzędne geodezyjne, gdzie:

B - szerokość geodezyjna, program zwraca tą wartość w stopniach dziesiętnych,

L - długość geodezyjna, również zwracana w stopniach dziesiętnych 

H - wysokość, odległość od elipsoidy tą wartość otrzymujemy w metrach.

_______________________________________________________________________


__BLH ==> XYZ__

_Transformacja odwrotna, przekształca współrzędne geodezyjne na ortokartezjańskie_:

Do programu wprowadzamy zmienne B, L podawane w stopniach dziesiętnych oraz H w metrach, w wyniku otrzymujemy współrzędne X, Y, Z w metrach.

___________________________________________________________


__XYZ ==> NEU__

_Transformacja ze współrzędnych ortokartezjańskich do topocentrycznych NEU_ (North, East, Up), 
w wyniku tej transformacji otrzymujemy tablicę z wartościami NEU, które są podane w metrach. 


_____________________________________________________________

__BL ==> PL2000__

_Transformacja współrzędnych geodezyjnych na współrzędne w układzie PL2000_, które program zwróci nam w metrach.Program obsługuje tutaj wszytskie trzy elipsoidy.


____________________________________________________________

__BL ==> PL1992__

_Transformacja_ analogiczna do powyższej, wprowadzając współrzędne geodezyjne (szerokość oraz długość) program zwróci współrzędne w układzie PL1992 podane w metrach. Program obsługuje tutaj wszystkie trzy elipsoidy.


 
# WYMAGANIA PROGRAMU

Do poprawnego działania programu należy skorzystać z pythona w wersji 3.6 lub nowszych a także zainstalowaną bibliotekę math, numpy oraz sys. Program został napisany dla systemu operacyjnego Windows, macOS, Linux.

# KORZYSTANIE Z PROGRAMU

_Dane wejściowe_ - dane wprowadzane do programu powinny mieć format pliku tekstowego, w którym podane będą współrzędne w odpowiednich jednostakch zależnie od typu transformacji, który chcemy wykonać


_Rezulat_ - program zwraca wartości przeprowadzonych transformacji w postaci pliku tekstowego do folderu, w którym znajduje się nasz plik wejściowy
________________________________________________________________

_Wywołanie programu_ - do poprawnego zadziałania programu potrzebujemy sklonować to repozytorium na swoje urządzenie oraz uruchomić CMD  w skolonwanym folderze.


1. Aby program zadziałał poprawnie musimy użyć komendy " __python kod.py__" ,
2. Program poprosi nas o podanie nazwy elipsoidy, zatem wybieramy jedną z trzech dostępnych < __wgs84_,  __GRS80_, __elipsoidę Krasowskiego_ > wpisując odpowiednio < __wgs84__, __grs80__, __krasowskiego__ >
3. Następnie wprowadzamy nazwę transformacji jaką chcemy wykonać (musimy to zrobić przy pomocy flagi), więc należy wpisać jedną z podanych :
   " --xyz2flh"
   
   "--flh2xyz"
   
   "--fl22000"
   
   "--fl21992"
   
   "--xyz2neu" po wpisaniu tej transformacji musimy podać wartości współrzędnych (x_0, y_0, z_0) punktu odniesienia
5. Ostatnim elementem naszej komendy jest podanie nazwy pliku z danymi wejściowymi
   

