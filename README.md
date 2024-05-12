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

_Dane wejściowe_ - dane wprowadzane do programu powinny mieć format pliku tekstowego, w którym podane będą współrzędne w odpowiednich jednostakch zależnie od typu transformacji, który chcemy wykonać. Dane w pliku powinny być oddzielone znakiem ",". 


_Rezulat_ - program zwraca wartości przeprowadzonych transformacji w postaci pliku tekstowego do folderu, w którym znajduje się nasz plik wejściowy
________________________________________________________________

__Wywołanie programu__

do poprawnego zadziałania programu potrzebujemy sklonować to repozytorium na swoje urządzenie oraz uruchomić CMD  w skolonwanym folderze.


1. Aby program zadziałał poprawnie musimy na samym początku komendy użyć zwrotu :  " __python kod.py__" ,
2. Program poprosi nas o podanie nazwy elipsoidy, zatem wybieramy jedną z trzech dostępnych < __wgs84_,  __GRS80_, __elipsoidę Krasowskiego_ > wpisując odpowiednio < __wgs84__, __grs80__, __krasowskiego__ >
3. Następnie wprowadzamy nazwę transformacji jaką chcemy wykonać (musimy to zrobić przy pomocy flagi), więc należy wpisać jedną z podanych :
   
   " --xyz2flh"
   
   "--flh2xyz"
   
   "--fl22000"
   
   "--fl21992"
   
   "--xyz2neu" po wpisaniu tej transformacji musimy podać wartości współrzędnych (x_0, y_0, z_0) punktu odniesienia
4. Ostatnim elementem naszej komendy jest podanie nazwy pliku z danymi wejściowymi

# PRZYKŁADOWE UŻYCIE 
Na początku będziemy potrzebować pliku wejściowego, istotne jest aby dane w pliku były oddzielone od siebie znakiem ",".

Nasz przykładowy plik nosi nazwę wsp_inp.txt, a oto zawarte w nim dane:

Współrzedne geocentryczny ECEF stacji pemanentnej GNSS
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu
  X[m]         Y[m]        Z[m]
# -----------------------------------------------------
3664940.500,1409153.590,5009571.170
3664940.510,1409153.580,5009571.167
3664940.520,1409153.570,5009571.167
3664940.530,1409153.560,5009571.168
3664940.520,1409153.590,5009571.170
3664940.514,1409153.584,5009571.166
3664940.525,1409153.575,5009571.166
3664940.533,1409153.564,5009571.169
3664940.515,1409153.590,5009571.170
3664940.514,1409153.584,5009571.169
3664940.515,1409153.595,5009571.169
3664940.513,1409153.584,5009571.171
   
Chcemy z podanych współrzędnych (tutaj ortokartezjańskich) przejść do współrzędnych geodezyjnych (f,l,h), skorzystamy z elipsoidy grs80; zatem wpiszemy odpowiednią komendę:

_python kod.py grs80 --xyz2flh wsp_inp.txt_

(komendę należy wpisywać w porządku podanym powyżej, w przeciwnym razie kod nie zadziała)

Otrzymamy taki plik z nastepującymi danymi:

f, l, h 
52.09727222,21.03153333,141.399
52.09727216,21.03153314,141.400
52.09727212,21.03153296,141.403
52.09727209,21.03153277,141.408
52.09727209,21.03153323,141.410
52.09727212,21.03153318,141.402
52.09727207,21.031533,141.406
52.09727206,21.03153281,141.411
52.09727212,21.03153325,141.407
52.09727214,21.03153318,141.405
52.0972721,21.03153332,141.408
52.09727215,21.03153318,141.406

# UWAGI 
Aby kod działał prawidłowo należy zadbać o:
> porządek wpisywania kodu
> podanie wszystkich potrzebnych informacji
> poprawność zapisu komendy
> dbałość o jednostki w kodzie wejściowym oraz o separator danych, kod przyjmuje __tylko__ separator określany znakiem __" ,"__ .

# BŁĘDY I ROZWIĄZYWANIE
> jeśli program będzie podawał informację "podaj elipsoide", wtedy należy sprawdzić czy wszystkie wymagane informacje są zawarte w komendzie
> jeśli po wpisaniu komendy nie znajdziemy w folderze pliku z rezultatami informacji, należy sprawdzić czy :
> > w komendzie nie występują błędy
> > dane wejściowe mają odpowiedni separator
> > porządek kodu jest poprawny 


