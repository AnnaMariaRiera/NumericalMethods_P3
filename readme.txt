Exercici 1:
El programa en C a utilitzar �s el que s'anomena ex1.c 
En primer lloc, per compilar el programa, utilitzem la instrucci� seg�ent:
gcc -Wall -o ex1 ex1.c -lm
A continuaci�, escrivim: 
./ex1
Despr�s ens demana el nombre n de l'enunciat; hem de posar 4 i tornar a executar el programa. Seguidament posem 8 i el tornem a executar. Fem aix� per 4, 8, 16, 32 i 64.
En cada cas, el programa ens imprimir� per pantalla:
- Els valors dels nodes equidistants i els valors dels nodes de Chebyshev.
- La matriu que representa la taula del m�tode de les difer�ncies dividides generada tant amb nodes equidistants com de Chebyshev.
- Per nodes equispaiats:
	- Els coeficients resultants del polinomi interpolador (que s�n els valors de la diagonal de la matriu anterior)
	- El valor absolut de les difer�ncies entre la funci� i el polinomi, ambd�s avaluats als punts xk de l'enunciat de l'apartat (a). Aquests �ltims valors sortiran en forma de coordenades per a facilitar fer-ne el gr�fic.
	- El valor m�xim dels valors absoluts anteriors.
- El mateix per nodes de Chebyshev.
----------------------------------------------------------------------------------------------------------------------------------------
Exercici 2: 
El programa en C a utilitzar �s el que s'anomena ex2.c 
En primer lloc, per compilar el programa, utilitzem la instrucci� seg�ent:
gcc -Wall -o ex2 ex2.c -lm
A continuaci�, escrivim: 
./ex2
El programa ens imprimeix:
    - Matriu de les difer�ncies dividides amb 6 decimals
    - Coeficients del polinomi interpolador en 16 decimals
    - Valor del polinomi interpolador avaluat en 0
Per cada apartat (interpolant valors positius, negatius i sim�trics) i per cadascun dels graus del polinomi (1, 3 i 5) en cada cas.

----------------------------------------------------------------------------------------------------------------------------------------
Exercici 3:

Apartat a:
El programa en C a utilitzar en aquest apartat �s el que s'anomena ex3.c
En primer lloc, compilem el programa amb la seg�ent instrucci�:
gcc -Wall -o ex3 ex3.c -lm
A continuaci� l'executem:
./ex3
El programa ens demana primerament que introdu�m la n, que al nostre cas seran n=2,...,16, per tant, anirem executant el programa per cadascuna de les n desitjades. Per exemple si volem n=2 escriurem: 2
A continuaci� ens demana quins nodes volem. Primer demanarem nodes equidistants, per tant: 1
Quan tornem a executar el programa altra volta seguint els mateixos passos, i amb la mateixa n, li demanarem nodes de Txebyshev, per tant: 2.
�s a dir, per una mateixa n haurem d'executar el programa dues vegades, una per obtenir les dades per nodes equidistants i altra per obtenir les dades per nodes de Txebyshev.
El programa ens imprimeix:
     -Els nodes (els x_j)
     -El valor de la funci� f(x) avaluada als nodes
     -El valor dels coeficients l_j
     -Els punts de la forma (x,log(|f(x)-p(x)|)) pels valors de x_k de l'enunciat


Apartat b:
El programa en C a utilitzar en aquest apartat �s el que s'anomena ex3b.c
En primer lloc, compilem el programa amb la seg�ent instrucci�:
gcc -Wall -o ex3b ex3b.c -lm
A continuaci� l'executem:
./ex3b
El programa ens demana primerament que introdu�m la n, que al nostre cas seran n majors que 18, per tant, anirem executant el programa per cadascuna de les n desitjades. Per exemple si volem n=18 escriurem: 18
A continuaci� ens demana quins nodes volem. Primer demanarem nodes equidistants, per tant: 1
Quan tornem a executar el programa altra volta seguint els mateixos passos, i amb la mateixa n, li demanarem nodes de Txebyshev, per tant: 2.
�s a dir, per una mateixa n haurem d'executar el programa dues vegades, una per obtenir les dades per nodes equidistants i altra per obtenir les dades per nodes de Txebyshev.
El programa ens imprimeix:
     -Els nodes (els x_j)
     -El valor de la funci� f tilde avaluada als nodes
     -El valor dels coeficients l_j
     -Els punts de la forma (x,log(|f(x)-p(x)|)) pels valors de x_k de l'enunciat


----------------------------------------------------------------------------------------------------------------------------------------
Exercici optatiu:

El programa en C a utilitzar en aquest apartat �s el que s'anomena ex4.c
En primer lloc, compilem el programa amb la instrucci� seg�ent:
gcc -Wall -o ex4 ex4.c -lm
A continuaci� l'executem:
./ex4
El programa ens demanar� que introdu�m la n, i en el nostre cas la primera vegada que l'executem introduirem 4. Despr�s anirem executant-lo per n=8,16,32,64.
El programa ens imprimir� la matriu que s'ha de resoldre per trobar els moments.
Posteriorment ens demanar� si volem trobar el resultat utilitzant pivotatge. Com que la resposta es afirmativa, introdu�m: 1
El programa ens imprimir�:
      -Els coeficients de cada spline
      -Els punts de la corba de forma (x,|f(x)-s(x)|) per cada x_k