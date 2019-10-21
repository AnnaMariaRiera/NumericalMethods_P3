//Interpolació per Newton (diferències dividides) de f(x)=1+1/(1+25*x*x) per x entre -1 i 1.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//Definim el nombre pi amb 17 decimals
#define P 3.14159265359879323

//Fem una funció void per a imprimir la matriu per pantalla
void PintaMatriu ( double** matriu, int n)
{
	int i, j;
	printf("\n");
	for (i=0 ; i<=n ; i++)
	{ 
		for ( j=0 ; j<=n ; j++)
		{
			printf("%f ", matriu[i][j]); 
			//Imprimim la matriu amb pocs decimals perquè només volem revisar que estigui bé
		}
		printf("\n");
	}
	printf("\n");
}


int main()
{
  //Totes les variables necessàries al problema les fem en precisió doble
  double *x; //Nodes equidistants (vector per a guardar tots els valors)
  double *y; //Nodes de Chebishev (vector per a guardar tots els valors)
  double **coefsX; //Definim una matriu que ens guardarà tots els coeficients del mètode de Newton amb nodes equidistants
  double **coefsY; //Definim una matriu que ens guardarà tots els coeficients del mètode de Newton amb nodes de Chebishev
  int n; //Nombre de nodes
  int i, j; //Variables auxiliars per al for i per referir-nos a les components de la matriu i els vectors
  double polinomiX, polinomiY; //Són els polinomis interpoladors de x i y resp. avaluats en les xk de l'apartat a
  double f; //Per a avaluar la funció f(x) en xk
  double x1, y1; //Variables auxiliars per a fer i avaluar els polinomis interpoladors
  double xk; //Valors per a estudiar l'error (donat a l'enunciat de l'apartat a)
  int k;
  double maxX, maxY; //Error màxim
  
  printf("Introdueix el nombre n (els de l'enunciat del problema):\n"); //El nombre de nodes que prendrem seran n+1 (ja que és de 0 fins a n)
  scanf("%d",&n);
  printf("Tindrem %d nodes, ja que van de 0 fins a %d.\n", n+1, n);
  
  //Reservem memòria per les matrius i els vectors (com tenim n+1 nodes, reservem per n+1 posicions)
  coefsX=(double**)malloc((n+1)*sizeof(double*)); //Primera fila
  for(i=0;i<=n;i++) //Una fila per cada element de la primera columna
  {
	  coefsX[i]=(double*)malloc((n+1)*sizeof(double));
  }
  coefsY=(double**)malloc((n+1)*sizeof(double*)); //Primera fila
  for(i=0;i<=n;i++) //Una fila per cada element de la primera columna
  {
	  coefsY[i]=(double*)malloc((n+1)*sizeof(double));
  }
  x=(double*)malloc((n+1)*sizeof(double));
  y=(double*)malloc((n+1)*sizeof(double));
  
  
  //Rellenem els vectors dels nodes
  printf("Els valors de xj amb els dos metodes son:\n");
  printf("Nodes equispaiats:\n");
  for(i=0;i<=n;i++)
  {
	  x[i]=-1+i*(2/((double)n)); //Obliguem al denominador a ser de tipus double, ja que si no ens divideix per un enter i no ens fa el que volem
	  printf("x[%d]=%.16f\n",i,x[i]);
  }
  printf("\nNodes de Chebishev:\n");
  for(i=0;i<=n;i++)
  {
	  y[i]=cos((2*i+1)/((double)(n+1))*P/2);
	  printf("y[%d]=%.16f\n",i,y[i]);
  }

  
  //Rellenem totes les entrades de les dos matrius mitjançant el mètode de les diferències dividides de Newton
  //NODES EQUISPAIATS:
  //Primera columna de la matriu: f(xi)
   for (i=0 ; i<=n; i++){
	   coefsX[i][0]=1+1/(1+25*x[i]*x[i]);
   }
  //La resta de columnes
  for (j=1 ; j<=n; j++){
	  for (i=0 ; i<=n ; i++){
			if (j>i) {
			  coefsX[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsX[i][j]= (coefsX[i][j-1]-coefsX[i-1][j-1])/(x[i]-x[i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //NODES DE CHEBYSHEV
   //Primera columna de la matriu: f(yi)
   for (i=0 ; i<=n; i++){
	   coefsY[i][0]=1+1/(1+25*y[i]*y[i]);
   }
  //La resta de columnes
  for (j=1 ; j<=n; j++){
	  for (i=0 ; i<=n ; i++){
			if (j>i) {
			  coefsY[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsY[i][j]= (coefsY[i][j-1]-coefsY[i-1][j-1])/(y[i]-y[i-j]); //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
    }
  
  
  //Imprimim les matrius per pantalla per revisar que tot estigui bé
   printf("\nMatriu de les diferencies dividides amb nodes equispaiats:\n");
   PintaMatriu(coefsX,n);
   printf("Matriu de les diferencies dividides amb nodes de Chebyshev:\n");
   PintaMatriu(coefsY,n);
   printf("\n");
  
  
  //Imprimim els coeficients del polinomi interpolador i fem que el programa ens imprimeixi el valor absolut de les diferències entre f(xk) i p(xk) per cada valor de xk donat per l'enunciat (apartat a)
  printf("NODES EQUISPAIATS:\n");
  printf("Els coeficients del polinomi interpolador resultant son:\n");
  for (i=0 ; i<=n ; i++){
	  printf("%.16f\n", coefsX[i][i]);
  }
  printf("\n");
  printf("El valor absolut de la diferencia entre f(xk) i p(xk) per cada xk es (en forma de punts, on el primer nombre es xk i el segon el val abs de la diferencia):\n");
  
  xk=0;
  maxX=0; 
  for (k=0; k<=180; k++){
	    polinomiX=0;
		xk=-0.989+k*0.011;
		x1=1;
		for (i=0; i<n; i++){
			x1=x1*(xk-x[i]);
			polinomiX = coefsX[i+1][i+1]*x1 + polinomiX;	
		}
		polinomiX = coefsX[0][0] + polinomiX; //Tenim el polinomi avaluat a xk
		f=1/(1+25*xk*xk); //f(xk)
		printf("[%f,%f],", xk, fabs(f-polinomiX));
		if ( (fabs(f-polinomiX)) > maxX){
			maxX=fabs(f-polinomiX);
		}
	}
  printf("\n\n");
  printf("El val abs de la diferencia maxima es %.16f\n", maxX);
  printf("\n\n\n");
  
  
  printf("NODES DE CHEBYSHEV:\n");
  printf("Els coeficients del polinomi interpolador resultant son:\n");
  for (i=0 ; i<=n ; i++){
	  printf("%.16f\n", coefsY[i][i]);
  }
  printf("\n");
  printf("El valor absolut de la diferencia entre f(xk) i p(xk) per cada xk es (en forma de punts, on el primer nombre es xk i el segon el val abs de la diferencia):\n");
  xk=0;
  maxY=0;
  for (k=0; k<=180; k++){
	    polinomiY=0;
		y1=1;
		xk=-0.989+k*0.011;
		for (i=0; i<n; i++){
			y1=y1*(xk-y[i]);
			polinomiY = coefsY[i+1][i+1]*y1 + polinomiY;	
		}
		polinomiY = coefsY[0][0] + polinomiY; //Tenim el polinomi avaluat a xk
		f=1/(1+25*xk*xk); //f(xk)
		printf("[%f,%f],", xk, fabs(f-polinomiY));
		if ( (fabs(f-polinomiY)) > maxY){
			maxY=fabs(f-polinomiY);
		}
	}
  printf("\n\n");
  printf("El val abs de la diferencia maxima es %.16f\n", maxY);
  printf("\n");
  
  
  //Alliberem memòria de les matrius
  for(i=0;i<=n;i++)
  {
     free(coefsX[i]);
     coefsX[i]=NULL;
  }
  free(coefsX);
  coefsX=NULL;
  for(i=0;i<=n;i++)
  {
     free(coefsY[i]);
     coefsY[i]=NULL;
  }
  free(coefsY);
  coefsY=NULL;
  
  //Alliberem memòria dels dos vectors
  free(x);
  x=NULL;
  free(y);
  y=NULL;
  
  return 0;
}