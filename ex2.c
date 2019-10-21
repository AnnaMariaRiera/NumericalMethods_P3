//Problema 2
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
  double x[12]={1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0}; 
  double J0[12]={0.281818559374385, 0.223890779141236, 0.166606980331990, 0.110362266922174, 0.055539784445602, 0.002507683297244, -0.048383776468198, -0.096804954397038, -0.142449370046012, -0.185036033364387, -0.224311545791968, -0.260051954901934};
  double coefsP1[2][2];
  double  coefsN1[2][2];
  double  coefsS1[2][2];
  double coefsP3[4][4];
  double  coefsN3[4][4];
  double  coefsS3[4][4];
  double coefsP5[6][6];
  double  coefsN5[6][6];
  double  coefsS5[6][6]; 
  //Matrius que guardaran els coeficients del mètode de les diferències dividides: interpolant valors positius, negatius i simètrics respectivament, per a polinomis de grau 1, 3 i 5
  double polinomi; //Variable per a avaluar el polinomi interpolador en 0
  double x1; //Variable auxiliar per fer el polinomi interpolador
  int n1=1, n3=3, n5=5; //Nombre de nodes
  int i, j;
  
  //Farem Newton invers: posem al revés les primeres dos columnes de la taula i apliquem l'algoritme és a dir, de f(x) trobem f-1(f(x))=x
  printf("\nApartat A: interpolacio de valors positius.\n");
  printf("Interpolacio inversa de grau 1:\n");  //Grau 1: 2 nodes
  //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n1; i++){
	   coefsP1[i][0]=J0[5-n1+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n1; j++){
	  for (i=0 ; i<=n1 ; i++){
			if (j>i) {
			  coefsP1[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsP1[i][j]= (coefsP1[i][j-1]-coefsP1[i-1][j-1])/(x[5-n1+i]-x[5-n1+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n1 ; i++)
	{ 
		for ( j=0 ; j<=n1 ; j++)
		{
			printf("%f ", coefsP1[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n1 ; i++){
	  printf("%.16f\n", coefsP1[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n1; i++){
		x1=x1*(-x[5-n1+i]);
		polinomi = coefsP1[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsP1[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
  
   printf("Interpolacio inversa de grau 3:\n");  //Grau 3: 4 nodes
   //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n3; i++){
	   coefsP3[i][0]=J0[5-n3+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n3; j++){
	  for (i=0 ; i<=n3 ; i++){
			if (j>i) {
			  coefsP3[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsP3[i][j]= (coefsP3[i][j-1]-coefsP3[i-1][j-1])/(x[5-n3+i]-x[5-n3+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n3 ; i++)
	{ 
		for ( j=0 ; j<=n3 ; j++)
		{
			printf("%f ", coefsP3[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n3 ; i++){
	  printf("%.16f\n", coefsP3[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n3; i++){
		x1=x1*(-x[5-n3+i]);
		polinomi = coefsP3[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsP3[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
   
   
   
   printf("Interpolacio inversa de grau 5:\n");  //Grau 5: 6 nodes 
   //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n5; i++){
	   coefsP5[i][0]=J0[5-n5+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n5; j++){
	  for (i=0 ; i<=n5; i++){
			if (j>i) {
			  coefsP5[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsP5[i][j]= (coefsP5[i][j-1]-coefsP5[i-1][j-1])/(x[5-n5+i]-x[5-n5+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n5 ; i++)
	{ 
		for ( j=0 ; j<=n5 ; j++)
		{
			printf("%f ", coefsP5[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n5 ; i++){
	  printf("%.16f\n", coefsP5[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n5; i++){
		x1=x1*(-x[5-n5+i]);
		polinomi = coefsP5[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsP5[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  

   
  
   //...................................................................................................................................................... 
   printf("\nApartat B: interpolacio de valors negatius.\n");
   printf("Interpolacio inversa de grau 1:\n");  //Grau 1: 2 nodes
  //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n1; i++){
	   coefsN1[i][0]=J0[6+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n1; j++){
	  for (i=0 ; i<=n1 ; i++){
			if (j>i) {
			  coefsN1[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsN1[i][j]= (coefsN1[i][j-1]-coefsN1[i-1][j-1])/(x[6+i]-x[6+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n1 ; i++)
	{ 
		for ( j=0 ; j<=n1 ; j++)
		{
			printf("%f ", coefsN1[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n1 ; i++){
	  printf("%.16f\n", coefsN1[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n1; i++){
		x1=x1*(-x[6+i]);
		polinomi = coefsN1[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsN1[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
  
  
  
   printf("Interpolacio inversa de grau 3:\n");  //Grau 3: 4 nodes
   //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n3; i++){
	   coefsN3[i][0]=J0[6+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n3; j++){
	  for (i=0 ; i<=n3 ; i++){
			if (j>i) {
			  coefsN3[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsN3[i][j]= (coefsN3[i][j-1]-coefsN3[i-1][j-1])/(x[6+i]-x[6+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n3 ; i++)
	{ 
		for ( j=0 ; j<=n3 ; j++)
		{
			printf("%f ", coefsN3[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n3 ; i++){
	  printf("%.16f\n", coefsN3[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n3; i++){
		x1=x1*(-x[6+i]);
		polinomi = coefsN3[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsN3[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
   
   
   
   printf("Interpolacio inversa de grau 5:\n");  //Grau 5: 6 nodes 
   //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n5; i++){
	   coefsN5[i][0]=J0[6+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n5; j++){
	  for (i=0 ; i<=n5; i++){
			if (j>i) {
			  coefsN5[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsN5[i][j]= (coefsN5[i][j-1]-coefsN5[i-1][j-1])/(x[6+i]-x[6+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n5 ; i++)
	{ 
		for ( j=0 ; j<=n5 ; j++)
		{
			printf("%f ", coefsN5[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n5 ; i++){
	  printf("%.16f\n", coefsN5[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n5; i++){
		x1=x1*(-x[6+i]);
		polinomi = coefsN5[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsN5[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
  
  
  
   //......................................................................................................................................................  
   printf("\nApartat C: interpolacio de valors simetrics.\n");
   printf("Interpolacio inversa de grau 1:\n");  //Grau 1: 2 nodes
  //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n1; i++){
	   coefsS1[i][0]=J0[5+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n1; j++){
	  for (i=0 ; i<=n1 ; i++){
			if (j>i) {
			  coefsS1[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsS1[i][j]= (coefsS1[i][j-1]-coefsS1[i-1][j-1])/(x[5+i]-x[5+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n1 ; i++)
	{ 
		for ( j=0 ; j<=n1 ; j++)
		{
			printf("%f ", coefsS1[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n1 ; i++){
	  printf("%.16f\n", coefsS1[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n1; i++){
		x1=x1*(-x[5+i]);
		polinomi = coefsS1[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsS1[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
  
  
  
   printf("Interpolacio inversa de grau 3:\n");  //Grau 3: 4 nodes
   //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n3; i++){
	   coefsS3[i][0]=J0[4+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n3; j++){
	  for (i=0 ; i<=n3 ; i++){
			if (j>i) {
			  coefsS3[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsS3[i][j]= (coefsS3[i][j-1]-coefsS3[i-1][j-1])/(x[4+i]-x[4+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n3 ; i++)
	{ 
		for ( j=0 ; j<=n3 ; j++)
		{
			printf("%f ", coefsS3[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n3 ; i++){
	  printf("%.16f\n", coefsS3[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n3; i++){
		x1=x1*(-x[4+i]);
		polinomi = coefsS3[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsS3[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
   
   
   
   printf("Interpolacio inversa de grau 5:\n");  //Grau 5: 6 nodes 
   //Rellenem totes les entrades de les matrius
  //Primera columna de la matriu: J0
   for (i=0 ; i<=n5; i++){
	   coefsS5[i][0]=J0[3+i];
   }
  //La resta de columnes
  for (j=1 ; j<=n5; j++){
	  for (i=0 ; i<=n5; i++){
			if (j>i) {
			  coefsS5[i][j] = 0; //Per sobre la diagonal tot són 0
			}
			else{
				coefsS5[i][j]= (coefsS5[i][j-1]-coefsS5[i-1][j-1])/(x[3+i]-x[3+i-j]);  //La diagonal i per baix guardem els valors de la taula de Newton
			}
		}
	}
   //Imprimim les matrius per pantalla per revisar que tot estigui bé, imprimim els coeficients que toquen per separat i avaluem el polinomi interpolador en 0
   printf("\nMatriu de les diferencies dividides:\n");
   printf("\n");
	for (i=0 ; i<=n5 ; i++)
	{ 
		for ( j=0 ; j<=n5 ; j++)
		{
			printf("%f ", coefsS5[i][j]);
		}
		printf("\n");
	}
	printf("\n");
    printf("Els coeficients del polinomi interpolador resultant son:\n");
    for (i=0 ; i<=n5 ; i++){
	  printf("%.16f\n", coefsS5[i][i]);
    }
   printf("\n");
   polinomi=0;
   x1=1;
   for (i=0; i<n5; i++){
		x1=x1*(-x[3+i]);
		polinomi = coefsS5[i+1][i+1]*x1 + polinomi;	
	}
   polinomi = coefsS5[0][0] + polinomi; //Tenim el polinomi avaluat a 0
   printf("El polinomi interpolador avaluat en 0 (i per tant, el valor x* que volem) es: %.16f\n\n", polinomi);
  
  return 0;
}