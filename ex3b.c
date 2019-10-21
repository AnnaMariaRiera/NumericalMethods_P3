#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
   int i,j,k;
   double n;
   double Pi;
   double resposta; // Definim una variable resposta en precisió doble que utilitzarem per saber com volem els nodes
   double* x; // Definim un vector x que serà el vector dels x_j
   double* fx; // Definim un vector fx que serà el vector de les solucions d'avaluar la funció 1/(3+x) als x_j
   double* l; // Definim un vector l que serà el vector dels l_j, els denominadors dels L_j a partir dels quals podem calcular el polinomi interpolador
   double* num; //Definim un vector dels numeradors dels L_j amb les x sent les x_k de l'enunciat.
   double* x1; // Definim un vector x1 que serà el vector dels x_k segons l'enunciat
   double* fx1; // Definim un vector fx1 que serà el vector de les solucions d'avaluar la funció (1/(3+x))+0.05sin(2Pinx) als x_k
   double* p;
   printf("Introdueix n\n");
   scanf("%lf",&n);
   Pi=3.14159265359879323;
   x=(double*)malloc((n+1)*sizeof(double));
   fx=(double*)malloc((n+1)*sizeof(double));
   l=(double*)malloc((n+1)*sizeof(double));
   num=(double*)malloc((n+1)*sizeof(double));
   x1=(double*)malloc(181*sizeof(double));
   fx1=(double*)malloc(181*sizeof(double));
   p=(double*)malloc(181*sizeof(double));
   
   printf("Vols nodes equiespaiats o de Txebyshev?\n 1-Equiespaiats\n 2-Txebyshev\n"); // Demanem a l'usuari quin tipus de nodes vol utilitzar
   scanf("%lf",&resposta);
   
   if(resposta==1)
   {
   
     for(j=0;j<=n;j++) // Calcularem els diferents x_j equiespaiats
     {
        x[j]=-1+j*(2/n);
		printf("x[%d]=%.16f\n",j,x[j]);
     }	
      
	  for(j=0;j<=n;j++) // Calcularem els diferents fx_j equiespaiats
     {
        fx[j]=(1/(3+x[j]))+0.05*sin(2*Pi*n*x[j]);
		printf("fx[%d]=%.16f\n",j,fx[j]);
     }	
	 
	 l[0]=(x[0]-x[1]); 
       for(i=2;i<=n;i++)  // Calculem el l_0, és a dir, el denominador de l'expressió de L_0
       {
          l[0]= l[0]*(x[0]-x[i]);		
       }
	   printf("l[0]=%.16f\n",l[0]);
	 for(j=1;j<=n;j++) // Calculem els diferents l_1 fins a l_n (l_j és el denominador de l'expressió de L_j)
	 { 
	   l[j]=(x[j]-x[0]);
       for(i=1;i<=n;i++)
       {
          if(i!=j)	
          {
             l[j]= l[j]*(x[j]-x[i]);
          }		
       }
	   printf("l[%d]=%.16f\n",j,l[j]);
	 }

	}
	
	
	
	
	if(resposta==2)
   {
   
     for(j=0;j<=n;j++) // Calcularem els diferents x_j de Txebyshev
     {
        x[j]=cos((((2*j)+1)*Pi)/((n+1)*2));
		printf("x[%d]=%.16f\n",j,x[j]);
     }	
	 
	 for(j=0;j<=n;j++) // Calcularem els diferents fx_j de Txebyshev
     {
        fx[j]=(1/(3+x[j]))+0.05*sin(2*Pi*n*x[j]);
		printf("fx[%d]=%.16f\n",j,fx[j]);
     }	
  
	 l[0]=(x[0]-x[1]); 
       for(i=2;i<=n;i++)  // Calculem el l_0, és a dir, el denominador de l'expressió de L_0
       {
          l[0]= l[0]*(x[0]-x[i]);		
       }
	   printf("l[0]=%.16f\n",l[0]);
	 for(j=1;j<=n;j++) // Calculem els diferents l_1 fins a l_n (l_j és el denominador de l'expressió de L_j)
	 { 
	   l[j]=(x[j]-x[0]);
       for(i=1;i<=n;i++)
       {
          if(i!=j)	
          {
             l[j]= l[j]*(x[j]-x[i]);
          }		
       }
	   printf("l[%d]=%.16f\n",j,l[j]);
	 }

	}
	
	
	
	printf("Amb això obtindrem el polinomi interpolador p(x), ara volem els punts (x,log(f(x)-p(x))) per determinats x[k] segons l'enunciat\n");
	
	for(k=0;k<=180;k++) // Calcularem els diferents x_k
     {
        x1[k]=-0.989+(0.011*k);
     }
	for(k=0;k<=180;k++) //Calculem els diferents f(x_k)
	{
	    fx1[k]=(1/(3+x1[k]))+(0.05*sin(2*Pi*n*x1[k]));
    }
	for(k=0;k<=180;k++) //Calculem els numeradors dels L_j en funció dels diferents x_k
	{  
	   num[0]=x1[k]-x[1];
	   for(i=2;i<=n;i++)  
       {
          num[0]= num[0]*(x1[k]-x[i]);		  
       }
	   
	   for(i=1;i<=n;i++)
	   {
	     num[i]=x1[k]-x[0];
		 for(j=1;j<=n;j++)
		 {
		   if(i!=j)
		   {
		     num[i]=num[i]*(x1[k]-x[j]);
		   }
		 }
	   }
	   p[k]=fx[0]*(num[0]/l[0]);  // A partir d'aquests numeradors, conjuntament amb els valors de la funció en els nodes interpoladors i els denominadors calculats abans, podem trobar el polinomi avaluat en cada x_k, per tant tindrem un vector de 181 posicions on cada component és el polinomi avaluat en un x_k diferent.
	   for(i=1;i<=n;i++)
	   {
	       p[k]=p[k]+(fx[i]*(num[i]/l[i]));
	   }
	}
	for(k=0;k<=180;k++)
	{
	  printf("[%lf,%lf],",x1[k],log(fabs(fx1[k]-p[k])));
	}
	
	
	
	
	
	
return 0;
}