//Interpolació per splines cúbics naturals de f(x)=1/(1+25*x*x) per x entre -1 i 1.
#include<stdio.h>
#include<stdlib.h>   //Per a la memòria dinàmica
#include<math.h>

#define TOL 0.000001

/*Treballem amb 2 o 3 asteriscs (a les matrius) segons si les volem passar per paràmetre o per referència. Una matriu és un apuntador (*)
d'apuntadors (**), per això, de manera natural sempre portarà dos asteriscs. No obstant, a l'hora de cridar la matriu mitjançant una funció, 
hem de tenir sempre present què volem fer amb ella. Si només la copiem i la usem per a altres coses, i no en volem modificar res, la passarem 
com a paràmetre(els dos asteriscs normals). Però si sí que la volem modificar (com per exemple, a l'hora de col·locar-li tots
els seus valors [com a la funció LlegeixMatriuISolucions]) l'haurem de passar per referència: apuntant a la matriu, és a dir, ***.  */

void PintaMatriu ( double** matriu, int n)
{
	int i, j;
	printf("\n");
	for (i=0 ; i<n ; i++)
	{ 
		for ( j=0 ; j<n ; j++)
		{
			printf("%f ", matriu[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void Pivotar ( double** matriu, double** P, double** Q, int on_som, int n)
/*El paràmetre "on_som" és perquè no sempre voldrem que ens busqui l'element (en valor absolut) més gran de tota la matriu, sinó que a partir del segon,
voldrem que ens el busqui sense mirar les files anteriors. */
{
	int i, j, max_i, max_j;
	double intermig;

	max_i = on_som;
	max_j = on_som;
	
	for (i=on_som; i<n; i++) {
		for (j=on_som; j<n ; j++) {
			if ( abs( matriu[i][j] ) > abs( matriu[max_i][max_j] ))  {
				max_i=i;    /* matriu[i][j] no és el que hi ha DINS la caixeta de la posició ij, sinó que ÉS la caixeta */
				max_j=j;    /* Amb max_i=i, el que fem sí que es canviar el de DINS la caixeta */
			}
		}
	}
	/*La funció "Pivotar" és una funció que creem per a pivotar EL QUE VOLGUEM i QUAN volguem. Per tant, no li podem dir que ho faci des de la primera fila,
	sinó que volem que ens ho faci des d'on nosaltres li direm: és a dir, "on som". */
	 
	for (i=0; i<n; i++) {
		intermig = matriu[on_som][i];
		matriu[on_som][i] = matriu[max_i][i];
		matriu[max_i][i] = intermig;
		
		intermig = P[on_som][i];
		P[on_som][i] = P[max_i][i];
		P[max_i][i] = intermig;
	}
	/*Com al joc del solitari, necessitem una variable intermitja que ens ajudi a canviar de lloc les files. */
	
	for (j=0; j<n; j++) {
		intermig = matriu[j][on_som];
		matriu[j][on_som] = matriu[j][max_j];
		matriu[j][max_j] = intermig;
		
		intermig = Q[j][on_som];
		Q[j][on_som] = Q[j][max_j];
		Q[j][max_j] = intermig;
	}

}			

double** LU_SensePivotatge ( double** matriu, double*** lu, int n)
{
	int i, j, k;
	double multiplicador;
	double **LU, **L;
	
	//Reservem memòria per a les dues matrius que acabem de declarar
	LU = (double**) malloc (n*sizeof(double*)); /*Apuntador (per exemple) a cadascuna de les files*/
	L = (double**) malloc (n*sizeof(double*));
	
	for (i = 0; i < n; i++){   //Demanem memòria per cada columna de les matrius
		LU[i] = (double*) malloc (n*sizeof(double)); /*Apuntador (seguint el mateix exemple que abans) a cada columna: cada MEMBRE de cada fila*/
		L[i] = (double*) malloc (n*sizeof(double));
		for (j = 0; j < n; j++){
			LU[i][j] = matriu[i][j];
		}
	}
	

	for (i=0 ; i<n ; i++)
	{
		if ( fabs(LU[i][i])<TOL) {
			return NULL;
		}
		for (j=i+1 ; j<n ; j++)
		{
			multiplicador = LU[j][i] / LU[i][i];
			for (k=0 ; k<n ; k++)
			{
				LU[j][k]= LU[j][k] - (LU[i][k]*multiplicador);
			}
			L[j][i]=multiplicador;
		}
	}
	
	/*Imprimirem L i U per separat per a què l'usuari les pugui veure i emplenem la part inferior de LU*/
	//Emplenem LU
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<i ; j++){
			LU[i][j] = L[i][j];
		}
	}
	
	//Imprimim U
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = LU[i][j];
			if (j<i) {
				L[i][j] = 0;
			}
		}
	}
	
	//Imprimim L
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = 0;
			if (j==i){
				L[i][i] = 1 ;
			}
			else if (j<i) {
				L[i][j] = LU[i][j];
			}
		}
	}
	
	//Alliberem la memòria de L
	for (i=0; i<n; i++){
		free(L[i]);
		L[i] = NULL;
	}
	free(L);
	L = NULL;
	
	*lu = LU;
	
return (*lu);
}

double** CreaIdentitat (int n)
{
	int i, j;
	double** Id;
	
	Id = (double**) malloc (n*sizeof(double*));
	
	for (i=0; i < n; i++){
		Id[i] = (double*) malloc (n*sizeof(double));
		
		for (j = 0; j < n; j++){
			if ( i==j ) {
				Id[i][j] = 1;
			}
			else {
				Id[i][j]=0;
			}
		}
	}
	
return(Id);
}

double** LU_AmbPivotatge ( double** matriu, double*** lu, double** P, double** Q, int n)
{
	int i, j, k;
	double multiplicador;
	double **LU, **L;
	
	LU = (double**) malloc (n*sizeof(double*));
	L = (double**) malloc (n*sizeof(double*));
	
	for (i = 0; i < n; i++){   //inicialitzem les dos matrius (la LU i una auxiliar _L)
		LU[i] = (double*) malloc (n*sizeof(double));
		L[i] = (double*) malloc (n*sizeof(double));
		for (j = 0; j < n; j++){
			LU[i][j] = matriu[i][j];
		}
	}
	
	for (i=0 ; i<n ; i++)
	{
		Pivotar ( LU, P, Q, i, n );
		if ( fabs(LU[i][i])<TOL ) {
			printf ("\nEl determinant de la matriu es nul o algun pivot es menor que %f.\nEl metode no es pot executar.\n", TOL);
			return NULL;
		}
		for (j=i+1 ; j<n ; j++)
		{
			multiplicador = LU[j][i] / LU[i][i];
			for (k=0 ; k<n ; k++)
			{
				LU[j][k]= LU[j][k] - (LU[i][k]*multiplicador);
			}
			L[j][i]=multiplicador;
		}
	}
	
	//Emplenem LU
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<i ; j++){
			LU[i][j] = L[i][j];
		}
	}
	
	//Imprimim U
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = LU[i][j];
			if (j<i) {
				L[i][j] = 0;
			}
		}
	}
	
	//Imprimim L
	for (i=0 ; i<n ; i++){
		for (j = 0 ; j<n; j++){
				L[i][j] = 0;
			if (j==i){
				L[i][i] = 1 ;
			}
			else if (j<i) {
				L[i][j] = LU[i][j];
			}
		}
	}
	
	*lu = LU;
return (*lu);
}

double* TriangularInferior (double** LU, double* vector, int n)
{
	double *y;
	double sumatori; 
	int i, j;
	
	y = (double*) malloc (n*sizeof(double));
	
	y[0] = vector[0];
	for ( i=1; i<n; i++ ) {
		sumatori=0;
		for ( j=0; j<i; j++ ) {
			sumatori += LU[i][j]*y[j];
		}
		y[i] = vector[i] - sumatori;
	}
	
return y;
}

double* TriangularSuperior (double** LU, double* y, int n)
{
	double *x;
	double sumatori; 
	int i, j;
	
	x = (double*) malloc (n*sizeof(double));
	
	x[n-1] = y[n-1] / LU[n-1][n-1];  /*Perquè tot es conta de 0 fins a n-1*/
	for (i=n-2; i>=0; i--) {
		sumatori=0;
		for (j=n-1; j>=i; j--) {
			sumatori += LU[i][j]*x[j];
		}
		x[i] = ( y[i] - sumatori ) / LU[i][i]; 
	}
	
return x;
}

double* MatriuPerVector (double** m, double* v, int n)
{
	int i, j;
	double* VectFinal;
	
	VectFinal = (double*) malloc (n*sizeof(double));
	
	for (i=0; i<n; i++) {
		VectFinal[i] = 0;
		for (j=0; j<n; j++) {
			VectFinal[i] += m[i][j]*v[j];
		}
	}
	
return VectFinal;
}

double* MetodeComplet_SensePivotatge (double** matriu, double* vector, int n)
{
	int i;
	double** LU; 
	double *y, *resultat;
	
	LU = LU_SensePivotatge (matriu, &LU, n);
	if (LU == NULL){
		return NULL;
	}
	
	y = TriangularInferior (LU, vector, n);
	resultat = TriangularSuperior (LU, y, n);
	
	//Alliberem memòria:
	for (i=0; i<n; i++){
		free(LU[i]);
		LU[i] = NULL;
	}
	free(LU);
	LU = NULL;
	
	free(y);
	y = NULL;
	
return resultat;
}

double* MetodeComplet_AmbPivotatge (double** matriu, double* vector, int n)
{
	int i;
	double **LU, **P, **Q; 
	double *y, *resultat2, *vector2, *resultat;
	
	P = CreaIdentitat (n);
	Q = CreaIdentitat (n);
	LU = LU_AmbPivotatge( matriu, &LU, P, Q, n);  /* Li passem l'adreça de LU perquè la passem per referència: a la funció és *** */
	
	if (LU == NULL){
		return NULL;
	}
	
	vector2 = MatriuPerVector ( P, vector, n);
	
	y = TriangularInferior (LU, vector2, n);
	resultat2 = TriangularSuperior (LU, y, n);
	
	resultat = MatriuPerVector (Q, resultat2, n);
	
	//Alliberem memòria:
	for (i=0; i<n; i++){
		free(LU[i]);
		LU[i] = NULL;
	}
	free(LU);
	LU = NULL;
	
	for (i=0; i<n; i++){
		free(P[i]);
		P[i] = NULL;
	}
	free(P);
	P = NULL;
	
	for (i=0; i<n; i++){
		free(Q[i]);
		Q[i] = NULL;
	}
	free(Q);
	Q = NULL;
	
	free(y);
	y = NULL;
	
	free(resultat2);
	resultat2 = NULL;
	
	free(vector2);
	vector2 = NULL;
	
return resultat;
}

double* PivotatgeCondicional (double** matriu, double* vector, int n)
{
	int resposta;
	double *resultat;
	
	printf ("\nEl determinant de la matriu es nul o algun pivot es menor que %f.\nProva amb pivotatge:\n", TOL);
	
	printf("\nQue vols fer?\n 1. Vull provar amb pivotatge.\n 2. Vull sortir del programa.\n");
	scanf("%d", &resposta);
				
	if (resposta==1) {
		resultat = MetodeComplet_AmbPivotatge (matriu, vector, n);
		if (resultat == NULL) {
			return NULL;
		}
	}
	else if (resposta==2) {
		return NULL;
	}
	else {
		printf ("Has triat una opcio que no es valida\n");
		return 0;
	}
	
return resultat;
}

int main()
{
	int resposta, n, i, j; //Variables auxiliars i n, el nombre de l'enunciat del problema
	double **matriu;
	double *vector, *resultat;
	double *x; //xj nodes equispaiats
	double *y; //vector amb els valors f(xj)
	double *h; //vector amb hi=x_i-x_{i-1}
	double *d;
	double *lambda;
	double *mu;
	double *M; //Vector dels moments
	double *alpha;
	double *beta;
	double *gamma;
	double *delta;
	double *s; //Vector per avaluar els splines en cada xk
	double xk; //Valors per a estudiar l'error (donat a l'enunciat de l'apartat a)
    int k;
    double maxX; //Error màxim
	double f; //Per a avaluar la funció f(x) en xk
	
	printf("Introdueix el nombre n (els de l'enunciat del problema):\n");
    scanf("%d",&n);
	
	x=(double*)malloc((n+1)*sizeof(double));
	y=(double*)malloc((n+1)*sizeof(double));
	h=(double*)malloc(n*sizeof(double)); //El vector en sí va de 0 fins n-1 (com sempre en C), però l'entenem com de 1 fins a n-1 (guardem en les posicions 0 fins n-1 els valors de h per 1 fins n-1). El mateix amb d, lambda i mu
	d=(double*)malloc((n-1)*sizeof(double));
	lambda=(double*)malloc((n-1)*sizeof(double));
	mu=(double*)malloc((n-1)*sizeof(double));
	M=(double*)malloc((n+1)*sizeof(double));
	matriu = (double**) malloc ((n-1)*sizeof(double*));
	for(i=0;i<(n-1);i++)
	{
	 matriu[i] = (double*) malloc ((n-1)*sizeof(double));
	}
	vector = (double*) malloc (n*sizeof(double));
	alpha=(double*)malloc(n*sizeof(double));
	beta=(double*)malloc(n*sizeof(double));
	gamma=(double*)malloc(n*sizeof(double));
	delta=(double*)malloc(n*sizeof(double));
	s=(double*)malloc(n*sizeof(double));
	
	//Rellenem els vectors: nodes equispaiats xj, f(xj), hi, di, lambdai, mui
    for(i=0;i<=n;i++)
    {
	   x[i]=-1+i*(2/((double)n));
	   y[i]=1/(1+25*x[i]*x[i]);
    }
	for(i=0;i<=n-1;i++)
    {
	   h[i]=x[i+1]-x[i];
    }
	for(i=0;i<=n-2;i++)
    {
	   d[i]=(6/(h[i]+h[i+1]))*(((y[i+2]-y[i+1])/h[i+1])-((y[i+1]-y[i]/h[i])));
	   lambda[i]=h[i+1]/(h[i]+h[i+1]);
	   mu[i]=h[i]/(h[i]+h[i+1]);
    }
	
	
	//Trobem els moments M: fem la matriu i resolem el sistema
	//Emplenem la matriu
	for (i=0 ; i<=n-2 ; i++)
	{ 
		for ( j=0 ; j<=n-2 ; j++)
		{
		  matriu[i][j]=0;
		}
		for ( j=0 ; j<=n-2 ; j++)
		{
			if (j==i) {
			  matriu[i][j] = 2; //A la diagonal tot 2
			}
		
	                if (j+1==i) {
			  matriu[i][j] = mu[i]; //La diagonal de baix tot mu
			}
		        if (j-1==i) {
			  matriu[i][j] = lambda[i]; //La diagonal de dalt tot lambda
			}
		}
		vector[i]=d[i];
	}
	PintaMatriu(matriu, n-1);
	//Resolem el sistema:
	printf("\nVols trobar el resultat utilitzant el metode amb pivotatge?\n 1. Si\n 2. No\n");
	scanf("%d", &resposta);
	if (resposta==1) {
		resultat = MetodeComplet_AmbPivotatge (matriu, vector, n-1);
		if (resultat == NULL) {
			return 1;
		}
	}
	else if (resposta==2) {
		resultat = MetodeComplet_SensePivotatge (matriu, vector, n-1);
		if (resultat == NULL) {
			resultat = PivotatgeCondicional (matriu, vector, n-1);
			if (resultat == NULL) {
			return 1;
			}
		}
	}
	//Ara tenim el resultat guardat a "resultat", així que el posem a "M" i hi afegim els termes que falten per acabar d'emplenar el vector de moments:
	M[0]=0;
	M[n]=0;
	for(i=0;i<(n-1);i++)
    {
	   M[i+1]=resultat[i];
    }
	
	//Construïm els coeficients dels splines, ara que ja tenim els moments
	for(i=0;i<=n-1;i++)
    {
	   alpha[i]=y[i];
	   beta[i]=((y[i+1]-y[i])/h[i])-(((2*M[i]+M[i+1])/((double)6))*h[i]);
	   gamma[i]=M[i]/((double)2);
	   delta[i]=((M[i+1]-M[i])/((double)6*h[i]));
	   printf("Els coeficients de cada spline son: \nspline %d: \nalpha=%.16f, beta=%.16f, gamma=%.16f, delta=%.16f\n", i, alpha[i], beta[i], gamma[i], delta[i]);
    }
	
	
	//Estudiem l'error d'interpolació: |f(xk)-p(xk)| amb p l'spline que toqui en cada cas.
	printf("\n");
    printf("El valor absolut de la diferencia entre f(xk) i p(xk) per cada xk es (en forma de punts, on el primer nombre es xk i el segon el val abs de la diferencia):\n");
    xk=0;
    maxX=0; 
    for (k=0; k<=180; k++){
		xk=-0.989+k*0.011;
		f=1/(1+25*xk*xk); //f(xk)
		for (i=0; i<=n-1; i++){
			if (xk>=x[i] && xk<=x[i+1]){
				s[i]=alpha[i] + beta[i]*(xk-x[i]) + gamma[i]*((xk-x[i])*(xk-x[i])) + delta[i]*((xk-x[i])*(xk-x[i])*(xk-x[i]));
				printf("[%f,%f],", xk, fabs(f-s[i]));
				if ((fabs(f-s[i])) > maxX){
					maxX=fabs(f-s[i]);
				}
			}
		}
	}
    printf("\n\n");
    printf("El val abs de la diferencia maxima es %.16f\n", maxX);
    printf("\n");
	
	
	//Alliberem memòria:
	for (i=0; i<(n-1); i++){
		free(matriu[i]);
		matriu[i] = NULL;
	}
	free(matriu);
	matriu = NULL;
	free(vector);
	vector = NULL;
	free(x);
	x = NULL;
	free(y);
	y = NULL;
	free(h);
	h = NULL;
	free(d);
	d = NULL;
	free(lambda);
	lambda = NULL;
	free(mu);
	mu = NULL;
	free(M);
	M = NULL;
	free(resultat);
	resultat = NULL;
	free(alpha);
	alpha = NULL;
	free(beta);
	beta = NULL;
	free(gamma);
	gamma = NULL;
	free(delta);
	delta = NULL;
	free(s);
	s = NULL;
return 0;
}