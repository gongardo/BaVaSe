#include <R.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


/*-------JZS-------*/


/*auxiliar functions for integration */

/*Structure for Parameters,S*/
struct par {
	double g;
	double n;
	double k_i;
	double k_0;
	double Q_i0;
};

double ezell_aux (double x, void *p){
	struct par * params=(struct par *)p;/*Defino un puntero a una estructura del tipo par*/
	/*Inicializo el puntero en la dirección de memoria en la que están los parametros que le estoy
	 pasando, p obligando a que sean de tipo struct par*/
	
	/*Defino los parametros que son los que estaran en la estructura que le pasamos*/
	double g=(params->g);
	double n=(params->n);
	double k=(params->k_i);/*it will be k2*/
	double kk0=(params->k_0);	
	double Q=(params->Q_i0);
	
	/*Calculo el valor de la función y lo devuelvo*/
	//double l=pow((1.0+x), (n-k)/2.0)*pow((1.0+Q*x), (kk0-n)/2)*pow((n/(2.0*M_PI)),0.5)*pow(x, -1.5)*exp(-n/(2.0*x));
	double l=exp(0.5*(n-k)*log(1.0+x) + 0.5*(kk0-n)*log(1.0+Q*x) + 0.5*log(g/(2.0*M_PI)) - 1.5*log(x) - g/(2.0*x));
	
	return l;
}


/*Integrated functions the arguments will be n,k,Qi0*/
double ezell (double g, double n, double k, double k0, double Q){
	/*allocate space por integration*/
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0, error=0.0;
	
	/*set parameters in the appropiate structure*/
	struct par  params={g, n, k, k0, Q};
	
	/*define the function and pass parameters*/
	gsl_function F;
	F.function = &ezell_aux;
	F.params = &params;
	
	/*integrate and save result and error*/
	gsl_integration_qagiu(&F, 0, 0, 1e-9,10000,w,&result,&error);
	
	/*free space*/
	gsl_integration_workspace_free (w);
	
	return result;
}



/* FUNCION QUE USAREMOS EN EL main.c*/
double eZSBF21fun(double g, int n, int k2, int k0, double Q)
{
    double ZSBF21 = ezell (g, (double) n, (double) k2, (double) k0, Q);
    if (!R_FINITE(ZSBF21)){error("A Bayes factor is infinite.");}
    return ZSBF21;
    
}

