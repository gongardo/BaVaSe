#include "R.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <R_ext/Utils.h>
#include<time.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>

#include "allBF.h"
//#include "allBF.c"


void eZSBF (double *pg, int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	double g=*pg;
	int k2=*pk2;
	int k0=*pk0;			
	double Q=*pQ;
	*B21=eZSBF21fun(g, n, k2, k0, Q);

}

