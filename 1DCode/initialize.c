/* function to compute the values of the other fields of Channel and set the initial condition */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"

double *A;
double *Q;

double *Fhat1L;
double *Fhat1R;
double *Fhat2L;
double *Fhat2R;
double *RHSA;
double *RHSQ;

double dt; 
double max_lamdba;
//const double nFriction = 0.0328;

const double g = 9.81;

#ifdef WDON
// For wetting-drying treatment
int *WD;
const double H0 = 1e-3;
const double VELZERO = 1e-1; 
//const double VELZERO = 1e-2; 
#endif


void initialize()
{
	A = malloc(NumNodes*sizeof(double));
	Q = malloc(NumNodes*sizeof(double));
	RHSA = calloc(NumNodes,sizeof(double));
	RHSQ = calloc(NumNodes,sizeof(double));
	//Fhat1L = calloc(NumNodes,sizeof(double));
	//Fhat1R = calloc(NumNodes,sizeof(double));
	//Fhat2L = calloc(NumNodes,sizeof(double));
	//Fhat2R = calloc(NumNodes,sizeof(double));
	
	#ifdef WDON
	WD = malloc(NumEl*sizeof(int));
	#endif

	for(int i = 1; i < NumNodes-1; i++)
	{
		double h; 
		double zval = Nodalz[i-1];
		double zeta = 2.0;
		//h = 1.35;
		h = zeta+zval;
		double bval = NodalB[i-1];
		double m1val = Nodalm1[i-1];
		double m2val = Nodalm2[i-1];
		A[i] = h*bval + 0.5*m1val*h*h + 0.5*m2val*h*h;
		Q[i] = 0;
	
	}
	A[0] = A[1];
	A[NumNodes-1] = A[NumNodes-2];
	Q[0] = Q[1];
	Q[NumNodes-1] = Q[NumNodes-2]; 


/*	for (int i =0; i < NumNodes; ++i)
	{
		for (int j =0; j < 2; ++j)
		{	
			double h;
			double zval = z[i];
			double zeta = 0.33;
			h = zeta+zval;

			// Emerged bump
			double zeta = 0.1;
			zeta = fmax(zeta, -zval);
			h = zeta + zval; 
*/
/*			// Parabolic bowl
			double xval = x[i];
			if ((xval > 0.5) && (xval < 2.5))
				h = -0.625+1.5*xval-0.5*xval*xval;
			else
				h = H0;
*/
			// Dam break on a dry bed
/*			if ( xval <= 5)
				h = 0.005;
			else
				h = H0;
*/		

			// MacDonald's short channel with smooth transition and shock
/*			double h_ex = 2.87912;
			double hval = fmax(h_ex -z[NumNodes-1] + zval, 0);

			double bval = b[i];
			A[i*2+j] = h*bval;
			Q[i*2+j] = 0;

			RHSA[i*2+j] = 0;
			RHSQ[i*2+j] = 0;
		}
	
		#ifdef WDON
		if (i < NumEl)
			WD[i] = -1;
		#endif

	}
*/
}
