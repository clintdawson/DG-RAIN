/* function to compute the values of the other fields of Channel and set the initial condition */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MeshAttributes.h"

double *zeta;
double *Qx;
double *Qy;

double *Fhat1dotn;
double *Fhat2dotn;
double *Fhat3dotn;
double *RHSZeta;
double *RHSQx;
double *RHSQy;

double *bzeta;
double *bQn;

const double g = 9.81;

// parameters and array for WD treatment
int *WD;
const double H0 = 1e-3;		// threshold below which an element is considered to be dry
const double VELZERO = 1e-2;	// threshold below which discharge won't be transferred from a dry node to a wet node

void initialize()
{
	zeta = malloc(3*NumEl*sizeof(double));
	Qx = malloc(3*NumEl*sizeof(double));
	Qy = malloc(3*NumEl*sizeof(double));
	RHSZeta = malloc(3*NumEl*sizeof(double));
	RHSQx = malloc(3*NumEl*sizeof(double));
	RHSQy = malloc(3*NumEl*sizeof(double));

	Fhat1dotn = malloc(3*NumEl*sizeof(double));
	Fhat2dotn = malloc(3*NumEl*sizeof(double));
	Fhat3dotn = malloc(3*NumEl*sizeof(double));

	bzeta = calloc(NumEdges,sizeof(double));
	bQn = calloc(NumEdges,sizeof(double));
	
	WD = malloc(NumEl*sizeof(int));

	for (int el=0; el < NumEl; ++el)
	{
		double NodalZeta[3];	
		for (int k=0; k<3; ++k)
		{
	
			// Emerged bump
		/*	double zetaval = 0.1;
			double zval = z[EltoVert[el*3+k]];	
			zetaval = fmax(zetaval, -zval);
			
			double H = zetaval + zval;
			if (H <= H0)
				H = H0;

			NodalZeta[k] = H - zval;
*/

/*			// Dam break on a dry bed
			double xval = x_coord[EltoVert[el*3+k]];
			double H;
			if (xval <= 5)
				H = 0.005;
			else
				H = H0;
			double zval = z[EltoVert[el*3+k]];
			NodalZeta[k] = H - zval;		
*/

			// Shinatro Parabolic bowl
/*			int nodeNum = EltoVert[el*3+k];
			double xval = x_coord[nodeNum];
			double yval = y_coord[nodeNum];
			double zval = z[nodeNum];
			double XC = 1;
			double YC = -0.41884;
			double alpha = 1e-7;

			double H = 1/(XC+YC) + alpha*(YC*YC-XC*XC)*(xval*xval+yval*yval)/((XC+YC)*(XC+YC));
			if (H < H0)
				H = H0;
			NodalZeta[k] = H - zval;	
*/

			//SW analytic sol parabolic bowl
	/*		double a = 1;
			double r0 = 0.8;
			double h0 = 0.1;
			double L = 4;
			double A = (a*a-r0*r0)/(a*a+r0*r0);
			int nodeNum = EltoVert[el*3+k];
			double xval = x_coord[nodeNum];
			double yval = y_coord[nodeNum];
			double zval = z[nodeNum];
			double rsq = (xval - L/2)*(xval - L/2) + (yval-L/2)*(yval-L/2);

			NodalZeta[k] = h0*(sqrt(1-A*A)/(1-A) - 1 - rsq/(a*a)*((1-A*A)/((1-A)*(1-A)) -1));			
			double H = NodalZeta[k] + zval;
			if (H < H0)
				NodalZeta[k] = H0 - zval;
	*/
			// MacDonald's short channel with smooth transition and shock
	/*		int nodeNum = EltoVert[el*3+k];
			double zval = z[nodeNum];
			double h_ex = 2.87912;
			double z_last = -0.0008;
			double hval = fmax(h_ex - z_last + zval,0);
			NodalZeta[k] = hval - zval;
	*/

			NodalZeta[k] = 0.1;
			Qx[el*3+k] = 0;
			Qy[el*3+k] =0;
		}
		zeta[el*3] = 0.5*(NodalZeta[0]+NodalZeta[1]);
		zeta[el*3+1] = 0.5*(NodalZeta[1]+NodalZeta[2]);
		zeta[el*3+2] = 0.5*(NodalZeta[2] + NodalZeta[0]);


		WD[el] = -1;
	}

}
