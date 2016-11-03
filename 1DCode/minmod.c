#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"
#include "constitutive_equations.h"
#include "globals.h"

extern void GetLGLWeights(int N, double* w);
void minmod()
{
	
	// average value of Zeta and Q
	double avgZeta[NumEl];
	double avgQ[NumEl];
	double Zeta[NumEl][Np];

	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i=0; i < NumEl; ++i)
	{
		double avgZetaVal = 0; 
		double avgQVal = 0;
		for (int j = 0; j < Np; j++)
		{
			int nodeNum = i*Np+j;
			double Aval = A[nodeNum+1];
			double bval = NodalB[nodeNum];
			double zval = Nodalz[nodeNum];
			double m1val = Nodalm1[nodeNum];
			double m2val = Nodalm2[nodeNum];
			double Hval = getH(Aval, bval, m1val, m2val);
			double Qval = Q[nodeNum+1];
			double zetaVal = Hval - zval;
			Zeta[i][j] = zetaVal;

			avgZetaVal += LGLWeight[j]*zetaVal;
			avgQVal += LGLWeight[j]*Qval;
		}
		avgZeta[i] = 0.5*avgZetaVal;
		avgQ[i] = 0.5*avgQVal;

	}

	double a1[NumEl], a2[NumEl], b1[NumEl], b2[NumEl];
	for (int i =0; i < NumEl; ++i)
	{	
		if (i > 0)
		{
			a1[i] = 2*(avgZeta[i] - avgZeta[i-1])/(dh[i] + dh[i-1]);
			a2[i] = 2*(avgQ[i] - avgQ[i-1])/(dh[i]+dh[i-1]);

		}
		if (i < NumEl-1)
		{	
			b1[i] = 2*(avgZeta[i+1] - avgZeta[i])/(dh[i]+dh[i+1]);
			b2[i] = 2*(avgQ[i+1] - avgQ[i])/(dh[i]+dh[i+1]);
		}


	}

	a1[0] = 0;
	a2[0] = 0;
	b1[NumEl-1] = 0;
	b2[NumEl-1] = 0;

	double sigma1[NumEl], sigma2[NumEl];	
	for (int i =0; i < NumEl; ++i)
	{
		sigma1[i] = sign(a1[i])*0.5*(1+sign(a1[i]*b1[i]))*fmin(fabs(a1[i]),fabs(b1[i]));
		sigma2[i] = sign(a2[i])*0.5*(1+sign(a2[i]*b2[i]))*fmin(fabs(a2[i]),fabs(b2[i]));

	}

	double Zetatilda[NumEl][2];
	double Qtilda[NumEl][2];

	for (int i = 0 ; i < NumEl; ++i)
	{
		Zetatilda[i][0] = avgZeta[i] - 0.5*sigma1[i]*dh[i];
		Zetatilda[i][1] = avgZeta[i] + 0.5*sigma1[i]*dh[i];
		Qtilda[i][0] = avgQ[i] - 0.5*sigma2[i]*dh[i];
		Qtilda[i][1] = avgQ[i] + 0.5*sigma2[i]*dh[i];
	}

	double eps0 = 1e-8;
	for (int i = 0; i < NumEl; ++i)
	{
		// limit only if limiting is required
		if (fabs(Zetatilda[i][0]-Zeta[i][0]) > eps0 || fabs(Zetatilda[i][1] - Zeta[i][Np-1]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Zeta[i][j] = Zetatilda[i][0] + j*(Zetatilda[i][1] - Zetatilda[i][0])/(Np-1);
				double hval = Zeta[i][j] + Nodalz[i*Np+j];
				double bval = NodalB[i*Np+j];

				double m1val = Nodalm1[i*Np+j];
				double m2val = Nodalm2[i*Np+j];
				A[i*Np+j+1] = hval*bval + 0.5*m1val*hval*hval+ 0.5*m2val*hval*hval;
			}
		}
	
		// limit only if limiting is required
		if (fabs(Qtilda[i][0]-Q[i*Np+1]) > eps0 || fabs(Qtilda[i][1] - Q[i*Np+Np]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Q[i*Np+j+1] = Qtilda[i][0] + j*(Qtilda[i][1] - Qtilda[i][0])/(Np-1);
			}
			
		}	

	}	


/*	for (int i = 0; i < NumNodes; ++i)
	{
		int modify;
		#ifdef WDON
		modify = (height[i][0] > H0)*(height[i][1]>H0)*((Zeta[i][0]+zEl[i][0]) > 0)*((Zeta[i][1]+zEl[i][1])>0);
		#else
		modify = 1;
		#endif

		for (int k =0; k < 2; ++k)
		{
			if (modify)
			{
				double zeta = Zeta[i][k];
				height[i][k] = zeta + zEl[i][k];
				A[i*2+k+1] = height[i][k] * bEl[i][k];
				Q[i*2+k+1] = QEl[i][k];
			}
		}
	}
*/
	A[0] = A[1];
	Q[0] = Q[1];
	A[NumNodes-1] = A[NumNodes-2];
	Q[NumNodes-1] = Q[NumNodes-2];


}
