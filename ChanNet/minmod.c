/***********************************************************************//**
* @file minmod.c
*
* This file contains code to apply the minmod slope limiter on the water
* surface height and the volumetric discharge obtained from the 1-D 
* RKDG scheme.
*
* **************************************************************************/

#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"
#include "ChannelsAndJunctions.h"

#ifdef WDON
extern const double H0;
#endif


/**********************************************************************//**
*
* Function to apply the minmod slope limiter on the 1-D conserved variables
* @param [in] Chan a pointer to the channel structure that we are currently
* working on
*
* ************************************************************************/
void minmod(struct channel *Chan)
{
	int TotalNumNodes = 2*Chan->NumNodes;
	int NumEl = Chan->NumEl;

	double height[NumEl][2];
	double zEl[NumEl][2];
	double Zeta[NumEl][2];
	double QEl[NumEl][2];
	double bEl[NumEl][2];
	double dh[NumEl];

	for (int i=0; i < NumEl; ++i)
	{
		double Aval1 = Chan->A[2*i+1];
		double Aval2 = Chan->A[2*i+2];
		double bval1 = Chan->b[i];
		double bval2 = Chan->b[i+1];
		height[i][0] = Aval1/bval1;
		height[i][1] = Aval2/bval2;
		zEl[i][0] = Chan->z[i];
		zEl[i][1] = Chan->z[i+1];
		bEl[i][0] = Chan->b[i];
		bEl[i][1] = Chan->b[i+1];	
		Zeta[i][0] = height[i][0]-zEl[i][0];
		Zeta[i][1] = height[i][1]-zEl[i][1];
		QEl[i][0] = Chan->Q[2*i+1];
		QEl[i][1] = Chan->Q[2*i+2];
	
		double x1 = Chan->x[i];
		double x2 = Chan->x[i+1];
		double y1 = Chan->y[i];
		double y2 = Chan->y[i+1];
		dh[i] = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)); 


	}

	// average value of Zeta and Q
	double avgZeta[NumEl];
	double avgQ[NumEl];
	for (int i =0; i < NumEl; ++i)
	{
		avgZeta[i] = (Zeta[i][0]+Zeta[i][1])/2;
		avgQ[i] = (QEl[i][0]+QEl[i][1])/2;
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
		if (fabs(Zetatilda[i][0]-Zeta[i][0]) > eps0 || fabs(Zetatilda[i][1] - Zeta[i][1]) > eps0)
		{
			Zeta[i][0] = Zetatilda[i][0];
			Zeta[i][1] = Zetatilda[i][1];
		}
	
		// limit only if limiting is required
		if (fabs(Qtilda[i][0]-QEl[i][0]) > eps0 || fabs(Qtilda[i][1] - QEl[i][1]) > eps0)
		{
			QEl[i][0] = Qtilda[i][0];
			QEl[i][1] = Qtilda[i][1];
		}	

	}	

	for (int i = 0; i < Chan->NumNodes; ++i)
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
				Chan->A[i*2+k+1] = height[i][k] * bEl[i][k];
				Chan->Q[i*2+k+1] = QEl[i][k];
			}
		}
	}

	Chan->A[0] = Chan->A[1];
	Chan->Q[0] = Chan->Q[1];
	Chan->A[TotalNumNodes-1] = Chan->A[TotalNumNodes-2];
	Chan->Q[TotalNumNodes-1] = Chan->Q[TotalNumNodes-2];


}
