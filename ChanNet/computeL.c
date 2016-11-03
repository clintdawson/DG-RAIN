#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ChannelsAndJunctions.h"
#include "mathfunctions.h"

/***************************************************************************************//**
* @file computeL.c
*
* This file contains code to evaluate the right hand side of the following discrete ODE
* obtained from the DG discretization of the 1-D St. Venant equations:
*
* \f$\frac{\partial \hat{\mathbf{w}}}{\partial t} = M^{-1} L(\hat{\mathbf{w}}, t)\f$
*
********************************************************************************************/

/* @cond FUNCTION_PROTOTYPES */
extern const double g;
extern void Compute1DInnerProducts(struct channel *Chan, int el, double *IP);
extern void RoeFlux1D(double A_L, double A_R, double Q_L, double Q_R, double b, double g, double beta, double *Fhat);
extern void LF(double A_L, double A_R, double Q_L, double Q_R, double b, double g, double *Fhat);
extern double getBeta(struct channel *Chan, int edge, int channelNumber, double time);
/* @endcond */

/***************************************************************************************//**
* Function for evaluating the right hand side of the discrete ODE obtained from the DG 
* discretization of the 1-D St. Venant equations
* @param [in] Chan a pointer to the channel structure corresponding to the channel that is 
* currently being worked on
* @param [in] time a double representing the current time of the simulation
* @param [in] channelNumber an integer number designated for the channel that is currently
* being worked on. This is only there for debugging purposes and we might not need it anymore.
* @param [out] RHSA a pointer to an array of size NumEl x 2 in which the right hand side for
* the wet cross-sectional area will be stored
* @param [out] RHSQ a pointer to an array of size NumEl x 2 in which the right hand side for
* the volumetric discharge will be stored
*
* *********************************************************************************************/

void computeL(struct channel *Chan, double time, int channelNumber, double *RHSA, double *RHSQ)
{
		
/************** Compute the flux (Roe's flux) at the faces *****************/
	int NumNodes = Chan->NumNodes;
	int NumEl = Chan->NumEl;	

	double Fhat1L[NumNodes];
	double Fhat2L[NumNodes];

	double Fhat1R[NumNodes];
	double Fhat2R[NumNodes];

	double beta = 1;

	//getBeta(Chan, NumNodes-1,channelNumber, time);

	for (int i=0; i < NumNodes; ++i)
	{
		double A_L, Q_L, A_R, Q_R, b;
		b = Chan->b[i];

		double tmpF[2];
		
		#ifdef WDON
		// Check to see if the elements separated by this boundary are both dry
		if (i > 0 && i < NumNodes-1)
		{
			if (Chan->WD[i-1] == 0 && Chan->WD[i] == 0)
			{
				// Reflection flux for the left element
				A_L = Chan->A[2*i];
				A_R = A_L;
				Q_L = Chan->Q[2*i];
				Q_R = -Q_L;
				RoeFlux1D(A_L, A_R, Q_L, Q_R, b,0,beta,tmpF);
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];
			
				// Reflection flux for the right element
				A_R = Chan->A[2*i+1];
				A_L = A_R;
				Q_R = Chan->Q[2*i+1];
				Q_L = -Q_R;
				RoeFlux1D(A_L, A_R, Q_L, Q_R, b, 0,beta,tmpF);
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];

				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("1D both element dry flux not a number, edge %d\n",i);
					exit(EXIT_FAILURE);
				}
			}
		}

		// if the elements are not both dry
		if (i == 0 || i == NumEl || Chan->WD[i-1] == 1 || Chan->WD[i] == 1)
		{
			A_L = Chan->A[2*i];
			A_R = Chan->A[2*i+1];

			Q_L = Chan->Q[2*i];
			Q_R = Chan->Q[2*i+1];

			RoeFlux1D(A_L, A_R, Q_L, Q_R, b, g,beta,tmpF);
			
			if (isnan(tmpF[0]) || isnan(tmpF[1]))
			{
				printf("1D both elements wet flux not a number, edge %d\n", i);
				exit(EXIT_FAILURE);
			}

			if (i == 0 || Chan->WD[i-1] == 1)
			{
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];
			}
			else
			{
				double newtmpF[2];
				RoeFlux1D(A_L, A_R, Q_L, Q_R, b, 0,beta, newtmpF);
				Fhat1L[i] = newtmpF[0];
				Fhat2L[i] = newtmpF[1];
				if (isnan(newtmpF[0]) || isnan(newtmpF[1]))
				{
					printf("1D Left element dry flux not a number, edge %d\n", i);
					exit(EXIT_FAILURE);
				}
			}
			
			if (i == NumEl || Chan->WD[i] == 1)
			{
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];
			}
			else
			{
				RoeFlux1D(A_L, A_R, Q_L, Q_R, b, 0, beta, tmpF);
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];
				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("1D right element dry flux not a number, edge %d\n", i);
					exit(EXIT_FAILURE);
				}
			}
		}

		#else
		A_L = Chan->A[2*i];
		A_R = Chan->A[2*i+1];

		Q_L = Chan->Q[2*i];
		Q_R = Chan->Q[2*i+1];

		//beta = getBeta(Chan,i,channelNumber,time);
		//Chan->beta = 1.5;
		//double beta = Chan->beta[i];

		//	printf("beta = %e\t channelNumber = %d\n",beta, channelNumber);

		RoeFlux1D(A_L, A_R, Q_L, Q_R, b, g,beta,tmpF);
		//LF(A_L, A_R, Q_L, Q_R, b, g, tmpF);

		if (isnan(tmpF[0]) || isnan(tmpF[1]))
		{
			printf("1D numerical flux not a number, channelNumber %d  edge %d\n", channelNumber,i);
			printf("A_L = %e, A_R = %e, Q_L = %e, Q_R = %e \n", A_L, A_R, Q_L, Q_R);
			exit(EXIT_FAILURE);
		}

		Fhat1L[i] = tmpF[0];
		Fhat2L[i] = tmpF[1];

		Fhat1R[i] = tmpF[0];
		Fhat2R[i] = tmpF[1];

		#endif

		/*****************************************************************************/
		double u_L = Q_L/A_L;
		double u_R = Q_R/A_R;
		double c_L = sqrt(g*A_L/b);
		double c_R = sqrt(g*A_R/b);

		if (i==0)
			Chan->max_lambda = fmax((fabs(u_L)+c_L), (fabs(u_R)+c_R));			
		// Compute maximum eigenvalue for the next time step
		double current_max =fmax((fabs(u_L) + c_L), (fabs(u_R) + c_R));
		Chan->max_lambda = max((Chan->max_lambda),(current_max));
	}
		
	// Compute the right hand side

	for (int k=0; k < NumEl; ++k) 
	{
		double IP[10];
		Compute1DInnerProducts(Chan, k, IP);
			
		double x1 = Chan->x[k];
		double x2 = Chan->x[k+1];
		double y1 = Chan->y[k];
		double y2 = Chan->y[k+1];
		
		double h = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	
		double F_el1[2] = {Fhat1R[k], Fhat1L[k+1]};
		double F_el2[2] = {Fhat2R[k], Fhat2L[k+1]};

		double invMrow1[2] = {4/h, -2/h};
		double invMrow2[2] = {-2/h, 4/h};
		//double invM[2][2] = {{4/(b-a), 2/(a-b)},{2/(a-b), 4/(b-a)}};

		double S1 = IP[0] + F_el1[0];
		double S2 = IP[1] - F_el1[1];
		double S3 = IP[2] + F_el2[0]  + IP[3] + IP[4] - IP[5];
		double S4 = IP[6] - F_el2[1] + IP[7] + IP[8] - IP[9];
	
	
		// Compute RHSA = invM*[S1;S2] and RHSQ = invM*[S3;S4]
		RHSA[2*k + 1] = invMrow1[0]*S1 + invMrow1[1]*S2;  
		RHSA[2*k + 2] = invMrow2[0]*S1 + invMrow2[1]*S2;
		RHSQ[2*k + 1] = invMrow1[0]*S3 + invMrow1[1]*S4;
		RHSQ[2*k + 2] = invMrow2[0]*S3 + invMrow2[1]*S4;

	}

	// The value at the ghost nodes will be the same as the values at the other side of the face

	RHSA[0] = RHSA[1];
	RHSA[2*Chan->NumNodes - 1] = RHSA[2*Chan->NumNodes - 2];

	RHSQ[0] = RHSQ[1];
	RHSQ[2*Chan->NumNodes - 1] = RHSQ[2*Chan->NumNodes - 2];
		

}


