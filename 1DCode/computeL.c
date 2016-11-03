#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "mathfunctions.h"
#include "globals.h"
#include "constitutive_equations.h"

/*****************************************************************************/
/***************************************************************************/
/************ Compute the spatial discretization of the equation ***********/
/***************************************************************************/


extern const double g;
extern void Compute_InnerProducts(double IP1[], double IP2[], double IP3[], double IP4[], double IP5[], double IP6[], double IP7[], double IP8[], double IP9[], double IP10[]);
extern void RoeFlux( double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double m1val, double m2val, double g);
extern void LF(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double g);

extern void* xcalloc(int items, int size);
extern void GetLGLWeights(int N, double *w);

void computeL()
{
	double* Fhat1L = xcalloc(NumNodes, sizeof(double));
	double* Fhat2L = xcalloc(NumNodes, sizeof(double));
	double* Fhat1R = xcalloc(NumNodes, sizeof(double));
	double* Fhat2R = xcalloc(NumNodes, sizeof(double));


#ifdef WDON
/**** Compute the mass over an element to use later for ensuring ********/
/**** positive mass with the flux in a wetting and drying treatment *************/
	double mass[NumEl];
	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i = 0; i < NumEl; ++i)
	{
		double avgArea = 0;
		int begNode = i*Np + 1;
		for (int j =0; j < Np; j++)
		{
			avgArea += LGLWeight[j]*A[begNode+1];
		}
		mass[i] = avgArea;
	}
#endif
		
/************** Compute the numerical flux at the faces *****************/
	for (int i=0; i < NumEdges; ++i)
	{
		double A_L, Q_L, A_R, Q_R, m1val, m2val;

		int leftNode = i*Np;
		int rightNode = leftNode+1;
	
		m1val = m1[i];
		m2val = m2[i];

		double tmpF[2];
		
		#ifdef WDON
		// Check to see if the elements separated by this boundary are both dry
		if (i > 0 && i < NumNodes-1)
		{
			if (WD[i-1] == 0 && WD[i] == 0)
			{
				// Reflection flux for the left element
				A_L = A[leftNode];
				A_R = A_L;
				Q_L = Q[leftNode];
				Q_R = -Q_L;
				//LF(tmpF, A_L, A_R, Q_L, Q_R,b[i], m1val, m2val, 0);
				RoeFlux(tmpF, A_L, A_R, Q_L, Q_R,b[i], m1val, m2val,0);
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];
				
				// Reflection flux for the right element
				A_R = A[rightNode];
				A_L = A_R;
				Q_R = Q[rightNode];
				Q_L = -Q_R;
				//LF(tmpF, A_L, A_R, Q_L, Q_R, b[i], m1val, m2val, 0);
				RoeFlux(tmpF, A_L, A_R, Q_L, Q_R, b[i], m1val, m2val, 0);
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];

				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("both elements dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
				
			}

		}

		// if the elements are not both dry
		if ((i == 0) || (i == NumEl) || WD[i-1] == 1 || WD[i] == 1)
		{
			A_L = A[leftNode];
			A_R = A[rightNode];
	
			Q_L = Q[leftNode];
			Q_R = Q[rightNode];

			RoeFlux(tmpF, A_L, A_R, Q_L,Q_R,b[i], m1val, m2val, g);
			//LF(tmpF, A_L, A_R, Q_L,Q_R,b[i], m1val, m2val, g);
	
			if (isnan(tmpF[0]) || isnan(tmpF[1]))
			{
				printf(" both elements wet flux not a number, edge %d \n",i);
				exit(EXIT_FAILURE);
			}
			
			if (i==0 || WD[i-1] == 1)
			{
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];	

			}
			else
			{
				double newtmpF[2];
				//LF(newtmpF, A_L, A_R, Q_L, Q_R, b[i], m1val, m2val, 0);
				RoeFlux(newtmpF, A_L, A_R, Q_L, Q_R, b[i], m1val, m2val, 0);
				Fhat1L[i] = newtmpF[0];
				Fhat2L[i] = newtmpF[1];
				if (isnan(newtmpF[0]) || isnan(newtmpF[1]))
				{
					printf("left element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
	
			}
			if (i == NumEl || WD[i]==1)
			{
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];
			}
			else
			{
				//LF(tmpF, A_L,A_R,Q_L,Q_R,b[i], m1val, m2val, 0);
				RoeFlux(tmpF, A_L,A_R,Q_L,Q_R,b[i], m1val, m2val, 0);
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];		
				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("right element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
	
			}
		}

		#else
		A_L = A[leftNode];
		A_R = A[rightNode];
	
		Q_L = Q[leftNode];
		Q_R = Q[rightNode];

		RoeFlux(tmpF, A_L, A_R, Q_L,Q_R,b[i],m1val,m2val,g);
		//LF(tmpF, A_L, A_R, Q_L,Q_R,b[i],m1val,m2val,g);
	
		if (isnan(tmpF[0]) || isnan(tmpF[1]))
		{
			printf("flux not a number, edge %d \n",i);
			printf("A_L = %lf A_R = %lf Q_L = %lf Q_R = %lf, b = %lf, m1 = %lf, m2 = %lf \n", A_L, A_R, Q_L, Q_R, b[i], m1val, m2val);
			exit(EXIT_FAILURE);
		}
			
		Fhat1L[i] = tmpF[0];
		Fhat2L[i] = tmpF[1];	

		Fhat1R[i] = tmpF[0];
		Fhat2R[i] = tmpF[1];
		
		#endif

		/**************************************************************************************************/
		double u_L = Q_L/A_L;
		double u_R = Q_R/A_R;
		double c_L = sqrt(g*A_L/b[i]);
		double c_R = sqrt(g*A_R/b[i]);

		if (i==0)
			max_lambda = fmax((fabs(u_L)+c_L), (fabs(u_R)+c_R));			
		// Compute maximum eigenvalue for the next time step
		double current_max =fmax((fabs(u_L) + c_L), (fabs(u_R) + c_R));
		max_lambda = max((max_lambda),(current_max));
	}

	printf("max_lambda = %lf\n", max_lambda);
	#ifdef WDON
	for (int i =1; i < NumNodes-1; ++i)
	{
		int leftNode = i*Np;
		int rightNode = leftNode+1;
		double maxBetaOverAlpha = 2;
		// Check to see if this flux might possibly result in negative mass
		// If the mass will be negative on the left side
	      	if (Fhat1L[i]*maxBetaOverAlpha*dt > mass[i-1])
		{
	    double A_L = A[leftNode];
			double A_R = A_L;
			double Q_L = Q[leftNode];
	    double Q_R = -Q_L;
	    double tmpF[2];
			double localG = g;
			
			
			if (WD[i] == 0)
				localG = 0;
			
			//LF(tmpF, A_L, A_R, Q_L, Q_R, b[i],m1val,m2val,localG);
			RoeFlux(tmpF, A_L, A_R, Q_L, Q_R, b[i],m1val,m2val,localG);
	      		Fhat1L[i] = tmpF[0];
	      		Fhat2L[i] = tmpF[1];
			printf("element %d will be dry\n",i);
	    	}
		
	
		// If the mass will be negative on the right side
		if (-Fhat1R[i]*maxBetaOverAlpha*dt > mass[i])
		{
			double A_R = A[rightNode];
			double A_L = A_R;
			double Q_R = Q[rightNode];
			double Q_L = -Q_R;
			double tmpF[2];
			double localG = g;
			if (WD[i] == 0)
				localG = 0;
			//LF(tmpF, A_L, A_R, Q_L, Q_R, b[i],m1val,m2val,localG);
			RoeFlux(tmpF, A_L, A_R, Q_L, Q_R, b[i],m1val,m2val,localG);
			Fhat1R[i] = tmpF[0];
			Fhat2R[i] = tmpF[1];
			printf("element %d will be dry \n",i);
				
		}
	}
	#endif

	// Compute the right hand side
//	double IP1[NumEl], IP2[NumEl], IP3[NumEl], IP4[NumEl], IP5[NumEl], IP6[NumEl], IP7[NumEl], IP8[NumEl], IP9[NumEl], IP10[NumEl];	
//	Compute_InnerProducts(IP1, IP2, IP3, IP4, IP5, IP6, IP7, IP8, IP9, IP10);

	int k;
	for (k=0; k < NumEl; ++k) 
	{
		
		double h = dh[k];	
		gsl_vector *F1 = gsl_vector_alloc(Np);
		gsl_vector *F2 = gsl_vector_alloc(Np);
		gsl_vector *ST21 = gsl_vector_alloc(Np);
		gsl_vector *ST22 = gsl_vector_alloc(Np);
		gsl_vector *ST23 = gsl_vector_alloc(Np);

		int begNode = k*Np;
		for (int i = 0; i < Np; i++)
		{
			double Aval = A[begNode+i+1];
			double Qval = Q[begNode+i+1];
			double bval = NodalB[begNode+i];
			double S0 = dz[begNode+i];
			double m1val = Nodalm1[begNode+i];
			double m2val = Nodalm2[begNode+i];
			double dm1val = dm1[begNode+i];
			double dm2val = dm2[begNode+i];
			double hval = getH(Aval, bval, m1val, m2val);
			double I1val = getI1(Aval, bval, m1val, m2val);
			double I2val = getI2(Aval, bval, db[k*Np+i], m1val, dm1val, m2val, dm2val);
			double Sfval = getS_f(Aval, Qval, bval, m1val, m2val, NodalnFriction[begNode+i]);
			double beta = 1.0;
			double localG;
			#ifdef WDON
			if(WD[k]==1)
				localG = g;
			else
				localG = 0;
			#else
				localG = g;
			#endif
			
			gsl_vector_set(F1, i, Qval);
			gsl_vector_set(F2, i, beta*Qval*Qval/Aval + g*I1val); 
			gsl_vector_set(ST21, i, g*I2val);
			//gsl_vector_set(ST22, i, g*hval*bval*S0);
			gsl_vector_set(ST22, i, g*Aval*S0);
			gsl_vector_set(ST23, i, g*Aval*Sfval);

		}

		gsl_vector *localFhat1 = gsl_vector_alloc(2);
		gsl_vector_set(localFhat1, 0, Fhat1R[k]);
		gsl_vector_set(localFhat1, 1, Fhat1L[k+1]);

		gsl_vector *localFhat2 = gsl_vector_alloc(2);
		gsl_vector_set(localFhat2, 0, Fhat2R[k]);
		gsl_vector_set(localFhat2, 1, Fhat2L[k+1]);

		// cacluate the volume integral of the flux
		gsl_vector *localRHS1 = gsl_vector_calloc(Np);
		gsl_vector *localRHS2 = gsl_vector_calloc(Np);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, VolMat, F1, 1.0, localRHS1);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, VolMat, F2, 1.0, localRHS2);

		// calculate the surface integral of the flux
		gsl_vector *SurfPart1 = gsl_vector_calloc(Np);
		gsl_vector *SurfPart2 = gsl_vector_calloc(Np);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, LIFT, localFhat1, 1.0, SurfPart1);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, LIFT, localFhat2, 1.0, SurfPart2);

		// calculate the RHS
		gsl_vector_add(localRHS1, SurfPart1);
		gsl_vector_add(localRHS2, SurfPart2);
		gsl_vector_add(localRHS2, ST21);
		gsl_vector_add(localRHS2, ST22);
		gsl_vector_sub(localRHS2, ST23);

		gsl_vector_free(F1);
		gsl_vector_free(F2);
		gsl_vector_free(localFhat1);
		gsl_vector_free(localFhat2);
		gsl_vector_free(SurfPart1);
		gsl_vector_free(SurfPart2);
		gsl_vector_free(ST21);
		gsl_vector_free(ST22);
		gsl_vector_free(ST23);
	
		for(int i = 0; i < Np; i++)
		{
			RHSA[begNode+i+1] = gsl_vector_get(localRHS1, i);
			RHSQ[begNode+i+1] = gsl_vector_get(localRHS2, i);
		}

		gsl_vector_free(localRHS1);
		gsl_vector_free(localRHS2);

	/*	
		double F_el1[2] = {Fhat1R[k], Fhat1L[k+1]};
		double F_el2[2] = {Fhat2R[k], Fhat2L[k+1]};

		double invMrow1[2] = {4/h, -2/h};
		double invMrow2[2] = {-2/h, 4/h};
		//double invM[2][2] = {{4/(b-a), 2/(a-b)},{2/(a-b), 4/(b-a)}};

		double S1 = IP1[k] + F_el1[0];
		double S2 = IP2[k] - F_el1[1];
		double S3 = IP3[k] + F_el2[0]  + IP4[k] + IP5[k] - IP6[k];
		double S4 = IP7[k] - F_el2[1] + IP8[k] + IP9[k] - IP10[k];
	
		// Compute RHSA = invM*[S1;S2] and RHSQ = invM*[S3;S4]
		RHSA[2*k + 1] = invMrow1[0]*S1 + invMrow1[1]*S2;  
		RHSA[2*k + 2] = invMrow2[0]*S1 + invMrow2[1]*S2;
		RHSQ[2*k + 1] = invMrow1[0]*S3 + invMrow1[1]*S4;
		RHSQ[2*k + 2] = invMrow2[0]*S3 + invMrow2[1]*S4;
*/
	}

	// The value at the ghost nodes will be the same as the values at the other side of the face

	RHSA[0] = RHSA[1];
	RHSA[NumNodes - 1] = RHSA[NumNodes - 2];

	RHSQ[0] = RHSQ[1];
	RHSQ[NumNodes - 1] = RHSQ[NumNodes - 2];

	free(Fhat1L);
	free(Fhat1R);
	free(Fhat2L);
	free(Fhat2R);

}


