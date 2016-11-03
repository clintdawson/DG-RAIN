/*******************************************************************************//**
*
* @file numericalFlux1D.c
*
* This file contains code to compute different types of numerical fluxes for the
* 1-D St. Venant equations.
***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"

/* @cond FUNCTION_PROTOTPYES */ 
extern double g;
extern double getI1(double, double);
extern double getI2(double, double, double);

/* @endcond */

/****************************************************************************************//**
*
* Function to calculate Roe's flux for the 1-D St. Venant equations at an interface (node).
* @param[in] A_L value of the wet cross-sectional area at the node coming from the 
* left element.
* @param [in] A_R value of the wet cross-sectional area at the node coming from the right
* element
* @param [in] Q_L value of the volumetric discharge at the node coming from the left element
* @param [in] Q_R value of the volumetric discharge at the node coming from the right element
* @param [in] b value of the width of the channel at the node. We assume that the width of the 
* channel is prescribed at each node and thus the value will be continuous at nodes
* @param [in] localG value of the gravitational acceleration constant. This will be the same
* as g, except when we set it to zero to handle wetting and drying.
* @param [in] beta value of the momentum correction coefficient at the node
* @param [out] Fhat a pointer to an array of size 2 where the numerical flux computed will be stored
*
* ************************************************************************************/
void RoeFlux1D(double A_L, double A_R, double Q_L, double Q_R, double b, double localG, double beta, double *Fhat)
{
	double I1_L, I1_R;
	double F1_L, F1_R;
	double F2_L, F2_R;
	double c_L, c_R;
	double u_L, u_R;
	double lambda1_L, lambda1_R, lambda2_L, lambda2_R;
	double uhat, chat, alpha1, alpha2;
	double lambdahat1, lambdahat2;
	double epsilon1, epsilon2;

	I1_L = getI1(A_L, b);
	I1_R = getI1(A_R, b);

	F1_L = Q_L;
	F1_R = Q_R;

	F2_L = pow(Q_L,2)/A_L + localG*I1_L;
	F2_R = pow(Q_R,2)/A_R + localG*I1_R;
		
	c_L = sqrt(g*A_L/b);
	c_R = sqrt(g*A_R/b);
		
	u_L = Q_L/A_L;
	u_R = Q_R/A_R;
	
	lambda1_L = u_L + sqrt(c_L*c_L - beta*u_L*u_L + beta*beta*u_L*u_L);
	lambda1_R = u_R + sqrt(c_R*c_R - beta*u_R*u_R + beta*beta*u_R*u_R);
	lambda2_L = u_L - sqrt(c_L*c_L - beta*u_L*u_L + beta*beta*u_L*u_L);
	lambda2_R = u_R - sqrt(c_R*c_R - beta*u_R*u_R + beta*beta*u_R*u_R);

	/*lambda1_L = u_L + c_L;
	lambda1_R = u_R + c_R;
	lambda2_L = u_L - c_L;
	lambda2_R = u_R - c_R;
	*/
	
	uhat = (Q_L/sqrt(A_L) + Q_R/sqrt(A_R))/(sqrt(A_L) + sqrt(A_R));
	chat = sqrt(g/2*(A_L + A_R)/b);
		
	/*lambdahat1 = uhat + chat;
	lambdahat2 = uhat - chat;
*/
	lambdahat1 = uhat + sqrt(chat*chat - beta*uhat*uhat + beta*beta*uhat*uhat);
	lambdahat2 = uhat - sqrt(chat*chat - beta*uhat*uhat + beta*beta*uhat*uhat);

	epsilon1 = max(0, (lambdahat1 - lambda1_L));
	epsilon1 = max(epsilon1, (lambda1_R - lambdahat1));
	lambdahat1 = max(epsilon1,fabs(lambdahat1));
		
	epsilon2 = max(0, (lambdahat2 - lambda2_L));
	epsilon2 = max(epsilon2, (lambda2_R - lambdahat2));
	lambdahat2 = max(epsilon2,fabs(lambdahat2));

	// Eigenvectors
	//double Rhat[2][2] = {{1,1},{uhat+chat,uhat-chat}};

	double Rhat[2][2];
	Rhat[0][0] = 1;
	Rhat[0][1] = 1;
	Rhat[1][0] = (beta*uhat*uhat - chat*chat)/(beta*uhat - sqrt(chat*chat - beta*uhat*uhat + beta*beta*uhat*uhat));
	Rhat[1][1] = (beta*uhat*uhat - chat*chat)/(beta*uhat + sqrt(chat*chat - beta*uhat*uhat + beta*beta*uhat*uhat));	

	double dtm = 1/(Rhat[0][0]*Rhat[1][1]-Rhat[0][1]*Rhat[1][0]);
	double Rhatinv[2][2] = {{dtm*Rhat[1][1], -dtm*Rhat[0][1]},{-dtm*Rhat[1][0], dtm*Rhat[0][0]}};

	double RoeJac[2][2];
	RoeJac[0][0] = lambdahat1*Rhat[0][0]*Rhatinv[0][0]+lambdahat2*Rhat[0][1]*Rhatinv[1][0]; 
	RoeJac[0][1] = lambdahat1*Rhat[0][0]*Rhatinv[0][1]+lambdahat2*Rhat[0][1]*Rhatinv[1][1]; 
	RoeJac[1][0] = lambdahat1*Rhat[1][0]*Rhatinv[0][0]+lambdahat2*Rhat[1][1]*Rhatinv[1][0]; 
	RoeJac[1][1] = lambdahat1*Rhat[1][0]*Rhatinv[0][1]+lambdahat2*Rhat[1][1]*Rhatinv[1][1]; 

	double jump[2] = {A_L - A_R, Q_L - Q_R};

	Fhat[0] = 0.5*(F1_L + F1_R) + 0.5*(RoeJac[0][0]*jump[0]+RoeJac[0][1]*jump[1]);
	Fhat[1] = 0.5*(F2_L + F2_R) + 0.5*(RoeJac[1][0]*jump[0]+RoeJac[1][1]*jump[1]);

}

/****************************************************************************************//**
*
* Function to calculate local Lax-Friedrich's flux for the 1-D St. Venant equations at an
*  interface (node).
* @param[in] A_L value of the wet cross-sectional area at the node coming from the 
* left element.
* @param [in] A_R value of the wet cross-sectional area at the node coming from the right
* element
* @param [in] Q_L value of the volumetric discharge at the node coming from the left element
* @param [in] Q_R value of the volumetric discharge at the node coming from the right element
* @param [in] b value of the width of the channel at the node. We assume that the width of the 
* channel is prescribed at each node and thus the value will be continuous at nodes
* @param [in] localG value of the gravitational acceleration constant. This will be the same
* as g, except when we set it to zero to handle wetting and drying.
* @param [out] Fhat a pointer to an array where the numerical flux computed will be stored
*
* ************************************************************************************/
void LF(double A_L, double A_R, double Q_L, double Q_R, double b, double localG, double *Fhat)
{
	double I1_L, I1_R;
	double F1_L, F1_R;
	double F2_L, F2_R;
	double c_L, c_R;
	double u_L, u_R;

	I1_L = getI1(A_L, b);
	I1_R = getI1(A_R, b);

	F1_L = Q_L;
	F1_R = Q_R;

	F2_L = pow(Q_L,2)/A_L + localG*I1_L;
	F2_R = pow(Q_R,2)/A_R + localG*I1_R;
		
	c_L = sqrt(g*A_L/b);
	c_R = sqrt(g*A_R/b);
		
	u_L = Q_L/A_L;
	u_R = Q_R/A_R;

	double C = max(fabs(u_L + c_L), fabs(u_R+c_R));
	C = max(C, fabs(u_L - c_L));
	C = max(C, fabs(u_R - c_R));
	double jump[2] = {A_L - A_R, Q_L - Q_R};
	Fhat[0] = 0.5*(F1_L+F1_R) + 0.5*C*jump[0];
	Fhat[1] = 0.5*(F2_L+F2_R) + 0.5*C*jump[1];
	

}


