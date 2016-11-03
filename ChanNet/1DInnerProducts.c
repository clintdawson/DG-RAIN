#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"
#include "ChannelsAndJunctions.h"

/**************************************************************************//**
* @file 1DInnerProducts.c 
*****************************************************************************/

/* @cond FUNCTION_PROTOTYPES */
extern const double g;
extern double getI1(double A, double b);
extern double getI2(double A, double b, double db);
extern double getS_f(double A, double Q, double b, double n);
/* @endcond */

/********************************************************************//**
* Function for evaluating the inner products that occur on the right hand side
* of the DG equations for the 1-D channels
* @param [in] Chan the channel structure corresponding to the channel that is 
* currently being worked on
* @param [in] el integer element number of the channel element that is currently 
* being worked on
* @param [out] IP Pointer to a double array of size 10 in which the ten inner 
* products are returned
*************************************************************************/


void Compute1DInnerProducts(struct channel *Chan, int el, double *IP)
{

	// Gauss nodes for a 2-point Gauss quadrature rule
	double t[] = {-0.577350269189626,  0.577350269189626};
	//Jacobian of the transformation 
	double x1 = Chan->x[el];
	double x2 = Chan->x[el+1];
	double y1 = Chan->y[el];
	double y2 = Chan->y[el+1];
	double h = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	double jac = 0.5*h;
	
	// Get the value of the basis functions and the derivatives of the basis functions at the Gauss points for all elements.
	int sizeoft =  sizeof(t)/sizeof(t[0]);
	double psi1[sizeoft], psi2[sizeoft], dpsi1, dpsi2;
	
	for (int i =0; i < sizeoft; ++i){
		psi1[i] = 0.5*(1 - t[i]);
		psi2[i] = 0.5*(1 + t[i]);
	}

	// Gradient along the direction of the channel
	dpsi1 = -0.5/jac;
	dpsi2 = 0.5/jac;

	// bathymetry and width data for the channel
	double z1 = Chan->z[el];
	double z2 = Chan->z[el+1];
	double b1 = Chan->b[el];
	double b2 = Chan->b[el+1];

	//double beta1 = Chan->beta[el];
	//double beta2 = Chan->beta[el+1];	

	double A1, A2, Q1, Q2, b_t1, b_t2, db,nFric1,nFric2;
	double I1_1, I1_2, I2_1, I2_2, S_0;
	double S_f1, S_f2;
	double betaquad1, betaquad2;

	// Get the value of Area and Discharge at the Gauss node
	// Also get the breadth and its derivatives at the Gauss points 
	A1 = psi1[0]*Chan->A[2*el+1] + psi2[0]*Chan->A[2*el+2];	
	A2 = psi1[1]*Chan->A[2*el+1] + psi2[1]*Chan->A[2*el+2];

	// Linear interpolation of the Q values at the nodes
	Q1 = psi1[0]*Chan->Q[2*el+1] + psi2[0]*Chan->Q[2*el+2];
	Q2 = psi1[1]*Chan->Q[2*el+1] + psi2[1]*Chan->Q[2*el+2];
	
	// Linear interpolatio of the breadth values at the nodes
	b_t1 = psi1[0]*b1 + psi2[0]*b2;
	b_t2 = psi1[1]*b1 + psi2[1]*b2;
	
	//betaquad1 = psi1[0]*beta1 + psi2[0]*beta2;
	//betaquad2 = psi1[1]*beta1 + psi2[1]*beta2;

	betaquad1 = 1;
	betaquad2 = 1;

	nFric1 = psi1[0]*Chan->nFriction[el] + psi2[0]*Chan->nFriction[el+1];
	nFric2 = psi1[1]*Chan->nFriction[el] + psi2[1]*Chan->nFriction[el+1];
	
	// value of I1 at the Gauss points
	I1_1 = getI1(A1, b_t1);
	I1_2 = getI1(A2, b_t2);

	// value of S_f at the Gauss points
	S_f1 = getS_f(A1, Q1, b_t1, nFric1);
	S_f2 = getS_f(A2, Q2, b_t2, nFric2);

	db = (b2-b1)/h;
	
	// S_0 will be piecewise constant over the element
	S_0 = (z2-z1)/h;

	I2_1 = getI2(A1, b_t1,db);
	I2_2 = getI2(A2, b_t2,db);
	
	// Compute the inner products we are interested in 
	
	double localG;
	#ifdef WDON
	if (Chan->WD[el] == 1)
		localG = g;
	else
		localG = 0;
	#else
		localG = g;
	#endif


	IP[0] = jac*dpsi1*(Q1+Q2);
	IP[1] = jac*dpsi2*(Q1+Q2);

	IP[2] = jac*dpsi1*(betaquad1*Q1*Q1/A1+localG*I1_1 + Q2*Q2*betaquad2/A2 + localG*I1_2);
	IP[3] = localG*jac*(psi1[0]*I2_1 + psi1[1]*I2_2);
	IP[4] = localG*jac*S_0*(psi1[0]*A1+psi1[1]*A2);
	IP[5] = localG*jac*(psi1[0]*A1*S_f1 + psi1[1]*A2*S_f2);		

	IP[6] = jac*dpsi2*(betaquad1*Q1*Q1/A1+localG*I1_1 + Q2*Q2*betaquad2/A2 + localG*I1_2);
	IP[7] = localG*jac*(psi2[0]*I2_1 + psi2[1]*I2_2);
	IP[8] = localG*jac*S_0*(psi2[0]*A1+psi2[1]*A2);
	IP[9] = localG*jac*(psi2[0]*A1*S_f1 + psi2[1]*A2*S_f2);		

}


