
/****************** Boundary Conditions ************************************/	

#include <math.h>
#include "globals.h"

extern const double g;

void boundary_conditions()
{
	// propagation of a shock
//	Q[0] = 140;
//	A[0] = 12;


	// smooth subcritical flow
/*	Q[0] = 4.42;
	A[0] = A[1];
	double h = 2.0;
	A[NumNodes-1] = b[NumEdges-1]*h;
	Q[NumNodes-1] = Q[NumNodes-2];
*/

// smooth transcritical flow
Q[0] = 1.53;
A[0] = A[1];
double h = 0.66;
A[NumNodes-1] = b[NumEdges-1]*h;
Q[NumNodes-1] = Q[NumNodes-2];


/*	// transcritical with shock
	Q[0] = 0.18;
	double h = 0.33;
	A[0] = A[1];
	Q[NumNodes-1] = Q[NumNodes-2];
	A[NumNodes-1] = b[NumEdges-1]*h;	
*/	
	// Dam break boundary conditions
/*	A[0] = 0.005;
	Q[0] = Q[1];
	A[2*NumNodes-1] = A[2*NumNodes-2];
	Q[2*NumNodes-1] = Q[2*NumNodes-2];
*/
	// No flow in or out
//	Q[0] = 0;
//	Q[2*NumNodes-1] = 0;
//	//A[0] = 0.1;
//	A[0] = A[1];
//	A[2*NumNodes-1] = A[2*NumNodes-2];


	// Parabolic bowl
/*	double Hdry = 1e-3;
	A[0] = Hdry*b[0];
	Q[0] = 0;
	A[2*NumNodes-1] = Hdry*b[NumNodes-1];
	Q[2*NumNodes-1] = 0;
*/

	// Macdonalds short channel with smooth transition and shock
	//A[0] = A[1];
	//Q[0] = 2;
	//double hDownstream = 2.87844;
	//double bDownstream = b[NumNodes-1];
	//A[2*NumNodes-1] = bDownstream*hDownstream;
	//Q[2*NumNodes-1] = Q[2*NumNodes-2];
	
/*	// Crossley_Wright converging/diverging channel (S1)
	Q[0] = 20;
	A[0] = A[1];
	double h = 0.1;
	double height = A[NumNodes-1]/b[NumEdges-1];
	double H_w = height - h;
	double u = 2.0/3*0.6*pow(2*g*H_w,0.5);
	Q[NumNodes-1] = u*A[NumNodes-1];
	A[NumNodes-1] = A[NumNodes-2];

*/

	// Crossley_Wright trapezoidal channel (S2)  (MacDonald's test case)
/*	Q[0] = 20;
	A[0] = A[1];
	double height = 1.3449963;
	double m1val = m1[NumEdges-1];
	double m2val = m2[NumEdges-1];
	A[NumNodes-1] = height*b[NumEdges-1] + 0.5*m1val*height*height + 0.5*m2val*height*height;
	Q[NumNodes-1] = Q[NumNodes-2];
*/

	
}




