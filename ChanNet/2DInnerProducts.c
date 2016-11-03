#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mathfunctions.h"
#include "MeshAttributes.h"

/**********************************************************************//**
* @file 2DInnerProducts.c
* *************************************************************************/

/* @cond GLOBAL_VARS */
extern const double g;
/* @endcond */

/************************************************************************************//**
* Function for evaluating the inner products that occur on the right hand side of the
* DG equations for the 2-D junctions
* @param [in] junc a pointer to the junction structure corresponding to the junction that
*  is currently being worked on.
* @param [in] Fhat1dotn a pointer to the array that contains the first component of the
* numerical normal flux vector at the edges. Size = NumEl x 3
* @param [in] Fhat2dotn a pointer to the array that contains the second component of the
* numerical normal flux vector at the edges. Size = NumEl x 3
* @param [in] Fhat3dotn a pointer to the array that contains the third component of the 
* numerical normal flux vector at the edges. Size = NumEl x 3
* @param [in] el integer element number of the element that is currently being worked on
* @param [out] SI pointer to a double array of size 9 in which all the surface integrals
* are returned
* @param [out] VI pointer to a double array of size 21 in which all the area integrals
* are returned
* ***************************************************************************************/ 

void Compute2DInnerProducts(struct junction *junc, double *Fhat1dotn, double *Fhat2dotn, 
							double *Fhat3dotn, int el, double *SI, double *VI)
{

	int v0 = junc->EltoVert[el*3];
	int v1 = junc->EltoVert[el*3+1];
	int v2 = junc->EltoVert[el*3+2];
	
	double xv1, xv2, xv3, yv1, yv2, yv3;
	double x1, x2, x3, y1, y2, y3;
	double edg_length[3], edg_jac[3];

	xv1 = junc->x[v0];
	xv2 = junc->x[v1];
	xv3 = junc->x[v2];
	yv1 = junc->y[v0];
	yv2 = junc->y[v1];
	yv3 = junc->y[v2];

	// midpoints of edges
	x1 = 0.5*(xv1+xv2);
	x2 = 0.5*(xv2+xv3);
	x3 = 0.5*(xv3+xv1);
	y1 = 0.5*(yv1+yv2);
	y2 = 0.5*(yv2+yv3);
	y3 = 0.5*(yv3+yv1);

	edg_length[0] = sqrt((xv2-xv1)*(xv2-xv1)+(yv2-yv1)*(yv2-yv1));
	edg_length[1] = sqrt((xv3-xv2)*(xv3-xv2)+(yv3-yv2)*(yv3-yv2));
	edg_length[2] = sqrt((xv3-xv1)*(xv3-xv1)+(yv3-yv1)*(yv3-yv1));
	
	// Store the jacobian for each edge
	for (int e = 0; e < 3; e++)
	{
		edg_jac[e] = edg_length[e]/2;
	}

	// quadrature points on the edges
	double x_e1[2], x_e2[2], x_e3[2];
	double y_e1[2], y_e2[2], y_e3[2];

	double sqrt3 = sqrt(3);
	x_e1[0] = (sqrt3*x1-x2+x3)/sqrt3;
	y_e1[0] = (sqrt3*y1-y2+y3)/sqrt3;
	x_e1[1] = (sqrt3*x1+x2-x3)/sqrt3;
	y_e1[1] = (sqrt3*y1+y2-y3)/sqrt3;

	x_e2[0] = (x1+sqrt3*x2-x3)/sqrt3;
	y_e2[0] = (y1+sqrt3*y2-y3)/sqrt3;
	x_e2[1] = (-x1+sqrt3*x2+x3)/sqrt3;
	y_e2[1] = (-y1+sqrt3*y2+y3)/sqrt3;

	x_e3[0] = (-x1+x2+sqrt3*x3)/sqrt3;
	y_e3[0] = (-y1+y2+sqrt3*y3)/sqrt3;
	x_e3[1] = (x1-x2+sqrt3*x3)/sqrt3;
	y_e3[1] = (y1-y2+sqrt3*y3)/sqrt3;


	// Value of the basis functions at the quadrature points on the edges
	double Psi1_e1[2], Psi1_e2[2], Psi1_e3[2];
	double Psi2_e1[2], Psi2_e2[2], Psi2_e3[2];
	double Psi3_e1[2], Psi3_e2[2], Psi3_e3[2];

	double denom = x1*y2 - x1*y3 + x2*y3 - x2*y1 + x3*y1 - x3*y2;
	for (int n=0; n<2; ++n)
	{
		Psi1_e1[n] = -(y_e1[n]*x2 - y_e1[n]*x3 - x_e1[n]*y2 + x3*y2 + x_e1[n]*y3 - x2*y3)/denom;
		Psi1_e2[n] = -(y_e2[n]*x2 - y_e2[n]*x3 - x_e2[n]*y2 + x3*y2 + x_e2[n]*y3 - x2*y3)/denom;
		Psi1_e3[n] = -(y_e3[n]*x2 - y_e3[n]*x3 - x_e3[n]*y2 + x3*y2 + x_e3[n]*y3 - x2*y3)/denom;
		
		Psi2_e1[n] = -(-x1*y_e1[n] + x3*y_e1[n] -x3*y1 + x_e1[n]*y1 -x_e1[n]*y3 + x1*y3)/denom;
		Psi2_e2[n] = -(-x1*y_e2[n] + x3*y_e2[n] -x3*y1 + x_e2[n]*y1 -x_e2[n]*y3 + x1*y3)/denom;
		Psi2_e3[n] = -(-x1*y_e3[n] + x3*y_e3[n] -x3*y1 + x_e3[n]*y1 -x_e3[n]*y3 + x1*y3)/denom;

		Psi3_e1[n] = -(y_e1[n]*x1 - y_e1[n]*x2 - x_e1[n]*y1 + x2*y1 + x_e1[n]*y2 - x1*y2)/denom;
		Psi3_e2[n] = -(y_e2[n]*x1 - y_e2[n]*x2 - x_e2[n]*y1 + x2*y1 + x_e2[n]*y2 - x1*y2)/denom;
		Psi3_e3[n] = -(y_e3[n]*x1 - y_e3[n]*x2 - x_e3[n]*y1 + x2*y1 + x_e3[n]*y2 - x1*y2)/denom;
	}

	// The value of Fhat at the Gauss nodes (which are the midpoints). For the first edge, our interpolation node is not the midpoint. So we need to compute the value at the midpoint. f(1/2,0) = f1 + f2(1-2*r) + f3*(-1+2*r);
	double F1_edg0, F1_edg1, F1_edg2;
	double F2_edg0, F2_edg1, F2_edg2;
	double F3_edg0, F3_edg1, F3_edg2;


	F1_edg1 = Fhat1dotn[el*3+1];
	F1_edg2 = Fhat1dotn[el*3+2];
	F1_edg0 = Fhat1dotn[el*3];

	F2_edg1 = Fhat2dotn[el*3+1];
	F2_edg2 = Fhat2dotn[el*3+2];
	F2_edg0 = Fhat2dotn[el*3]; 

	F3_edg1 = Fhat3dotn[el*3+1];
	F3_edg2 = Fhat3dotn[el*3+2];
	F3_edg0 = Fhat3dotn[el*3]; 

	// Solve for a, b and c in F = [a+bx; a+cy] by solving F.n = Fhatdotn

	// Get the normal vectors
	double jac = 4*(x2-x3)*(y2-y1)-4*(x2-x1)*(y2-y3);
	double nx1 = sign(jac)*2*(y2-y3)/edg_length[0];
	double ny1 = sign(jac)*2*(x3-x2)/edg_length[0];
	double nx2 = sign(jac)*2*(y3-y1)/edg_length[1];
	double ny2 = sign(jac)*2*(x1-x3)/edg_length[1];
	double nx3 = sign(jac)*2*(y1-y2)/edg_length[2];
	double ny3 = sign(jac)*2*(x2-x1)/edg_length[2];


	double M11 = nx1; double M12 = x1*nx1+y1*ny1; double M13 = ny1;
	double M21 = nx2; double M22 = x2*nx2+y2*ny2; double M23 = ny2;
	double M31 = nx3; double M32 = x3*nx3+y3*ny3; double M33 = ny3;

	double R11 = F1_edg0; double R12= F1_edg1; double R13 = F1_edg2; 
	double R21 = F2_edg0; double R22= F2_edg1; double R23 = F2_edg2; 
	double R31 = F3_edg0; double R32= F3_edg1; double R33 = F3_edg2; 

	denom = M13*M22*M31 - M12*M23*M31 - M13*M21*M32 + M11*M23*M32 + M12*M21*M33 - M11*M22*M33;

	double a1 = (R13*M13*M22 - R13*M12*M23 - R12*M13*M32 + R11*M23*M32 + R12*M12*M33 - R11*M22*M33)/denom;
	double b1 = -(R13*M13*M21 - R13*M11*M23 - R12*M13*M31 + R11*M23*M31 + R12*M11*M33 - R11*M21*M33)/denom;
	double c1 = (R13*M12*M21 - R13*M11*M22 - R12*M12*M31 + R11*M22*M31 + R12*M11*M32 - R11*M21*M32)/denom;

	double a2 = (R23*M13*M22 - R23*M12*M23 - R22*M13*M32 + R21*M23*M32 + R22*M12*M33 - R21*M22*M33)/denom;
	double b2 = -(R23*M13*M21 - R23*M11*M23 - R22*M13*M31 + R21*M23*M31 + R22*M11*M33 - R21*M21*M33)/denom;
	double c2 = (R23*M12*M21 - R23*M11*M22 - R22*M12*M31 + R21*M22*M31 + R22*M11*M32 - R21*M21*M32)/denom;

	double a3 = (R33*M13*M22 - R33*M12*M23 - R32*M13*M32 + R31*M23*M32 + R32*M12*M33 - R31*M22*M33)/denom;
	double b3 = -(R33*M13*M21 - R33*M11*M23 - R32*M13*M31 + R31*M23*M31 + R32*M11*M33 - R31*M21*M33)/denom;
	double c3 = (R33*M12*M21 - R33*M11*M22 - R32*M12*M31 + R31*M22*M31 + R32*M11*M32 - R31*M21*M32)/denom;


	// The value of Fdotn at the interpolation points
	double F1dotn_e1[2], F1dotn_e2[2], F1dotn_e3[2];
	double F2dotn_e1[2], F2dotn_e2[2], F2dotn_e3[2];
	double F3dotn_e1[2], F3dotn_e2[2], F3dotn_e3[2];

	for (int n = 0; n < 2; ++n)
	{
		F1dotn_e1[n] = (a1+b1*x_e1[n])*nx1 + (c1+b1*y_e1[n])*ny1; 
		F1dotn_e2[n] = (a1+b1*x_e2[n])*nx2 + (c1+b1*y_e2[n])*ny2; 
		F1dotn_e3[n] = (a1+b1*x_e3[n])*nx3 + (c1+b1*y_e3[n])*ny3;	

		F2dotn_e1[n] = (a2+b2*x_e1[n])*nx1 + (c2+b2*y_e1[n])*ny1; 
		F2dotn_e2[n] = (a2+b2*x_e2[n])*nx2 + (c2+b2*y_e2[n])*ny2; 
		F2dotn_e3[n] = (a2+b2*x_e3[n])*nx3 + (c2+b2*y_e3[n])*ny3;	

		F3dotn_e1[n] = (a3+b3*x_e1[n])*nx1 + (c3+b3*y_e1[n])*ny1; 
		F3dotn_e2[n] = (a3+b3*x_e2[n])*nx2 + (c3+b3*y_e2[n])*ny2; 
		F3dotn_e3[n] = (a3+b3*x_e3[n])*nx3 + (c3+b3*y_e3[n])*ny3; 
	}


	SI[0] = edg_jac[0]*(F1dotn_e1[0]*Psi1_e1[0] + F1dotn_e1[1]*Psi1_e1[1]) + 
	        edg_jac[1]*(F1dotn_e2[0]*Psi1_e2[0] + F1dotn_e2[1]*Psi1_e2[1]) + 
	        edg_jac[2]*(F1dotn_e3[0]*Psi1_e3[0] + F1dotn_e3[1]*Psi1_e3[1]);

	SI[1] = edg_jac[0]*(F1dotn_e1[0]*Psi2_e1[0] + F1dotn_e1[1]*Psi2_e1[1]) + 
	        edg_jac[1]*(F1dotn_e2[0]*Psi2_e2[0] + F1dotn_e2[1]*Psi2_e2[1]) + 
	        edg_jac[2]*(F1dotn_e3[0]*Psi2_e3[0] + F1dotn_e3[1]*Psi2_e3[1]);

	SI[2] = edg_jac[0]*(F1dotn_e1[0]*Psi3_e1[0] + F1dotn_e1[1]*Psi3_e1[1]) + 
	        edg_jac[1]*(F1dotn_e2[0]*Psi3_e2[0] + F1dotn_e2[1]*Psi3_e2[1]) + 
	        edg_jac[2]*(F1dotn_e3[0]*Psi3_e3[0] + F1dotn_e3[1]*Psi3_e3[1]);

	SI[3] = edg_jac[0]*(F2dotn_e1[0]*Psi1_e1[0] + F2dotn_e1[1]*Psi1_e1[1]) + 
	        edg_jac[1]*(F2dotn_e2[0]*Psi1_e2[0] + F2dotn_e2[1]*Psi1_e2[1]) + 
	        edg_jac[2]*(F2dotn_e3[0]*Psi1_e3[0] + F2dotn_e3[1]*Psi1_e3[1]);

	SI[4] = edg_jac[0]*(F2dotn_e1[0]*Psi2_e1[0] + F2dotn_e1[1]*Psi2_e1[1]) + 
	        edg_jac[1]*(F2dotn_e2[0]*Psi2_e2[0] + F2dotn_e2[1]*Psi2_e2[1]) + 
	        edg_jac[2]*(F2dotn_e3[0]*Psi2_e3[0] + F2dotn_e3[1]*Psi2_e3[1]);

	SI[5] = edg_jac[0]*(F2dotn_e1[0]*Psi3_e1[0] + F2dotn_e1[1]*Psi3_e1[1]) + 
	        edg_jac[1]*(F2dotn_e2[0]*Psi3_e2[0] + F2dotn_e2[1]*Psi3_e2[1]) + 
	        edg_jac[2]*(F2dotn_e3[0]*Psi3_e3[0] + F2dotn_e3[1]*Psi3_e3[1]);

	SI[6] = edg_jac[0]*(F3dotn_e1[0]*Psi1_e1[0] + F3dotn_e1[1]*Psi1_e1[1]) + 
	        edg_jac[1]*(F3dotn_e2[0]*Psi1_e2[0] + F3dotn_e2[1]*Psi1_e2[1]) + 
	        edg_jac[2]*(F3dotn_e3[0]*Psi1_e3[0] + F3dotn_e3[1]*Psi1_e3[1]);

	SI[7] = edg_jac[0]*(F3dotn_e1[0]*Psi2_e1[0] + F3dotn_e1[1]*Psi2_e1[1]) + 
	        edg_jac[1]*(F3dotn_e2[0]*Psi2_e2[0] + F3dotn_e2[1]*Psi2_e2[1]) + 
	        edg_jac[2]*(F3dotn_e3[0]*Psi2_e3[0] + F3dotn_e3[1]*Psi2_e3[1]);

	SI[8] = edg_jac[0]*(F3dotn_e1[0]*Psi3_e1[0] + F3dotn_e1[1]*Psi3_e1[1]) + 
	        edg_jac[1]*(F3dotn_e2[0]*Psi3_e2[0] + F3dotn_e2[1]*Psi3_e2[1]) + 
		   edg_jac[2]*(F3dotn_e3[0]*Psi3_e3[0] + F3dotn_e3[1]*Psi3_e3[1]);


	/* Area Integrals */

	// weights associated with the Gauss nodes
	double weight = 1.0/6;

	
	// Store the value of F at the three Gauss nodes for each element
	// the three Gauss nodes are the midpoints of the edges, two of which are the points at which we know the values of our variables.

	double zeta_rearr[3];
	double Qx_rearr[3];
	double Qy_rearr[3];
	double H[3];
	double tau[3];
	double gradz[2];
	double gradPsi1[2], gradPsi2[2], gradPsi3[2];
	double midz[3];
	double nFriction[3];

	int n0 = junc->EltoVert[el*3];
	int n1 = junc->EltoVert[el*3+1];
	int n2 = junc->EltoVert[el*3+2];

	midz[0] = 0.5*(junc->z[n0]+junc->z[n1]);
	midz[1] = 0.5*(junc->z[n1]+junc->z[n2]);
	midz[2] = 0.5*(junc->z[n2]+junc->z[n0]);

	nFriction[0] = 0.5*(junc->nFriction[n0]+junc->nFriction[n1]);
	nFriction[1] = 0.5*(junc->nFriction[n1]+junc->nFriction[n2]);
	nFriction[2] = 0.5*(junc->nFriction[n2]+junc->nFriction[n0]);

	// Calculate the cartesian derivatives of the basis functions
	gradPsi1[0] = 4/jac*(y2-y3);
	gradPsi1[1] = 4/jac*(x3-x2);
	gradPsi2[0] = 4/jac*(y3-y1);
	gradPsi2[1] = 4/jac*(x1-x3);
	gradPsi3[0] = 4/jac*(y1-y2);
	gradPsi3[1] = 4/jac*(x2-x1);

	for (int e=0; e < 3; ++e)
	{
		zeta_rearr[e] = junc->zeta[el*3+e];
		Qx_rearr[e] = junc->Qx[el*3+e];
		Qy_rearr[e] = junc->Qy[el*3+e];
		H[e] = zeta_rearr[e]+midz[e];
		double u = Qx_rearr[e]/H[e];
		double v = Qy_rearr[e]/H[e];
		tau[e] = g*nFriction[e]*nFriction[e]*sqrt(u*u+v*v)/(pow(H[e],4./3));
	}

	double F1x[3], F2x[3], F3x[3];
	double F1y[3], F2y[3], F3y[3];

	double localG = g;

	#ifdef WDON
	if (junc->WD[el] == 0)
		localG = 0;
	#endif

	for (int e=0; e<3; ++e)
	{
		F1x[e] = Qx_rearr[e];
		F2x[e] = pow(Qx_rearr[e],2)/H[e]+0.5*localG*(pow(H[e],2)-pow(midz[e],2));
		F3x[e] = Qx_rearr[e]*Qy_rearr[e]/H[e];
		
		F1y[e] = Qy_rearr[e];
		F2y[e] = Qx_rearr[e]*Qy_rearr[e]/H[e];	
		F3y[e] = pow(Qy_rearr[e],2)/H[e]+0.5*localG*(pow(H[e],2)-pow(midz[e],2));
	}

	// gradz constant over an element
	double z1 = midz[0];
	double z2 = midz[1];
	double z3 = midz[2];
	
	for (int j=0; j<2; ++j)
	{
		gradz[j] = z1*gradPsi1[j] + z2*gradPsi2[j] + z3*gradPsi3[j]; 				
	
	}
	
	
	 double VF1dotn_e1 = F1x[0]*nx1 + F1y[0]*ny1;
	 double VF1dotn_e2 = F1x[1]*nx2 + F1y[1]*ny2;
	 double VF1dotn_e3 = F1x[2]*nx3 + F1y[2]*ny3;
	
	 double VF2dotn_e1 = F2x[0]*nx1 + F2y[0]*ny1;
	 double VF2dotn_e2 = F2x[1]*nx2 + F2y[1]*ny2;
	 double VF2dotn_e3 = F2x[2]*nx3 + F2y[2]*ny3;
	
	 double VF3dotn_e1 = F3x[0]*nx1 + F3y[0]*ny1;
	 double VF3dotn_e2 = F3x[1]*nx2 + F3y[1]*ny2;
	 double VF3dotn_e3 = F3x[2]*nx3 + F3y[2]*ny3;
	
	 M11 = nx1;  M12 = x1*nx1+y1*ny1;  M13 = ny1;
	 M21 = nx2;  M22 = x2*nx2+y2*ny2;  M23 = ny2;
	 M31 = nx3;  M32 = x3*nx3+y3*ny3;  M33 = ny3;
	
	 R11 = VF1dotn_e1;  R12 = VF1dotn_e2;  R13 = VF1dotn_e3; 
	 R21 = VF2dotn_e1;  R22 = VF2dotn_e2;  R23 = VF2dotn_e3; 
	 R31 = VF3dotn_e1;  R32 = VF3dotn_e2;  R33 = VF3dotn_e3; 
	
	 denom = M13*M22*M31 - M12*M23*M31 - M13*M21*M32 + M11*M23*M32 + M12*M21*M33 - M11*M22*M33;
	
	 a1 = (R13*M13*M22 - R13*M12*M23 - R12*M13*M32 + R11*M23*M32 + R12*M12*M33 - R11*M22*M33)/denom;
	 b1 = -(R13*M13*M21 - R13*M11*M23 - R12*M13*M31 + R11*M23*M31 + R12*M11*M33 - R11*M21*M33)/denom;
	 c1 = (R13*M12*M21 - R13*M11*M22 - R12*M12*M31 + R11*M22*M31 + R12*M11*M32 - R11*M21*M32)/denom;
	
	 a2 = (R23*M13*M22 - R23*M12*M23 - R22*M13*M32 + R21*M23*M32 + R22*M12*M33 - R21*M22*M33)/denom;
	 b2 = -(R23*M13*M21 - R23*M11*M23 - R22*M13*M31 + R21*M23*M31 + R22*M11*M33 - R21*M21*M33)/denom;
	 c2 = (R23*M12*M21 - R23*M11*M22 - R22*M12*M31 + R21*M22*M31 + R22*M11*M32 - R21*M21*M32)/denom;
	
	 a3 = (R33*M13*M22 - R33*M12*M23 - R32*M13*M32 + R31*M23*M32 + R32*M12*M33 - R31*M22*M33)/denom;
	 b3 = -(R33*M13*M21 - R33*M11*M23 - R32*M13*M31 + R31*M23*M31 + R32*M11*M33 - R31*M21*M33)/denom;
	 c3 = (R33*M12*M21 - R33*M11*M22 - R32*M12*M31 + R31*M22*M31 + R32*M11*M32 - R31*M21*M32)/denom;


	VI[0] = fabs(jac)*(0.5*a1*gradPsi1[0]+0.5*c1*gradPsi1[1] + 
	     	   b1*weight*(gradPsi1[0]*(x1+ x2 + x3) + 
	     	   gradPsi1[1]*(y1+y2+y3)));
	VI[1] = fabs(jac)*(0.5*a1*gradPsi2[0]+0.5*c1*gradPsi2[1] + 
	     	   b1*weight*(gradPsi2[0]*(x1+ x2 + x3) + 
	     	   gradPsi2[1]*(y1+y2+y3)));
	VI[2] = fabs(jac)*(0.5*a1*gradPsi3[0]+0.5*c1*gradPsi3[1] + 
	     	   b1*weight*(gradPsi3[0]*(x1+ x2 + x3) + 
	     	   gradPsi3[1]*(y1+y2+y3)));
	
	VI[3] = fabs(jac)*(0.5*a2*gradPsi1[0]+0.5*c2*gradPsi1[1] + 
	     	   b2*weight*(gradPsi1[0]*(x1+ x2 + x3) + 
	     	   gradPsi1[1]*(y1+y2+y3)));
	VI[4] = localG*fabs(jac)*gradz[0]*weight*zeta_rearr[0];
	VI[5] = fabs(jac)*weight*tau[0]*Qx_rearr[0];
	
	VI[6] = fabs(jac)*(0.5*a2*gradPsi2[0]+0.5*c2*gradPsi2[1] + 
	     	   b2*weight*(gradPsi2[0]*(x1+ x2 + x3) + 
	     	   gradPsi2[1]*(y1+y2+y3)));
	VI[7] = localG*fabs(jac)*gradz[0]*weight*zeta_rearr[1];
	VI[8] = fabs(jac)*weight*tau[1]*Qx_rearr[1];
	
	VI[9] = fabs(jac)*(0.5*a2*gradPsi3[0]+0.5*c2*gradPsi3[1] + 
	     	   b2*weight*(gradPsi3[0]*(x1+ x2 + x3) + 
	     	   gradPsi3[1]*(y1+y2+y3)));
	VI[10] = localG*fabs(jac)*gradz[0]*weight*zeta_rearr[2];
	VI[11] = fabs(jac)*weight*tau[2]*Qx_rearr[2];
	
	VI[12] = fabs(jac)*(0.5*a3*gradPsi1[0]+0.5*c3*gradPsi1[1] + 
	     	   b3*weight*(gradPsi1[0]*(x1+ x2 + x3) + 
	     	   gradPsi1[1]*(y1+y2+y3)));
	VI[13] = localG*fabs(jac)*gradz[1]*weight*zeta_rearr[0];
	VI[14] = fabs(jac)*weight*tau[0]*Qy_rearr[0];
	
	VI[15] = fabs(jac)*(0.5*a3*gradPsi2[0]+0.5*c3*gradPsi2[1] + 
	     	   b3*weight*(gradPsi2[0]*(x1+ x2 + x3) + 
	     	   gradPsi2[1]*(y1+y2+y3)));
	VI[16] = localG*fabs(jac)*gradz[1]*weight*zeta_rearr[1];
	VI[17] = fabs(jac)*weight*tau[1]*Qy_rearr[1];
	
	VI[18] = fabs(jac)*(0.5*a3*gradPsi3[0]+0.5*c3*gradPsi3[1] + 
	     	   b3*weight*(gradPsi3[0]*(x1+ x2 + x3) + 
	     	   gradPsi3[1]*(y1+y2+y3)));
	VI[19] = localG*fabs(jac)*gradz[1]*weight*zeta_rearr[2];
	VI[20] = fabs(jac)*weight*tau[2]*Qy_rearr[2];



}


