#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathfunctions.h"
#include "MeshAttributes.h"

/*****************************************************************************************//**
* @file compute2DL.c
*
* This file contains code to evaluate the right hand side of the following discrete ODE obtained
* from the DG discretization of the 2-D shallow water equations:
* 
* \f$\frac{\partial \hat{\mathbf{w}}}{\partial t} = M^{-1}L(\hat{\mathbf{w}},t)\f$
*
* *******************************************************************************************/

/* @cond FUNCTION_PROTOTYPES */

extern double RoeFlux2D(double zeta_in, double zeta_ex, double Qx_in, double Qx_ex, 
	double Qy_in, double Qy_ex, double z_edge, double nx, double ny, double localG, double *Fhatdotn);
extern void Compute2DInnerProducts(struct junction *junc, double *Fhat1dotn, double *Fhat2dotn, double *Fhat3dotn, int el, double *SI, double *VI);

extern const double g;

/* @endcond */

/***********************************************************************************************//**
* Function for evaluating the right hand side of the discrete ODE obtained from the DG discretization
* of the 2-D shallow water equations
* @param [in] junc a pointer to the junction structure corresponding to the junction that is 
* currently being worked on 
* @param [in] time  a double representing the current time of the simulation
* @param [out] RHSZeta a pointer to an array of size NumEl x 3 in which the right hand side for
* the water surface elevation will be stored
* @param [out] RHSQx a pointer to an array of size NumEl x 3 in which the right hand side for the
* momentum in the x-direction is stored
* @param [out] RHSQy a pointer to an array of size NumEl x 3 in which the right hand side for the
* momentum in the y-direction is stored
*
***************************************************************************************************/

void compute2DL(struct junction *junc, double time, double *RHSZeta, double *RHSQx, double *RHSQy)
{

/********************************* Compute Roe's flux at the face **********************************************/
	int NumEdges = junc->TotalNumEdges;
	int NumEl = junc->NumEl;
	
	double Fhat1dotn[NumEl*3], Fhat2dotn[NumEl*3], Fhat3dotn[NumEl*3];

	for (int i=0; i<NumEdges; ++i)
	{
		double zeta_in, Qx_in, Qy_in;
		double zeta_ex, Qx_ex, Qy_ex;
		double z_edge;
		double nx, ny;
		double tx, ty;	
	
		double Fn_in[3], Fn_ex[3];
		
		int el_in = junc->EdgtoEls[index(i,0,2)];
		int el_ex = junc->EdgtoEls[index(i,1,2)];

		int int_el_local_edg;
		int ex_el_local_edg;
	
		int v0 = junc->EdgtoVert[i*2];
		int v1 = junc->EdgtoVert[i*2+1];
	
		int el_in_v0 = junc->EltoVert[el_in*3];
		int el_in_v1 = junc->EltoVert[el_in*3+1];
		int el_in_v2 = junc->EltoVert[el_in*3+2];
	
		int el_ex_v0 = junc->EltoVert[el_ex*3];
		int el_ex_v1 = junc->EltoVert[el_ex*3+1];
		int el_ex_v2 = junc->EltoVert[el_ex*3+2];
	
		if ((v0 == el_in_v0 && v1 == el_in_v1) || (v0 == el_in_v1 && v1 == el_in_v0))
			int_el_local_edg = 0;
		else if ((v0 == el_in_v1 && v1 == el_in_v2) || (v0 == el_in_v2 && v1 == el_in_v1))
			int_el_local_edg = 1;
		else if ((v0 == el_in_v2 && v1 == el_in_v0) || (v0 == el_in_v0 && v1 == el_in_v2))
			int_el_local_edg = 2;
		else
		{
			printf("error in calculating interior element local edge number \n");
			exit(EXIT_FAILURE);
		}
	
		if ((v0 == el_ex_v0 && v1 == el_ex_v1) || (v0 == el_ex_v1 && v1 == el_ex_v0))
			ex_el_local_edg = 0;
		else if ((v0 == el_ex_v1 && v1 == el_ex_v2) || (v0 == el_ex_v2 && v1 == el_ex_v1))
			ex_el_local_edg = 1;
		else if ((v0 == el_ex_v2 && v1 == el_ex_v0) || (v0 == el_ex_v0 && v1 == el_ex_v2))
			ex_el_local_edg = 2;
		else
		{
			printf("error in calculating exterior element local edge number \n");
			exit(EXIT_FAILURE);
		}
	
	
		// we'll just pick the normal vector to the "in" element to work with. the normal vector to the adjoining element will just be the negative of this vector
		double edg_length = sqrt((junc->x[v1]-junc->x[v0])*(junc->x[v1]-junc->x[v0]) + 
				(junc->y[v1]-junc->y[v0])*(junc->y[v1]-junc->y[v0]));
		double jac = junc->x[el_in_v0]*(junc->y[el_in_v1]-junc->y[el_in_v2]) + 
			junc->x[el_in_v1]*(junc->y[el_in_v2]-junc->y[el_in_v0]) + 
			junc->x[el_in_v2]*(junc->y[el_in_v0] - junc->y[el_in_v1]);
		double x0, x1;
		double y0, y1;
		if (int_el_local_edg == 0)
		{
			x0 = junc->x[el_in_v0];
			x1 = junc->x[el_in_v1];
			y0 = junc->y[el_in_v0];
			y1 = junc->y[el_in_v1];
	
		}
	
		else if (int_el_local_edg == 1)
		{
			x0 = junc->x[el_in_v1];
			x1 = junc->x[el_in_v2];
			y0 = junc->y[el_in_v1];
			y1 = junc->y[el_in_v2];
	
		}
	
		else if (int_el_local_edg == 2)
		{
			x0 = junc->x[el_in_v2];
			x1 = junc->x[el_in_v0];
			y0 = junc->y[el_in_v2];
			y1 = junc->y[el_in_v0];
	
		}
		nx = sign(jac)*(y1-y0)/edg_length;
		ny = sign(jac)*(x0-x1)/edg_length;
	
		// nx = junc->normvecx[index(el_in,int_el_local_edg,3)];
		// ny = junc->normvecy[index(el_in,int_el_local_edg,3)];
	
		tx = -ny;;
		ty = nx;	
	
		zeta_in = junc->zeta[index(el_in,int_el_local_edg,3)];
		Qx_in = junc->Qx[index(el_in,int_el_local_edg,3)];
		Qy_in = junc->Qy[index(el_in,int_el_local_edg,3)];
		z_edge = (junc->z[v0]+junc->z[v1])/2;
	
		// if the edge is an exterior edge connected to a channel 
		int bdrypres = junc->BdryPrescribed[i];
		// if inflow or outflow boundary, i.e. if bdrypres = 1 or 2
		if (bdrypres)	
		{
			zeta_ex = junc->bzeta[i];
			double Q_N_ex = junc->bQn[i];
			double Q_T_ex = Qx_in*tx + Qy_in*ty;
			double denom = 1./(nx*ty-ny*tx);
			Qx_ex = (ty*Q_N_ex - ny*Q_T_ex)*denom;
			Qy_ex = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
		}
	
		// if edge i is a boundary edge but not connected to a channel, implement no flux boundary condition
		else if (el_in == el_ex)
		{
			zeta_ex = junc->zeta[index(el_ex, ex_el_local_edg,3)];
			// Compute the velocity in the normal direction
			double Q_N_in = Qx_in*nx + Qy_in*ny;
			double Q_T_in = Qx_in*tx + Qy_in*ty;
	
			// Reflect the velocity in the normal direction
			double Q_N_ex = -Q_N_in;
			double Q_T_ex = Q_T_in;
	
			// Compute the x and y components of the external state flow
			double denom = 1./(nx*ty - ny*tx);
			Qx_ex = (ty*Q_N_ex - ny*Q_T_ex)*denom;
			Qy_ex = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
		
		}
			
		// if the edge is not a boundary edge
		else
		{
			zeta_ex = junc->zeta[index(el_ex,ex_el_local_edg,3)];
			Qx_ex = junc->Qx[index(el_ex,ex_el_local_edg,3)];
			Qy_ex = junc->Qy[index(el_ex,ex_el_local_edg,3)];
		}
	
		#ifdef WDON
		// Check to see if both of the elements separated by this edge are dry
		if (junc->WD[el_in] == 0 && junc->WD[el_ex] == 0)
		{
			// Reflection flux for the interior element
			double zeta_ex_ref = zeta_in;
			double Q_N_in = Qx_in*nx + Qy_in*ny;
			double Q_T_in = Qx_in*tx + Qy_in*ty;
		
			double Q_N_ex = -Q_N_in;
			double Q_T_ex = Q_T_in;
			
			double denom = 1./(nx*ty - ny*tx);
			double Qx_ex_ref = (ty*Q_N_ex - ny*Q_T_ex)*denom;
			double Qy_ex_ref = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
		
			double max_lam_in = RoeFlux2D(zeta_in, zeta_ex_ref, Qx_in,Qx_ex_ref,
				Qy_in, Qy_ex_ref, z_edge, nx, ny, 0, Fn_in);
		
			// Reflection flux for the exterior element
			double zeta_in_ref = zeta_ex;
			Q_N_ex = Qx_ex*nx + Qy_ex*ny;
			Q_T_ex = Qx_ex*tx + Qy_ex*ty;
		
			Q_N_in = -Q_N_ex;
			Q_T_in = Q_T_ex;
			
			double Qx_in_ref = (ty*Q_N_in - ny*Q_T_in)*denom;
			double Qy_in_ref = (-tx*Q_N_in + nx*Q_T_ex)*denom;
			double max_lam_ex = RoeFlux2D(zeta_in_ref, zeta_ex, Qx_in_ref, Qx_ex, Qy_in_ref, Qy_ex, z_edge, nx, ny, 0, Fn_ex);
		
			if (i ==0)
				junc->max_lambda = max(max_lam_in, max_lam_ex);
			else
			{
				junc->max_lambda = fmax(junc->max_lambda, max_lam_in);
				junc->max_lambda = fmax(junc->max_lambda, max_lam_ex);
			}
		
		}
		else // if the elements aren't both dry
		{
			double Fhatdotn[3];
			double current_max_lam= RoeFlux2D(zeta_in, zeta_ex, Qx_in, Qx_ex, Qy_in, Qy_ex, z_edge, nx, ny, g,Fhatdotn);
			if (isnan(Fhatdotn[0]) || isnan(Fhatdotn[1]) || isnan(Fhatdotn[2]))
			{
				printf("2D both elements wet flux not a number, edge %d, time %e \n",i,time);
				exit(EXIT_FAILURE);
			}
			
			if (junc->WD[el_in] == 1)
			{
				for (int j =0; j < 3; ++j)
				{
					Fn_in[j] = Fhatdotn[j];
				}
			}
			else
			{
				double max_lam_in = RoeFlux2D(zeta_in, zeta_ex, Qx_in, Qx_ex, Qy_in, Qy_ex,
					z_edge, nx, ny, 0, Fn_in);
				if (isnan(Fn_in[0]) || isnan(Fn_in[1]) || isnan(Fn_in[2]))
				{
					printf("2D interior element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
				current_max_lam = max(current_max_lam, max_lam_in);
			}
	
			if (junc->WD[el_ex] == 1)
			{
				for (int j = 0; j < 3; ++j)
				{
					Fn_ex[j] = Fhatdotn[j];
				}
			}
			else
			{
				double max_lam_ex = RoeFlux2D(zeta_in, zeta_ex, Qx_in, Qx_ex,
					Qy_in, Qy_ex, z_edge, nx, ny, 0, Fn_ex);
				if (isnan(Fn_ex[0]) || isnan(Fn_ex[1]) || isnan(Fn_ex[2]))
				{
					printf("2D exterior element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
				current_max_lam = max(current_max_lam, max_lam_ex);
			}
				
			if (i == 0)
				junc->max_lambda = current_max_lam;
			else
				junc->max_lambda = fmax(junc->max_lambda, current_max_lam);
		}
			
		#else
		double Fhatdotn[3];
		double current_max_lam= RoeFlux2D(zeta_in, zeta_ex, Qx_in, Qx_ex, Qy_in, Qy_ex, z_edge, nx, ny, g, Fhatdotn);
		if (isnan(Fhatdotn[0]) || isnan(Fhatdotn[1]) || isnan(Fhatdotn[2]))
		{
			printf("2D numerical flux not a number, edge %d \n",i);
			printf("zeta_in = %e, zeta_ex = %e, BdryPrescribed = %d \n", zeta_in, zeta_ex, bdrypres);
			printf("el1 = %d, el2 = %d vrt1 = %d, vrt2 = %d\n", el_in, el_ex, v1, v0);
			exit(EXIT_FAILURE);
		}
			
		for (int j =0; j < 3; ++j)
		{
			Fn_in[j] = Fhatdotn[j];
			Fn_ex[j] = Fhatdotn[j];
		}
	
		if (i == 0)
			junc->max_lambda = current_max_lam;
		else
			junc->max_lambda = max(junc->max_lambda, current_max_lam);
	
		#endif
	
	
		// store Fhatdotn for the two elements connected by this edge
		Fhat1dotn[index(el_in,int_el_local_edg,3)] = Fn_in[0];
		Fhat2dotn[index(el_in,int_el_local_edg,3)] = Fn_in[1];
		Fhat3dotn[index(el_in,int_el_local_edg,3)] = Fn_in[2];
	
		// store the value for the exterior element, only if the edge is not a boundary edge
		if (el_ex != el_in)
		{
			Fhat1dotn[index(el_ex,ex_el_local_edg,3)] = -Fn_ex[0];
			Fhat2dotn[index(el_ex,ex_el_local_edg,3)] = -Fn_ex[1];
			Fhat3dotn[index(el_ex,ex_el_local_edg,3)] = -Fn_ex[2];
		}
	
		
	}

	//double SI[NumEl][9],VI[NumEl][21];
	
	//Compute2DInnerProducts(junc, SI, VI,time);
	
	
	for (int k=0; k<NumEl; ++k)
	{

		double SI[9], VI[21];
		Compute2DInnerProducts(junc, Fhat1dotn, Fhat2dotn, Fhat3dotn, k, SI, VI);

		double zeta1 = -SI[0] + VI[0];
		double zeta2 = -SI[1] + VI[1];
		double zeta3 = -SI[2] + VI[2];
	
		double Qx1 = -SI[3] + VI[3] + VI[4] - VI[5];
		double Qx2 = -SI[4] + VI[6] + VI[7] - VI[8];
		double Qx3 = -SI[5] + VI[9] + VI[10] - VI[11];
	
		double Qy1 = -SI[6] + VI[12] + VI[13] - VI[14];
		double Qy2 = -SI[7] + VI[15] + VI[16] - VI[17];
		double Qy3 = -SI[8] + VI[18] + VI[19] - VI[20];
	
		double zeta_el[3] = {zeta1, zeta2, zeta3};
		double Qx_el[3] = {Qx1, Qx2, Qx3};
		double Qy_el[3] = {Qy1, Qy2, Qy3};
	
		double x0 = junc->x[junc->EltoVert[k*3]];
		double x1 = junc->x[junc->EltoVert[k*3+1]];
		double x2 = junc->x[junc->EltoVert[k*3+2]];
		double y0 = junc->y[junc->EltoVert[k*3]];
		double y1 = junc->y[junc->EltoVert[k*3+1]];
		double y2 = junc->y[junc->EltoVert[k*3+2]];
	
		double jac = x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1);
	
		double M11 = 6/fabs(jac);
		double M12 = 0;
		double M13 = 0;
		double M21 = 0;
		double M22 = M11;
		double M23 = 0;
		double M31 = 0;
		double M32 = 0;
		double M33 = M11;
		double Minv[3][3] = {{M11,M12,M13},{M21,M22,M23},{M31,M32,M33}};
	
		double RHSZeta_el[3];
		double RHSQx_el[3];
		double RHSQy_el[3];
	
		// RHSZeta = invM*[zeta1;zeta2;zeta33] and so on
		for (int j = 0; j<3; ++j)
		{
			RHSZeta_el[j] = Minv[j][0]*zeta_el[0]+Minv[j][1]*zeta_el[1] + Minv[j][2]*zeta_el[2];
			RHSQx_el[j] = Minv[j][0]*Qx_el[0]+Minv[j][1]*Qx_el[1] + Minv[j][2]*Qx_el[2];
			RHSQy_el[j] = Minv[j][0]*Qy_el[0]+Minv[j][1]*Qy_el[1] + Minv[j][2]*Qy_el[2];
	
		}

		for (int j=0; j<3; ++j){
			RHSZeta[index(k,j,3)] = RHSZeta_el[j];
			RHSQx[index(k,j,3)] = RHSQx_el[j];
			RHSQy[index(k,j,3)] = RHSQy_el[j];
		}
	
	}

}

 



