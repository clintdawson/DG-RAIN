#include <stdio.h>
#include <math.h>
#include "math_functions.h"
#include "constitutive_equations.h"
#include "Globals.h"
#include "Watershed.h"

extern void GetLGLWeights(int N, double* w);

void minmodChan(struct channel* Chan)
{
	int NumEl = Chan->NumEl;
	double *dh = Chan->dh;
	int Np = Chan->Np;
	int P = Chan->P;
	int NumNodes = Chan->NumNodes;

	// average value of Zeta and Q
	double* avgZeta = xcalloc(NumEl, sizeof(double));
	double* avgQ = xcalloc(NumEl, sizeof(double));
	double** Zeta = xcalloc(NumEl, sizeof(double*));
	
	for (int k = 0; k < NumEl; k++)
		Zeta[k] = xcalloc(Np, sizeof(double));

	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i=0; i < NumEl; ++i)
	{
		double avgZetaVal = 0; 
		double avgQVal = 0;
		for (int j = 0; j < Np; j++)
		{
			int nodeNum = i*Np+j;
			double Aval = Chan->A[nodeNum+1];
			double bval = Chan->NodalB[nodeNum];
			double zval = Chan->NodalZ[nodeNum];
			double m1val = Chan->Nodalm1[nodeNum];
			double m2val = Chan->Nodalm2[nodeNum];
			double Hval = getH(Aval, bval, m1val, m2val);
			double Qval = Chan->Q[nodeNum+1];
			double zetaVal = Hval - zval;
			Zeta[i][j] = zetaVal;

			avgZetaVal += LGLWeight[j]*zetaVal;
			avgQVal += LGLWeight[j]*Qval;
		}
		avgZeta[i] = 0.5*avgZetaVal;
		avgQ[i] = 0.5*avgQVal;

	}

	double eps0 = 1e-8;
	for (int i =0; i < NumEl; ++i)
	{	
		double a1, a2, b1, b2, sigma1, sigma2;
		double Zetatilda[2], Qtilda[2];
		if (i > 0)
		{
			a1 = 2*(avgZeta[i] - avgZeta[i-1])/(dh[i] + dh[i-1]);
			a2 = 2*(avgQ[i] - avgQ[i-1])/(dh[i]+dh[i-1]);
		}
		else
		{	
			a1 = 0;
			a2 = 0;
		}
		if (i < NumEl-1)
		{	
			b1 = 2*(avgZeta[i+1] - avgZeta[i])/(dh[i]+dh[i+1]);
			b2 = 2*(avgQ[i+1] - avgQ[i])/(dh[i]+dh[i+1]);
		}
		else
		{
			b1 = 0; 
			b2 = 0;
		}

		sigma1 = sign(a1)*0.5*(1+sign(a1*b1))*fmin(fabs(a1),fabs(b1));
		sigma2 = sign(a2)*0.5*(1+sign(a2*b2))*fmin(fabs(a2),fabs(b2));

		Zetatilda[0] = avgZeta[i] - 0.5*sigma1*dh[i];
		Zetatilda[1] = avgZeta[i] + 0.5*sigma1*dh[i];
		Qtilda[0] = avgQ[i] - 0.5*sigma2*dh[i];
		Qtilda[1] = avgQ[i] + 0.5*sigma2*dh[i];

		// limit only if limiting is required
		if (fabs(Zetatilda[0]-Zeta[i][0]) > eps0 || fabs(Zetatilda[1] - Zeta[i][Np-1]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Zeta[i][j] = Zetatilda[0] + j*(Zetatilda[1] - Zetatilda[0])/(Np-1);
				double hval = Zeta[i][j] + Chan->NodalZ[i*Np+j];
				double bval = Chan->NodalB[i*Np+j];

				double m1val = Chan->Nodalm1[i*Np+j];
				double m2val = Chan->Nodalm2[i*Np+j];
				Chan->A[i*Np+j+1] = hval*bval + 0.5*m1val*hval*hval+ 0.5*m2val*hval*hval;
			}
		}
	
		// limit only if limiting is required
		if (fabs(Qtilda[0]-Chan->Q[i*Np+1]) > eps0 || fabs(Qtilda[1] - Chan->Q[i*Np+Np]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Chan->Q[i*Np+j+1] = Qtilda[0] + j*(Qtilda[1] - Qtilda[0])/(Np-1);
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
	//Chan->A[0] = Chan->A[1];
	//Chan->Q[0] = Chan->Q[1];
	//Chan->A[NumNodes-1] = Chan->A[NumNodes-2];
	//Chan->Q[NumNodes-1] = Chan->Q[NumNodes-2];

	for (int k = 0; k < NumEl; k++)
		free(Zeta[k]);

	free(Zeta);
	free(avgZeta);
	free(avgQ);
}

void minmodHeightChan(struct channel* Chan)
{
	int NumEl = Chan->NumEl;
	double *dh = Chan->dh;
	int Np = Chan->Np;
	int P = Chan->P;
	int NumNodes = Chan->NumNodes;

	// average value of Zeta and Q
	double* avgHeight = xcalloc(NumEl, sizeof(double));
	double* avgQ = xcalloc(NumEl, sizeof(double));
	double** Height = xcalloc(NumEl, sizeof(double*));
	
	for (int k = 0; k < NumEl; k++)
		Height[k] = xcalloc(Np, sizeof(double));

	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i=0; i < NumEl; ++i)
	{
		double avgHeightVal = 0; 
		double avgQVal = 0;
		for (int j = 0; j < Np; j++)
		{
			int nodeNum = i*Np+j;
			double Aval = Chan->A[nodeNum+1];
			double bval = Chan->NodalB[nodeNum];
			double m1val = Chan->Nodalm1[nodeNum];
			double m2val = Chan->Nodalm2[nodeNum];
			double Hval = getH(Aval, bval, m1val, m2val);
			double Qval = Chan->Q[nodeNum+1];
			Height[i][j] = Hval;

			avgHeightVal += LGLWeight[j]*Hval;
			avgQVal += LGLWeight[j]*Qval;
		}
		avgHeight[i] = 0.5*avgHeightVal;
		avgQ[i] = 0.5*avgQVal;

	}

	double eps0 = 1e-8;
	for (int i =0; i < NumEl; ++i)
	{	
		double a1, a2, b1, b2, sigma1, sigma2;
		double Htilda[2], Qtilda[2];
		if (i > 0)
		{
			a1 = 2*(avgHeight[i] - avgHeight[i-1])/(dh[i] + dh[i-1]);
			a2 = 2*(avgQ[i] - avgQ[i-1])/(dh[i]+dh[i-1]);
		}
		else
		{	
			a1 = 0;
			a2 = 0;
		}
		if (i < NumEl-1)
		{	
			b1 = 2*(avgHeight[i+1] - avgHeight[i])/(dh[i]+dh[i+1]);
			b2 = 2*(avgQ[i+1] - avgQ[i])/(dh[i]+dh[i+1]);
		}
		else
		{
			b1 = 0; 
			b2 = 0;
		}

		sigma1 = sign(a1)*0.5*(1+sign(a1*b1))*fmin(fabs(a1),fabs(b1));
		sigma2 = sign(a2)*0.5*(1+sign(a2*b2))*fmin(fabs(a2),fabs(b2));

		Htilda[0] = avgHeight[i] - 0.5*sigma1*dh[i];
		Htilda[1] = avgHeight[i] + 0.5*sigma1*dh[i];
		Qtilda[0] = avgQ[i] - 0.5*sigma2*dh[i];
		Qtilda[1] = avgQ[i] + 0.5*sigma2*dh[i];

		// limit only if limiting is required
		if (fabs(Htilda[0]-Height[i][0]) > eps0 || fabs(Htilda[1] - Height[i][Np-1]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				double hval = Htilda[0] + j*(Htilda[1] - Htilda[0])/(Np-1);
				Height[i][j] = hval;
				double bval = Chan->NodalB[i*Np+j];

				double m1val = Chan->Nodalm1[i*Np+j];
				double m2val = Chan->Nodalm2[i*Np+j];
				Chan->A[i*Np+j+1] = hval*bval + 0.5*m1val*hval*hval+ 0.5*m2val*hval*hval;
			}
		}
	
		// limit only if limiting is required
		if (fabs(Qtilda[0]-Chan->Q[i*Np+1]) > eps0 || fabs(Qtilda[1] - Chan->Q[i*Np+Np]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Chan->Q[i*Np+j+1] = Qtilda[0] + j*(Qtilda[1] - Qtilda[0])/(Np-1);
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
	Chan->A[0] = Chan->A[1];
	Chan->Q[0] = Chan->Q[1];
	Chan->A[NumNodes-1] = Chan->A[NumNodes-2];
	Chan->Q[NumNodes-1] = Chan->Q[NumNodes-2];

	for (int k = 0; k < NumEl; k++)
		free(Height[k]);

	free(Height);
	free(avgHeight);
	free(avgQ);
}



/************************ Slope limiter reduces convergence rate to 1.5 ******************/


void minmodKinematicEls(int fp)
{
	int NumEl = FloodplainList[fp]->NumEl;

	double* avgH = xcalloc(NumEl, sizeof(double));

	for (int k =0; k < NumEl; k++)
	{
		struct kinematicEl* kinEl = KinematicElList[k];
		if (kinEl->isActive)
		{
			int P = kinEl->P;
			int Np = kinEl->Np;
			double LGLWeight[Np];
			GetLGLWeights(P, LGLWeight);
		
			double avgVal = 0;
			for (int i =0; i < Np; i++)
			{
				double A = kinEl->A[i];
				double height = A/kinEl->weq;
				avgVal += LGLWeight[i]*height;
			}
			avgH[k] = 0.5*avgVal;
			//printf("avgHt = %3.13f\n", avgHt[k]);
		}
	}

	double eps0 = 1e-8;
	
	for (int k =0; k < NumEl; k++)
	{
		struct kinematicEl* kinEl = KinematicElList[k];
		if (kinEl->isActive)
		{
			int Np = kinEl->Np;
			int numUpstreamEl = kinEl->numUpstreamEls;
			int numDownstreamEls = kinEl->numDownstreamEls;
			double dh = kinEl->dh;
			double a, b, sigma, Htilda[2];
			int upstreamEl, downstreamEl;

			double avgUpstreamH = 0;
			if (numUpstreamEl > 0)
			{
				for (int j = 0; j < numUpstreamEl; j++)
				{
					upstreamEl = kinEl->upstreamEls[j];
					avgUpstreamH += avgH[upstreamEl];
				}
				avgUpstreamH /= numUpstreamEl;
				a = 2*(avgH[k] - avgUpstreamH)/(dh+KinematicElList[upstreamEl]->dh); 
			}
			else
				a = 0;

			if (numDownstreamEls == 1)
			{
				downstreamEl = kinEl->el2;
				b = 2*(avgH[downstreamEl] - avgH[k])/(dh + KinematicElList[downstreamEl]->dh);

			}
			else 
				b = 0;

			sigma = sign(a)*0.5*(1+sign(a*b))*fmin(fabs(a),fabs(b));
		
			Htilda[0] = avgH[k] - 0.5*sigma*dh;
			Htilda[1] = avgH[k] + 0.5*sigma*dh;

			double weq = kinEl->weq;
			if (fabs(Htilda[0] - kinEl->A[0]/weq) > eps0 || fabs(Htilda[1] - kinEl->A[Np-1]/weq) > eps0)
			{
				kinEl->A[0] = Htilda[0]*weq;
				for (int j =0; j < Np-1; j++)
				{
					kinEl->A[j] = (Htilda[0] + j*(Htilda[1] - Htilda[0])/(Np-1))*weq;
				}
				kinEl->A[Np-1] = Htilda[1]*weq;
			}
		}

	}

	free(avgH);
}

/*************************************************************************//**
*
* Routine to apply the BDS slope limiter on the 2-D conserved variables
* @param[in] junc a pointer to the junction structure that we are currently
* working on
*
* **************************************************************************/
// Slopelimiter currently only works for P = 1
void SlopeLimiter2D(struct TwoDRegion *junc)
{	
	int NumVerts = junc->NumVerts;
	int NumEl = junc->NumEl;

	int *ElCount = junc->ElCount;
	int **VtoEl = junc->VtoEl;
	int **VtoNode = junc->VtoNode;

	// Find the maximum and minimum of the average of each variable over all elements sharing a node
	double *Zeta_Max = malloc(NumVerts*sizeof(double));
	double *Qx_Max = malloc(NumVerts*sizeof(double));
	double *Qy_Max = malloc(NumVerts*sizeof(double));
	double *Zeta_Min = malloc(NumVerts*sizeof(double));
	double *Qx_Min = malloc(NumVerts*sizeof(double));
	double *Qy_Min = malloc(NumVerts*sizeof(double));

	for (int i =0; i < NumVerts; ++i)
	{
		Zeta_Max[i] = -99999;
		Zeta_Min[i] = 99999;
		Qx_Max[i] = -99999;
		Qx_Min[i] = 99999;
		Qy_Max[i] = -99999;
		Qy_Min[i] = 99999;

		int NumConnEl = ElCount[i];
		// loop over all the elements that share node i
		for (int j = 0; j < NumConnEl; ++j)
		{
			int el = VtoEl[i][j];
			int lv1 = VtoNode[el][0];
			int lv2 = VtoNode[el][1];
			int lv3 = VtoNode[el][2];

			double avgZeta = (junc->zeta[el][lv1] + junc->zeta[el][lv2] + junc->zeta[el][lv3])/3;
			double avgQx = (junc->Qx[el][lv1] + junc->Qx[el][lv2] + junc->Qx[el][lv3])/3;
			double avgQy = (junc->Qy[el][lv1] + junc->Qy[el][lv2] + junc->Qy[el][lv3])/3;

 			if (avgZeta < Zeta_Min[i])
				Zeta_Min[i] = avgZeta;
			if (avgZeta > Zeta_Max[i])
				Zeta_Max[i] = avgZeta;
			if (avgQx < Qx_Min[i])
				Qx_Min[i] = avgQx;
			if (avgQx > Qx_Max[i])
				Qx_Max[i] = avgQx;
			if (avgQy < Qy_Min[i])
				Qy_Min[i] = avgQy;
			if (avgQy > Qy_Max[i])
				Qy_Max[i] = avgQy;

		}

	}	


	// Loop over elements to calculate the new vertex values
	for (int el = 0; el < NumEl; ++el)
	{	
		double ZeVertex[3], ZeMin[3], ZeMax[3];
		double ZeAvg;

		int v1 = junc->EltoVert[el*3];
		int v2 = junc->EltoVert[el*3+1];
		int v3 = junc->EltoVert[el*3+2];

		int n1 = junc->VtoNode[el][0];
		int n2 = junc->VtoNode[el][1];
		int n3 = junc->VtoNode[el][2];

		int modify = 1;

		// loop over variables
		for (int ivar = 0; ivar < 3; ++ivar)
		{

			if (ivar == 0) 
			{
				ZeVertex[0] = junc->zeta[el][n1]; 
				ZeVertex[1] = junc->zeta[el][n2];
				ZeVertex[2] = junc->zeta[el][n3];
				ZeMax[0] = Zeta_Max[v1];
				ZeMax[1] = Zeta_Max[v2];
				ZeMax[2] = Zeta_Max[v3];
				ZeMin[0] = Zeta_Min[v1];
				ZeMin[1] = Zeta_Min[v2];
				ZeMin[2] = Zeta_Min[v3];

				// fix this section later
				//#ifdef WDON
				//double z1 = junc->Vz[v1];
				//double z2 = junc->Vz[v2];
				//double z3 = junc->Vz[v3];
				//double H1 = ZeVertex[0] + z1;
				//double H2 = ZeVertex[1] + z2;
				//double H3 = ZeVertex[2] + z3;

				//double Hmin0 = ZeMin[0] + z1;
				//double Hmiv1 = ZeMin[1] + z2;
				//double Hmiv2 = ZeMin[2] + z3;
				//if (H1 > H0 && H2 > H0 && H3 > H0 && Hmin0 > H0 && Hmiv1 > H0 && Hmiv2 > H0)
				//	modify = 1;
				//else
				//	modify = 0;
				//#endif

			}

			else if (ivar == 1)
			{
				ZeVertex[0] = junc->Qx[el][n1];
				ZeVertex[1] = junc->Qx[el][n2];
				ZeVertex[2] = junc->Qx[el][n3];
				ZeMax[0] = Qx_Max[v1];
				ZeMax[1] = Qx_Max[v2];
				ZeMax[2] = Qx_Max[v3];
				ZeMin[0] = Qx_Min[v1];
				ZeMin[1] = Qx_Min[v2];
				ZeMin[2] = Qx_Min[v3];

			}		
	
			else
			{
				ZeVertex[0] = junc->Qy[el][n1];
				ZeVertex[1] = junc->Qy[el][n2];
				ZeVertex[2] = junc->Qy[el][n3];
				ZeMax[0] = Qy_Max[v1];
				ZeMax[1] = Qy_Max[v2];
				ZeMax[2] = Qy_Max[v3];
				ZeMin[0] = Qy_Min[v1];
				ZeMin[1] = Qy_Min[v2];
				ZeMin[2] = Qy_Min[v3];

			}	

			ZeAvg = (ZeVertex[0] + ZeVertex[1] +ZeVertex[2])/3;
			// Reset the vertex value to be less than or equal to the max
			// and greater than or equal to the min at that vertex
			ZeVertex[0] = max(min(ZeVertex[0],ZeMax[0]),ZeMin[0]);
			ZeVertex[1] = max(min(ZeVertex[1],ZeMax[1]),ZeMin[1]);
			ZeVertex[2] = max(min(ZeVertex[2],ZeMax[2]),ZeMin[2]);


			// Loop over the vertices 3 times
			// If the value at the vertex is above (below) the max (min) at that vertex
			// then subtract off the difference and add it to the other vertices

			double diff[3];
			for (int ll= 0; ll < 3; ++ll)
			{
				double sumloc = (ZeVertex[0] + ZeVertex[1] + ZeVertex[2])/3;
				double sumdiff = (sumloc - ZeAvg)*3;
				int signdiff = sign(sumdiff);
				diff[0] = (ZeVertex[0] - ZeAvg)*signdiff;
				diff[1] = (ZeVertex[1] - ZeAvg)*signdiff;
				diff[2] = (ZeVertex[2] - ZeAvg)*signdiff;
				int inc1 = 0;
				if (diff[0] > 0)
					inc1 = 1;
				int inc2 = 0;
				if (diff[1] > 0)
					inc2 = 1;
				int inc3 = 0;
				if (diff[2] > 0)
					inc3 = 1;

				double kdp = inc1 + inc2 + inc3;
				
				for (int k = 0; k < 3; ++k)
				{
					double div = max(1, kdp);
					double redfac;
					double redmax;

					if (diff[k] > 0)
					{
						redfac = sumdiff*signdiff/div;
						kdp -= 1;
					} 
					else
						redfac = 0;

					if (signdiff > 0)
						redmax = ZeVertex[k] - ZeMin[k];
					else
						redmax = ZeMax[k] - ZeVertex[k];
					
					redfac = min(redfac,redmax);
					sumdiff -= redfac*signdiff;
					ZeVertex[k] -= redfac*signdiff;

				} // end k loop

			} // end ll loop

			if (ivar == 0)
			{
				//#ifdef WDON
				//double H1 = ZeVertex[0] + junc->Vz[v1];
				//double H2 = ZeVertex[1] + junc->Vz[v2];
				//double H3 = ZeVertex[2] + junc->Vz[v3];
				//if (H1 < H0 || H2 < H0 || H3 < H0)
				//	modify = 0;
				//#endif
				if (modify)
				{
					junc->zeta[el][n1] = ZeVertex[0];
					junc->zeta[el][n2] = ZeVertex[1];
					junc->zeta[el][n3] = ZeVertex[2];
				}

			}	

			if (ivar == 1)
			{
				if (modify)
				{
					junc->Qx[el][n1] = ZeVertex[0];
					junc->Qx[el][n2] = ZeVertex[1];
					junc->Qx[el][n3] = ZeVertex[2];
				}
			}	

			if (ivar == 2)
			{
				if (modify)
				{
					junc->Qy[el][n1] = ZeVertex[0];
					junc->Qy[el][n2] = ZeVertex[1];
					junc->Qy[el][n3] = ZeVertex[2];
				}
			}
	
		} // end loop over variables

	} // end loop over elements

	free(Zeta_Max);
	free(Qx_Max);
	free(Qy_Max);
	free(Zeta_Min);
	free(Qx_Min);
	free(Qy_Min);

}


/*************************************************************************//**
*
* Routine to apply the BDS slope limiter on the 2-D conserved variables
* @param[in] junc a pointer to the junction structure that we are currently
* working on (water height instead of water surface elevation)
*
* **************************************************************************/
// Slopelimiter currently only works for P = 1
void SlopeLimiterHeightJunc(struct TwoDRegion *junc)
{	
	int NumVerts = junc->NumVerts;
	int NumEl = junc->NumEl;

	int *ElCount = junc->ElCount;
	int **VtoEl = junc->VtoEl;
	int **VtoNode = junc->VtoNode;

	// Find the maximum and minimum of the average of each variable over all elements sharing a node
	double *H_Max = malloc(NumVerts*sizeof(double));
	double *Qx_Max = malloc(NumVerts*sizeof(double));
	double *Qy_Max = malloc(NumVerts*sizeof(double));
	double *H_Min = malloc(NumVerts*sizeof(double));
	double *Qx_Min = malloc(NumVerts*sizeof(double));
	double *Qy_Min = malloc(NumVerts*sizeof(double));

	for (int i =0; i < NumVerts; ++i)
	{
		H_Max[i] = -99999;
		H_Min[i] = 99999;
		Qx_Max[i] = -99999;
		Qx_Min[i] = 99999;
		Qy_Max[i] = -99999;
		Qy_Min[i] = 99999;

		int NumConnEl = ElCount[i];
		// loop over all the elements that share node i
		for (int j = 0; j < NumConnEl; ++j)
		{
			int el = VtoEl[i][j];
			int lv1 = VtoNode[el][0];
			int lv2 = VtoNode[el][1];
			int lv3 = VtoNode[el][2];

			double h1 = junc->zeta[el][lv1] + junc->NodalZ[el][lv1];
			double h2 = junc->zeta[el][lv2] + junc->NodalZ[el][lv2];
			double h3 = junc->zeta[el][lv3] + junc->NodalZ[el][lv3];

			double avgH = (h1+h2+h3)/3;
			double avgQx = (junc->Qx[el][lv1] + junc->Qx[el][lv2] + junc->Qx[el][lv3])/3;
			double avgQy = (junc->Qy[el][lv1] + junc->Qy[el][lv2] + junc->Qy[el][lv3])/3;

 			if (avgH < H_Min[i])
				H_Min[i] = avgH;
			if (avgH > H_Max[i])
				H_Max[i] = avgH;
			if (avgQx < Qx_Min[i])
				Qx_Min[i] = avgQx;
			if (avgQx > Qx_Max[i])
				Qx_Max[i] = avgQx;
			if (avgQy < Qy_Min[i])
				Qy_Min[i] = avgQy;
			if (avgQy > Qy_Max[i])
				Qy_Max[i] = avgQy;

		}

	}	


	// Loop over elements to calculate the new vertex values
	for (int el = 0; el < NumEl; ++el)
	{	
		double ZeVertex[3], ZeMin[3], ZeMax[3];
		double ZeAvg;

		int v1 = junc->EltoVert[el*3];
		int v2 = junc->EltoVert[el*3+1];
		int v3 = junc->EltoVert[el*3+2];

		int n1 = junc->VtoNode[el][0];
		int n2 = junc->VtoNode[el][1];
		int n3 = junc->VtoNode[el][2];

		int modify = 1;

		// loop over variables
		for (int ivar = 0; ivar < 3; ++ivar)
		{

			if (ivar == 0) 
			{
				ZeVertex[0] = junc->zeta[el][n1] + junc->NodalZ[el][n1]; 
				ZeVertex[1] = junc->zeta[el][n2] + junc->NodalZ[el][n2];
				ZeVertex[2] = junc->zeta[el][n3] + junc->NodalZ[el][n3];
				ZeMax[0] = H_Max[v1];
				ZeMax[1] = H_Max[v2];
				ZeMax[2] = H_Max[v3];
				ZeMin[0] = H_Min[v1];
				ZeMin[1] = H_Min[v2];
				ZeMin[2] = H_Min[v3];

				// fix this section later
				#ifdef WDON
				double z1 = junc->Vz[v1];
				double z2 = junc->Vz[v2];
				double z3 = junc->Vz[v3];
				double H1 = ZeVertex[0] + z1;
				double H2 = ZeVertex[1] + z2;
				double H3 = ZeVertex[2] + z3;

				double Hmin0 = ZeMin[0] + z1;
				double Hmiv1 = ZeMin[1] + z2;
				double Hmiv2 = ZeMin[2] + z3;
				if (H1 > H0 && H2 > H0 && H3 > H0 && Hmin0 > H0 && Hmiv1 > H0 && Hmiv2 > H0)
					modify = 1;
				else
					modify = 0;
				#endif

			}

			else if (ivar == 1)
			{
				ZeVertex[0] = junc->Qx[el][n1];
				ZeVertex[1] = junc->Qx[el][n2];
				ZeVertex[2] = junc->Qx[el][n3];
				ZeMax[0] = Qx_Max[v1];
				ZeMax[1] = Qx_Max[v2];
				ZeMax[2] = Qx_Max[v3];
				ZeMin[0] = Qx_Min[v1];
				ZeMin[1] = Qx_Min[v2];
				ZeMin[2] = Qx_Min[v3];

			}		
	
			else
			{
				ZeVertex[0] = junc->Qy[el][n1];
				ZeVertex[1] = junc->Qy[el][n2];
				ZeVertex[2] = junc->Qy[el][n3];
				ZeMax[0] = Qy_Max[v1];
				ZeMax[1] = Qy_Max[v2];
				ZeMax[2] = Qy_Max[v3];
				ZeMin[0] = Qy_Min[v1];
				ZeMin[1] = Qy_Min[v2];
				ZeMin[2] = Qy_Min[v3];

			}	

			ZeAvg = (ZeVertex[0] + ZeVertex[1] +ZeVertex[2])/3;
			// Reset the vertex value to be less than or equal to the max
			// and greater than or equal to the min at that vertex
			ZeVertex[0] = max(min(ZeVertex[0],ZeMax[0]),ZeMin[0]);
			ZeVertex[1] = max(min(ZeVertex[1],ZeMax[1]),ZeMin[1]);
			ZeVertex[2] = max(min(ZeVertex[2],ZeMax[2]),ZeMin[2]);


			// Loop over the vertices 3 times
			// If the value at the vertex is above (below) the max (min) at that vertex
			// then subtract off the difference and add it to the other vertices

			double diff[3];
			for (int ll= 0; ll < 3; ++ll)
			{
				double sumloc = (ZeVertex[0] + ZeVertex[1] + ZeVertex[2])/3;
				double sumdiff = (sumloc - ZeAvg)*3;
				int signdiff = sign(sumdiff);
				diff[0] = (ZeVertex[0] - ZeAvg)*signdiff;
				diff[1] = (ZeVertex[1] - ZeAvg)*signdiff;
				diff[2] = (ZeVertex[2] - ZeAvg)*signdiff;
				int inc1 = 0;
				if (diff[0] > 0)
					inc1 = 1;
				int inc2 = 0;
				if (diff[1] > 0)
					inc2 = 1;
				int inc3 = 0;
				if (diff[2] > 0)
					inc3 = 1;

				double kdp = inc1 + inc2 + inc3;
				
				for (int k = 0; k < 3; ++k)
				{
					double div = max(1, kdp);
					double redfac;
					double redmax;

					if (diff[k] > 0)
					{
						redfac = sumdiff*signdiff/div;
						kdp -= 1;
					} 
					else
						redfac = 0;

					if (signdiff > 0)
						redmax = ZeVertex[k] - ZeMin[k];
					else
						redmax = ZeMax[k] - ZeVertex[k];
					
					redfac = min(redfac,redmax);
					sumdiff -= redfac*signdiff;
					ZeVertex[k] -= redfac*signdiff;

				} // end k loop

			} // end ll loop

			if (ivar == 0)
			{
				#ifdef WDON
				double H1 = ZeVertex[0] + junc->Vz[v1];
				double H2 = ZeVertex[1] + junc->Vz[v2];
				double H3 = ZeVertex[2] + junc->Vz[v3];
				if (H1 < H0 || H2 < H0 || H3 < H0)
					modify = 0;
				#endif
				if (modify)
				{
					junc->zeta[el][n1] = ZeVertex[0] - junc->NodalZ[el][n1];
					junc->zeta[el][n2] = ZeVertex[1] - junc->NodalZ[el][n2];
					junc->zeta[el][n3] = ZeVertex[2] - junc->NodalZ[el][n2];
				}

			}	

			if (ivar == 1)
			{
				if (modify)
				{
					junc->Qx[el][n1] = ZeVertex[0];
					junc->Qx[el][n2] = ZeVertex[1];
					junc->Qx[el][n3] = ZeVertex[2];
				}
			}	

			if (ivar == 2)
			{
				if (modify)
				{
					junc->Qy[el][n1] = ZeVertex[0];
					junc->Qy[el][n2] = ZeVertex[1];
					junc->Qy[el][n3] = ZeVertex[2];
				}
			}
	
		} // end loop over variables

	} // end loop over elements

	free(H_Max);
	free(Qx_Max);
	free(Qy_Max);
	free(H_Min);
	free(Qx_Min);
	free(Qy_Min);

}


