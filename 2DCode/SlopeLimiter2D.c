#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MeshAttributes.h"
#include "mathfunctions.h"

extern const double g;
extern void* xcalloc(int number, int size);

void SlopeLimiter()
{	
/*	double avgZeta[NumEl], avgQx[NumEl], avgQy[NumEl];

	// Calculate the average value of the variables over each element
	for (int el = 0; el < NumEl; ++el)
	{
		int ind0 = el*3;
		int ind1 = el*3 + 1;
		int ind2 = el*3 + 2;
		double zeta0 = zeta[ind0];
		double zeta1 = zeta[ind1];
		double zeta2 = zeta[ind2];

		double Qx0 = Qx[ind0];
		double Qx1 = Qx[ind1];
		double Qx2 = Qx[ind2];

		double Qy0 = Qy[ind0];
		double Qy1 = Qy[ind1];
		double Qy2 = Qy[ind2];

		avgZeta[el] = (zeta0+zeta1+zeta2)/3;
		avgQx[el] = (Qx0+Qx1+Qx2)/3;
		avgQy[el] = (Qy0+Qy1+Qy2)/3;

	}

*/	// For each node, count the number of elements that share this node
	int *ElCount = xcalloc(NumNodes, sizeof(int));
	//int ElCount[NumNodes];	
	for (int n = 0; n < NumNodes; ++n)
	{
		ElCount[n] = 0;
	}
	for (int el = 0; el<NumEl; ++el)
	{
		for (int n = 0; n < 3; ++n)
		{
			int ind = el*3+n;
			int node = EltoVert[ind];
			ElCount[node] = ElCount[node]+1;
		}
	}

	// In NtoEl, store all the elements that share the node
	// In NtoLocalNode store the local node number of the global node
	// corresponding to the element stored in N to El
	// In ZetaVrtVal, QxVrtVal and QyVrtVal, store the nodal values of 
	// those quantities.
	int **NtoEl = xcalloc(NumNodes,sizeof(int*));
//	int **NtoLocalNode = malloc(NumNodes*sizeof(int*));

	for (int n = 0; n < NumNodes; ++n)
	{	
		int NumConnEl = ElCount[n];
		NtoEl[n] = xcalloc(NumConnEl,sizeof(int));
//		NtoLocalNode[n] = malloc(NumConnEl*sizeof(int));

		for (int i = 0; i < NumConnEl; ++i)
			NtoEl[n][i] = -1;
	}
	
	for (int el =0; el<NumEl; ++el)
	{
		for (int n =0; n < 3; ++n)
		{	
			int ind = el*3+n;
			int node = EltoVert[ind];
			for (int k = 0; k < ElCount[node]; ++k)
			{
				if (NtoEl[node][k] < 0)
				{
					NtoEl[node][k] = el;
		//			NtoLocalNode[node][k] = n;
					break;
				} 
			}
		}
	}

	double *Zeta_Max = xcalloc(NumNodes, sizeof(double));
	double *Qx_Max = xcalloc(NumNodes, sizeof(double));
	double *Qy_Max = xcalloc(NumNodes, sizeof(double));
	double *Zeta_Min = xcalloc(NumNodes, sizeof(double));
	double *Qx_Min = xcalloc(NumNodes, sizeof(double));
	double *Qy_Min = xcalloc(NumNodes, sizeof(double));

	//double Zeta_Max[NumNodes];
	//double Qx_Max[NumNodes];
	//double Qy_Max[NumNodes];
	//double Zeta_Min[NumNodes];
	//double Qx_Min[NumNodes];
	//double Qy_Min[NumNodes];

	// Find the maximum and minimum of the average of each variable over all elements sharing a node
	for (int i =0; i < NumNodes; ++i)
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
			int el = NtoEl[i][j];
 			double avgZeta = (zeta[3*el]+zeta[3*el+1]+zeta[3*el+2])/3;
 			double avgQx = (Qx[3*el]+Qx[3*el+1]+Qx[3*el+2])/3;
 			double avgQy = (Qy[3*el]+Qy[3*el+1]+Qy[3*el+2])/3;

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

			/*
			if (avgZeta[el] < Zeta_Min[i])
				Zeta_Min[i] = avgZeta[el];
			if (avgZeta[el] > Zeta_Max[i])
				Zeta_Max[i] = avgZeta[el];
			if (avgQx[el] < Qx_Min[i])
				Qx_Min[i] = avgQx[el];
			if (avgQx[el] > Qx_Max[i])
				Qx_Max[i] = avgQx[el];
			if (avgQy[el] < Qy_Min[i])
				Qy_Min[i] = avgQy[el];
			if (avgQy[el] > Qy_Max[i])
				Qy_Max[i] = avgQy[el];
			*/
		}

	}	


	// Loop over elements to calculate the new vertex values
	for (int el = 0; el < NumEl; ++el)
	{	
		double ZeVertex[3], ZeMin[3], ZeMax[3];
		double ZeAvg;

		int n1 = EltoVert[el*3];
		int n2 = EltoVert[el*3+1];
		int n3 = EltoVert[el*3+2];

		int modify = 0;

		// loop over variables
		for (int ivar = 0; ivar < 3; ++ivar)
		{

			if (ivar == 0) 
			{
				double Zeta1 = zeta[el*3];
				double Zeta2 = zeta[el*3+1];
				double Zeta3 = zeta[el*3+2];
				ZeVertex[0] = Zeta1 - Zeta2 + Zeta3;
				ZeVertex[1] = Zeta1 + Zeta2 - Zeta3;
				ZeVertex[2] = -Zeta1 + Zeta2 + Zeta3;
				ZeMax[0] = Zeta_Max[n1];
				ZeMax[1] = Zeta_Max[n2];
				ZeMax[2] = Zeta_Max[n3];
				ZeMin[0] = Zeta_Min[n1];
				ZeMin[1] = Zeta_Min[n2];
				ZeMin[2] = Zeta_Min[n3];
				ZeAvg = (Zeta1+Zeta2+Zeta3)/3;
				//ZeAvg = avgZeta[el];
				
				double z1 = z[n1];
				double z2 = z[n2];
				double z3 = z[n3];
				double H1 = ZeVertex[0] + z1;
 				double H2 = ZeVertex[1] + z2;
				double H3 = ZeVertex[2] + z3;

				double Hmin0 = ZeMin[0] + z1;
				double Hmin1 = ZeMin[1] + z2;
				double Hmin2 = ZeMin[2] + z3;
				if (H1 > H0 && H2 > H0 && H3> H0 && Hmin0 > H0 && Hmin1 > H0 && Hmin2 > H0)
					modify = 1;

			}

			else if (ivar == 1)
			{
				double Qx1 = Qx[el*3];
				double Qx2 = Qx[el*3+1];
				double Qx3 = Qx[el*3+2];
				ZeVertex[0] = Qx1 - Qx2 + Qx3;
				ZeVertex[1] = Qx1 + Qx2 - Qx3;
				ZeVertex[2] = -Qx1 + Qx2 + Qx3;
				ZeMax[0] = Qx_Max[n1];
				ZeMax[1] = Qx_Max[n2];
				ZeMax[2] = Qx_Max[n3];
				ZeMin[0] = Qx_Min[n1];
				ZeMin[1] = Qx_Min[n2];
				ZeMin[2] = Qx_Min[n3];
				ZeAvg = (Qx1+Qx2+Qx3)/3;
				//ZeAvg = avgQx[el];

			}		
	
			else
			{
				double Qy1 = Qy[el*3];
				double Qy2 = Qy[el*3+1];
				double Qy3 = Qy[el*3+2];
				ZeVertex[0] = Qy1 - Qy2 + Qy3;
				ZeVertex[1] = Qy1 + Qy2 - Qy3;
				ZeVertex[2] = -Qy1 + Qy2 + Qy3;
				ZeMax[0] = Qy_Max[n1];
				ZeMax[1] = Qy_Max[n2];
				ZeMax[2] = Qy_Max[n3];
				ZeMin[0] = Qy_Min[n1];
				ZeMin[1] = Qy_Min[n2];
				ZeMin[2] = Qy_Min[n3];
				ZeAvg = (Qy1+Qy2+Qy3)/3;
				//ZeAvg = avgQy[el];

			}	

			// Reset the vertex value to be less than or equal to the max
			// and greater than or equal to the min at that vertex
			ZeVertex[0] = fmax(min(ZeVertex[0],ZeMax[0]),ZeMin[0]);
			ZeVertex[1] = fmax(min(ZeVertex[1],ZeMax[1]),ZeMin[1]);
			ZeVertex[2] = fmax(min(ZeVertex[2],ZeMax[2]),ZeMin[2]);


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
					double div = fmax(1, kdp);
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

			double ZeMid1 = (ZeVertex[0]+ZeVertex[1])/2;
			double ZeMid2 = (ZeVertex[1]+ZeVertex[2])/2;
			double ZeMid3 = (ZeVertex[2]+ZeVertex[0])/2;

			if (ivar == 0)
			{
				double H1 = ZeVertex[0] + z[n1];
				double H2 = ZeVertex[1] + z[n2];
				double H3 = ZeVertex[2] + z[n3];
				if (H1 < H0 || H2 < H0 || H3 < H0)
					modify = 0;
				if (modify)
				{
					zeta[el*3] = ZeMid1;
					zeta[el*3+1] = ZeMid2;
					zeta[el*3+2] = ZeMid3;
				}

			}	

			if (ivar == 1)
			{
				if(modify)
				{
					Qx[el*3] = ZeMid1;
					Qx[el*3+1] = ZeMid2;
					Qx[el*3+2] = ZeMid3;
				}
			}	

			if (ivar == 2)
			{
				if(modify)
				{
					Qy[el*3] = ZeMid1;
					Qy[el*3+1] = ZeMid2;
					Qy[el*3+2] = ZeMid3;
				}
			}
	
		} // end loop over variables

	} // end loop over elements	
	
	for (int n = 0; n < NumNodes; ++n)
	{	
		free(NtoEl[n]);
//		free(NtoLocalNode[n]);
	}
	free(NtoEl);
	free(ElCount);
//	free(NtoLocalNode);
	free(Zeta_Max);
	free(Qx_Max);
	free(Zeta_Min);
	free(Qx_Min);
	free(Qy_Max);
	free(Qy_Min);
}
/* Cockburn and Shu Slope Limiter 
void SlopeLimiter(struct junction *junc)
{
	double ZERO = 1e-15;
	int NumEl = junc->NumEl;

	double xb0[NumEl];	// store the x-coordinate of the barycenter of each element
	double yb0[NumEl];	// store the y-coordinate of the barycenter of each element

	int E[NumEl][3];	// for a patch corresponding to element i, store its neighbors
	double xb[NumEl][3];	// store the x-coordinate of the barycenter of neighbors
	double yb[NumEl][3];	// store the y-coordinate of the barycenter of neighbors

	double A0[NumEl];		// store the area of the triangle formed by the vertices of an edge and the centroid	
 
	// The values of cell averages of the conserved variables and the bathymetry
	double avgzeta[NumEl], avgz[NumEl], avgQx[NumEl], avgQy[NumEl];

	double midx[NumEl][3];
	double midy[NumEl][3];
	double z[NumEl][3];				// z at the midpoint of the edges
	double normvecx[NumEl][3];
	double normvecy[NumEl][3];

	double zeta[NumEl][3];				// zeta at the midpoint of the edges
	double Qx[NumEl][3];				// Qx at the midpoint of the edges
	double Qy[NumEl][3];				// Qy at the midpoint of the edges

	for (int i =0; i < NumEl; ++i)
	{
		for (int j = 0; j < 3;  ++j)
		{
			int ind = i*3+j;
			int edg = junc->EltoEdges[ind];
				
			z[i][j] = junc->interpz[edg];
			midx[i][j] = junc->interpx[edg];
			midy[i][j] = junc->interpy[edg];

			zeta[i][j] = junc->zeta[ind];
			Qx[i][j] = junc->Qx[ind];
			Qy[i][j] = junc->Qy[ind];

			normvecx[i][j] = junc->normvecx[ind];
			normvecy[i][j] = junc->normvecy[ind];

			E[i][j] = junc->Neighbors[ind];
			
		}
	}

	// Find the barycenter and the area of the triangle formed by the connecting an edge of the triangle with the barycenter. Also compute the averages of H, Qx and Qy over each element.
	double psi1 = 1-2./3;
	double psi2 = -1+2./3+2./3;
	double psi3 = 1-2./3;

	for (int i =0; i < NumEl; ++i)
	{
		xb0[i] = midx[i][0]*psi1 + midx[i][1]*psi2 + midx[i][2]*psi3;
		yb0[i] = midy[i][0]*psi1 + midy[i][1]*psi2 + midy[i][2]*psi3;
		A0[i] = 1./3*junc->TriangleArea[i];

		// averages of conserved quantities
		avgzeta[i] = zeta[i][0]*psi1 + zeta[i][1]*psi2 + zeta[i][2]*psi3;
		avgQx[i] = Qx[i][0]*psi1 + Qx[i][1]*psi2 + Qx[i][2]*psi3;
		avgQy[i] = Qy[i][0]*psi1 + Qy[i][1]*psi2 + Qy[i][2]*psi3;
		avgz[i] = z[i][0]*psi1 + z[i][1]*psi2 + z[i][2]*psi3;

	}

	// Store the barycenters of the neighbors
	for (int i=0; i<NumEl; ++i)
	{
		for (int j = 0; j<3; ++j)
		{
			// find the element opposite face j of element i
			int el = E[i][j];
			xb[i][j] = xb0[el];
			yb[i][j] = yb0[el];
			
			// if edge j is a boundary edge, then compute the location
			// of the center of the reflected ghost element
			int edg_ind = i*3+j;
			int edg = junc->EltoEdges[edg_ind];
			int el_ind1 = edg*2;
			int el_ind2 = edg*2 + 1;
			int el1 = junc->EdgtoEls[el_ind1];
			int el2 = junc->EdgtoEls[el_ind2];
			if (el1 == el2)
			{
				double H = 2*A0[i]/junc->EdgLength[edg_ind];
				xb[i][j] = xb[i][j] + 2*junc->normvecx[edg_ind]/H;
				yb[i][j] = yb[i][j] + 2*junc->normvecy[edg_ind]/H;
			}

		}
	}

	// Loop over elements
	for (int el = 0; el < NumEl; ++el)
	{	double alpha1[3], alpha2[3], delta_w[3];	// alpha and delta for each edge
		double delta_O[3][3], delta_c[3][3], delta_hat[3][3];	// delta in the original and
							// characterisitic space
		double WC_O[3][3];		// store the averages for the three conserved quantities at the element el, and its two neighbors
		double WM_O[3];			// the values of the conserved quantities at the midpoint of an edge;
		double WC_C[3][3];			// Transform WC_O to characteristic space
		double WM_C[3];			// Transform WM_O to characteristic space 
		double Wtilda[3];
		double LW_O[3][3]; 		// Store limited variables
		

		// loop over the edges to compute alpha1 and alpha2
		for (int i = 0; i < 3; ++i)
		{
			int el_n1 = i;
			int el_n2 = (i+1)%3;
			double A[2][2];
			double B[2];
			A[0][0] = xb[el][el_n1] - xb0[el];
			A[0][1] = xb[el][el_n2] - xb0[el];
			A[1][0] = yb[el][el_n1] - yb0[el];
			A[1][1] = yb[el][el_n2] - yb0[el];
			B[0] = midx[el][i]-xb0[el];
			B[1] = midy[el][i]-yb0[el];
			double dtm = A[0][0]*A[1][1] - A[0][1]*A[1][0];
			alpha1[i] = (A[1][1]*B[0] - A[0][1]*B[1])/dtm;
			alpha2[i] = (-A[1][0]*B[0] + A[0][0]*B[1])/dtm;

		}

		// Set the variable vectors at the barycenter of the element
		WC_O[0][0] = avgzeta[el];
		WC_O[1][0] = avgQx[el];
		WC_O[2][0] = avgQy[el];

		// Start computing deltas for each edge
		for (int i = 0; i < 3; ++i)
		{
			// Compute eigenvalues and eigenvectors at the midpoint using 
			// the variables at the barycenter
			double Ht = avgzeta[el] + avgz[el];
			double nx = normvecx[el][i];
			double ny = normvecy[el][i];
			double u = avgQx[el]/Ht;
			double v = avgQy[el]/Ht;

			// right and left eigenvectors. Left is just the inverse of right in our case
			double RI[3][3], LE[3][3];
			double c = sqrt(g*Ht);
			RI[0][0] = 1.;
			RI[0][1] = 0.;
			RI[0][2] = 1.;
			
			RI[1][0] = u - c*nx;
			RI[1][1] = -ny;
			RI[1][2] = u + c*nx;
			
			RI[2][0] = v - c*ny;
			RI[2][1] = nx;
			RI[2][2] = v + c*ny;

			double dtm = 1./(-RI[0][2]*RI[1][1]*RI[2][0] + RI[0][1]*RI[1][2]*RI[2][0]+RI[0][2]*RI[1][0]*RI[2][1]-RI[0][0]*RI[1][2]*RI[2][1]-RI[0][1]*RI[1][0]*RI[2][2]+RI[0][0]*RI[1][1]*RI[2][2]);
			LE[0][0] = dtm*(-RI[1][2]*RI[2][1]+RI[1][1]*RI[2][2]);
			LE[0][1] = dtm*(RI[0][2]*RI[2][1]-RI[0][1]*RI[2][2]);
			LE[0][2] = dtm*(-RI[0][2]*RI[1][1]+RI[0][1]*RI[1][2]);
		
			LE[1][0] = dtm*(RI[1][2]*RI[2][0]-RI[1][0]*RI[2][2]);
			LE[1][1] = dtm*(-RI[0][2]*RI[2][0]+RI[0][0]*RI[2][2]);
			LE[1][2] = dtm*(RI[0][2]*RI[1][0] - RI[0][0]*RI[1][2]);
		
			LE[2][0] = dtm*(-RI[1][1]*RI[2][0]+RI[1][0]*RI[2][1]);
			LE[2][1] = dtm*(RI[0][1]*RI[2][0]-RI[0][0]*RI[2][1]);
			LE[2][2] = dtm*(-RI[0][1]*RI[1][0]+RI[0][0]*RI[1][1]);

			// Set variable vectors
			WM_O[0] = zeta[el][i];
			WM_O[1] = Qx[el][i];
			WM_O[2] = Qy[el][i];

			int el_n1 = i;
			int el_n2 = (i+1)%3;
			int n1 = E[el][el_n1];
			int n2 = E[el][el_n2];

			WC_O[0][1] = avgzeta[n1];
			WC_O[1][1] = avgQx[n1];
			WC_O[2][1] = avgQy[n1];

			WC_O[0][2] = avgzeta[n2];
			WC_O[1][2] = avgQx[n2];
			WC_O[2][2] = avgQy[n2];

			// take the boundary value to be the average value of the neighboring element if the edge is an inflow or an outflow edge
			int edg1_num = junc->EltoEdges[el*3+el_n1];
			int edg2_num = junc->EltoEdges[el*3+el_n2];
			if (junc->BdryPrescribed[edg1_num])
			{
				WC_O[0][1] = junc->bzeta[edg1_num]; 
				double Qn = junc->bQn[edg1_num];
				double nxe1 = junc->normvecx[el*3+edg1_num];
				double nye1 = junc->normvecy[el*3+edg1_num];
				double txe1 = -nye1;
				double tye1 = nxe1;
				double Qt = Qx[el][i]*txe1+Qy[el][i]*tye1;
				double denom = 1./(nxe1*tye1-nye1*txe1);
				WC_O[1][1] = (tye1*Qn - nye1*Qt)*denom; 
				WC_O[2][1] = (-txe1*Qn + nxe1*Qt)*denom; 
			}

			if (junc->BdryPrescribed[edg2_num])
			{
				WC_O[0][2] = junc->bzeta[edg2_num]; 
				double Qn = junc->bQn[edg2_num];
				double nxe2 = junc->normvecx[el*3+edg2_num];
				double nye2 = junc->normvecy[el*3+edg2_num];
				double txe2 = -nye2;
				double tye2 = nxe2;
				double Qt = Qx[el][i]*txe2+Qy[el][i]*tye2;
				double denom = 1./(nxe2*tye2-nye2*txe2);
				WC_O[1][2] = (tye2*Qn - nye2*Qt)*denom; 
				WC_O[2][2] = (-txe2*Qn + nxe2*Qt)*denom; 
			}

			// Transform original variables into the characteristic space
			for (int j = 0; j < 3; ++j)
			{
				WC_C[j][0] = LE[j][0]*WC_O[0][0]+LE[j][1]*WC_O[1][0]+LE[j][2]*WC_O[2][0];
				WC_C[j][1] = LE[j][0]*WC_O[0][1]+LE[j][1]*WC_O[1][1]+LE[j][2]*WC_O[2][1];
				WC_C[j][2] = LE[j][0]*WC_O[0][2]+LE[j][1]*WC_O[1][2]+LE[j][2]*WC_O[2][2];

				WM_C[j] = LE[j][0]*WM_O[0]+LE[j][1]*WM_O[1]+LE[j][2]*WM_O[2];
			}
			
			// Compute W_tilda and delta_w
			for (int j=0; j<3; ++j)
			{
				Wtilda[j] = WM_C[j] - WC_C[j][0];
				delta_w[j] = alpha1[i]*(WC_C[j][1]-WC_C[j][0]) +
					 alpha2[i]*(WC_C[j][2]-WC_C[j][0]); 
			}

			// Apply the TVB modified minmod function
			double r_sq = (midx[el][i]-xb0[el])*(midx[el][i]-xb0[el]) +
					(midy[el][i]-yb0[el])*(midy[el][i]-yb0[el]);

			// Loop for components of variable vectors
			for (int k=0; k < 3; ++k)
			{
				// take M = 50 and nu = 1.5 like in Cockburn and Shu
				double M = 50;
				double nu = 1.5;
				if (fabs(Wtilda[k]) < M*r_sq)
					delta_c[k][i] = Wtilda[k];
				else
				{
					if ((sign(Wtilda[k]) == sign(delta_w[k])) && 
						(fabs(Wtilda[k]) > 0))
					{
						int s = sign(Wtilda[k]);
						delta_c[k][i] = s*min(fabs(Wtilda[k]),fabs(nu*delta_w[k]));
					}
					else
						delta_c[k][i] = 0;

				}

			} // end loop for components of variable vectors

			// Transform back to the original space
			for (int k =0; k < 3; ++k)
			{
			
				delta_O[k][i] = RI[k][0]*delta_c[0][i] + RI[k][1]*delta_c[1][i] + RI[k][2]*delta_c[2][i];
			
			}

		} // end loop for the edges
		
		// Set the variable vector at midpoint 0. This is used later to judge whether the variables should be limited.
		WM_O[0] = zeta[el][0];
		WM_O[1] = Qx[el][0];
		WM_O[2] = Qy[el][0];

		// Find the possibly limited delta
		double W[3][3];
		for (int k = 0 ; k < 3; ++k)
		{
			W[0][k] = zeta[el][k];
			W[1][k] = Qx[el][k];
			W[2][k] = Qy[el][k];
		}
	
		// loop for the variable components
	
		for (int i = 0; i <3; ++i)
		{
			double pos  = 0;
			double neg = 0;
			for (int j = 0; j<3; ++j)
			{	
				pos += fmax(0,delta_O[i][j]);
				neg += fmax(0,-delta_O[i][j]);
			}
			if ((pos < ZERO) || (neg < ZERO))
			{
				for (int j =0; j < 3; ++j)
					delta_hat[i][j] = delta_O[i][j]; 
					//delta_hat[i][j] = 0;
			}
			else
			{
				double theta_p = min(1, neg/pos);
				double theta_n = min(1, pos/neg);
				for (int j=0; j<3; ++j)
					delta_hat[i][j] = theta_p*fmax(0,delta_O[i][j]) -
								theta_n*fmax(0,-delta_O[i][j]);
			}
			
			//if ((delta_hat[i][0] < (WM_O[i]-WC_O[i][0]-ZERO)) || 
			//	(delta_hat[i][0] > (WM_O[i]-WC_O[i][0] + ZERO)))
		//	{
				for (int j = 0; j < 3; ++j)
				{
					LW_O[i][j] = WC_O[i][0] + delta_hat[i][j];
				}
		//	}
		//	else
		//	{
		//		for (int j = 0; j < 3; ++j)
		//		{
		//			LW_O[i][j] = W[i][j];
		//		}
		//	} 
			

		}
		
		for (int j = 0; j < 3; ++j)
		{
			int ind = el*3+j;
			junc->zeta[ind] = LW_O[0][j];
			junc->Qx[ind] = LW_O[1][j];
			junc->Qy[ind] = LW_O[2][j];

		}

	} // end loop for the elements


} // end function

*/
	
		
		


