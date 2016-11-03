#include "ChannelsAndJunctions.h"
#include "MeshAttributes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathfunctions.h"

extern const double g;

/**********************************************************************************************//**
* This function couples the channels and junctions with each other by applying the coupling 
* conditions on all of the channels and the junctions
*
* *************************************************************************************************/

void internal_BC()
{
	// loop over the channels
	for (int i = 0; i < NumChannels; ++i)
	{
	
		int NumInflowEdges = ChannelList[i]->NumInflowJunctionEdges;
		double totalInflowQn = 0;
		double totalInflowArea = 0;
		int NumBoundCondition;

		// loop over all the junction edges that are connected to the 
		// beginning of this channel
		for (int j = 0; j < NumInflowEdges; ++j)
		{
			int globalEdgNum = ChannelList[i]->InflowJunctionEdges[2*j];
			int juncNum = ChannelList[i]->InflowJunctionEdges[2*j+1];
			
			struct junction *junc = JunctionList[juncNum];
			int el = junc->EdgtoEls[2*globalEdgNum];	
			
			int edg;
			int v1 = junc->EltoVert[el*3];
			int v2 = junc->EltoVert[el*3+1];
			int v3 = junc->EltoVert[el*3+2];
			double xv1 = junc->x[v1];
			double xv2 = junc->x[v2];
			double xv3 = junc->x[v3];
			double yv1 = junc->y[v1];
			double yv2 = junc->y[v2];
			double yv3 = junc->y[v3];
			double edgLength;

			int edgv1 = junc->EdgtoVert[globalEdgNum*2];
			int edgv2 = junc->EdgtoVert[globalEdgNum*2+1];
			
			double x_coord[2], y_coord[2];

			if ((v1 == edgv1 && v2 == edgv2) || (v1 == edgv2 && v2 == edgv1))
			{
				edg = 0;
				x_coord[0] = xv1;
				x_coord[1] = xv2;
				y_coord[0] = yv1;
				y_coord[1] = yv2;
				edgLength = sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2));
			}
			else if ((v2 == edgv1 && v3 == edgv2) || (v2 == edgv2 && v3 == edgv1))
			{
				edg = 1;
				x_coord[0] = xv2;
				x_coord[1] = xv3;
				y_coord[0] = yv2;
				y_coord[1] = yv3;

				edgLength = sqrt((xv2-xv3)*(xv2-xv3)+(yv2-yv3)*(yv2-yv3));
			}		
			else if ((v3 == edgv1 && v1 == edgv2) || (v3 == edgv2 && v1 == edgv1))
			{
				edg = 2;
				x_coord[0] = xv3;
				x_coord[1] = xv1;
				y_coord[0] = yv3;
				y_coord[1] = yv1;

				edgLength = sqrt((xv1-xv3)*(xv1-xv3)+(yv1-yv3)*(yv1-yv3));
			}		

			double juncZeta = junc->zeta[el*3+edg];
			double Qx = junc->Qx[el*3+edg];
			double Qy = junc->Qy[el*3+edg];
			double juncZ = 0.5*(junc->z[edgv1] + junc->z[edgv2]);
			double juncH = juncZeta + juncZ;

			double jac = (xv2-xv1)*(yv3-yv1) - (xv3-xv1)*(yv2-yv1);

			double nx = sign(jac)*(y_coord[1] - y_coord[0])/edgLength;
			double ny = sign(jac)*(x_coord[0] - x_coord[1])/edgLength;

			double juncQn = Qx*nx+Qy*ny;
			
			totalInflowArea += edgLength*juncH;
			totalInflowQn += juncQn*edgLength;	
			
			// apply boundary condition to the junction
			double A = ChannelList[i]->A[1];
			double b = ChannelList[i]->b[0];
			double H = A/(NumInflowEdges*edgLength);
			double Qn = ChannelList[i]->Q[1];
		
			JunctionList[juncNum]->bzeta[globalEdgNum] = H - juncZ;
			JunctionList[juncNum]->bQn[globalEdgNum] = Qn/(NumInflowEdges*edgLength);
			
		} // end inflow junction edges loop
	
		if (NumInflowEdges != 0)
		{
			ChannelList[i]->Q[0] = totalInflowQn;
			ChannelList[i]->A[0] = totalInflowArea;
			
		}

		int NumOutflowEdges = ChannelList[i]->NumOutflowJunctionEdges;
		double totalOutflowQn = 0;
		double totalOutflowArea = 0;

		// loop over all the edges of the junction that are connected to the end
		// of this channel
		for (int j = 0; j < NumOutflowEdges; ++j)
		{
			int globalEdgNum = ChannelList[i]->OutflowJunctionEdges[2*j];
			int juncNum = ChannelList[i]->OutflowJunctionEdges[2*j+1];
			
			struct junction *junc = JunctionList[juncNum];
			int el = junc->EdgtoEls[2*globalEdgNum];	
			
			int edg;
			int v1 = junc->EltoVert[el*3];
			int v2 = junc->EltoVert[el*3+1];
			int v3 = junc->EltoVert[el*3+2];
			double xv1 = junc->x[v1];
			double xv2 = junc->x[v2];
			double xv3 = junc->x[v3];
			double yv1 = junc->y[v1];
			double yv2 = junc->y[v2];
			double yv3 = junc->y[v3];
			double edgLength;

			int edgv1 = junc->EdgtoVert[globalEdgNum*2];
			int edgv2 = junc->EdgtoVert[globalEdgNum*2+1];

			double x_coord[2], y_coord[2];

			if ((v1 == edgv1 && v2 == edgv2) || (v1 == edgv2 && v2 == edgv1))
			{
				edg = 0;
				x_coord[0] = xv1;
				x_coord[1] = xv2;
				y_coord[0] = yv1;
				y_coord[1] = yv2;
				edgLength = sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2));
			}
			else if ((v2 == edgv1 && v3 == edgv2) || (v2 == edgv2 && v3 == edgv1))
			{
				edg = 1;
				x_coord[0] = xv2;
				x_coord[1] = xv3;
				y_coord[0] = yv2;
				y_coord[1] = yv3;

				edgLength = sqrt((xv2-xv3)*(xv2-xv3)+(yv2-yv3)*(yv2-yv3));
			}		
			else if ((v3 == edgv1 && v1 == edgv2) || (v3 == edgv2 && v1 == edgv1))
			{
				edg = 2;
				x_coord[0] = xv3;
				x_coord[1] = xv1;
				y_coord[0] = yv3;
				y_coord[1] = yv1;

				edgLength = sqrt((xv1-xv3)*(xv1-xv3)+(yv1-yv3)*(yv1-yv3));
			}		

			double juncZeta = junc->zeta[el*3+edg];
			double Qx = junc->Qx[el*3+edg];
			double Qy = junc->Qy[el*3+edg];
			double juncZ = 0.5*(junc->z[edgv1] + junc->z[edgv2]);
			double juncH = juncZeta + juncZ;

			double jac = (xv2-xv1)*(yv3-yv1) - (xv3-xv1)*(yv2-yv1);

			double nx = sign(jac)*(y_coord[1] - y_coord[0])/edgLength;
			double ny = sign(jac)*(x_coord[0] - x_coord[1])/edgLength;

			double juncQn = Qx*nx+Qy*ny;
			
			totalOutflowArea += juncH*edgLength;
			totalOutflowQn += juncQn*edgLength;	

			// apply boundary condition to the junction
			int NumNodes = ChannelList[i]->NumNodes;
			double A = ChannelList[i]->A[2*NumNodes-2];
			double b = ChannelList[i]->b[NumNodes-1];
			double H = A/(NumOutflowEdges*edgLength);
			double Qn = ChannelList[i]->Q[2*NumNodes-2];
				
			JunctionList[juncNum]->bzeta[globalEdgNum] = H - juncZ;
			JunctionList[juncNum]->bQn[globalEdgNum] = -Qn/(NumOutflowEdges*edgLength);

		} // end outflow junction edges loop

		if (NumOutflowEdges != 0)
		{
		
			int NumNodes = ChannelList[i]->NumNodes;
	
			ChannelList[i]->Q[2*NumNodes-1] = -totalOutflowQn;
			//ChannelList[i]->Q[2*NumNodes-1] = ChannelList[i]->Q[2*NumNodes-2];
			ChannelList[i]->A[2*NumNodes-1] = totalOutflowArea;
		}

	} // end channel loop
} // end function

/*
double getBeta(struct channel *Chan, int edge, int channelNumber,double time)
{
	double beta = 1;

	int juncNode1 = Chan->JunctionNode[0];
	int juncNode2 = Chan->JunctionNode[1];

	int connJunction = -1;
	int goesIn = 0;
	int goesOut = 0;

	for (int i = 0; i < NumJunctions; ++i)
	{
		int juncNode = JunctionList[i]->AssocNode;
		if (juncNode == juncNode1)
		{
			connJunction = i;
			goesOut = 1;
			break;

		}

		if (juncNode == juncNode2)
		{
			connJunction = i;
			goesIn = 1;
			break;
		}

	}

	int NumEl = Chan->NumEl;
	
	if ((goesIn && edge >= NumEl - 10) || (goesOut && edge < 10))
	//if (goesIn || goesOut)
	{
		double A_L = Chan->A[2*edge];
		double A_R = Chan->A[2*edge+1];
		double Q_L = Chan->Q[2*edge];
		double Q_R = Chan->Q[2*edge+1];
		double u_L = Q_L/A_L;
		double u_R = Q_R/A_R;

		double Ahat = 0.5*(A_L+A_R);
		double uhat = (Q_L/sqrt(A_L) + Q_R/sqrt(A_R))/(sqrt(A_L)+sqrt(A_R));
		
		double Qhat = Ahat*uhat;

		// find out which element is connected to the channel
		int elNum = -1;
		struct junction *junc = JunctionList[connJunction];
		int NumJuncEl = junc->NumEl;
		for (int j = 0; j < NumJuncEl; ++j)
		{
			if (junc->AssocChannels[j] == channelNumber)
			{
				elNum = j;
				break;
			}
		}

		// find the value of the conserved variables at the vertices that make up the
		// edge connected to the channel. Note that these will be vertices 1 and 2 of the element
		// because vertex 0 is always the junctionNode

		double zeta1 = junc->zeta[3*elNum];
		double zeta2 = junc->zeta[3*elNum+1];
		double zeta3 = junc->zeta[3*elNum+2];
		
		double Qx1 = junc->Qx[3*elNum];
		double Qx2 = junc->Qx[3*elNum+1];
		double Qx3 = junc->Qx[3*elNum+2];
	
		double Qy1 = junc->Qy[3*elNum];
		double Qy2 = junc->Qy[3*elNum+1];
		double Qy3 = junc->Qy[3*elNum+2];

		double vrt1_zeta = zeta1 + zeta2 - zeta3;
		double vrt2_zeta = -zeta1 + zeta2 + zeta3;

		double vrt1_Qx = Qx1 + Qx2 - Qx3;
		double vrt2_Qx = -Qx1 + Qx2 + Qx3;

		double vrt1_Qy = Qy1 + Qy2 - Qy3;
		double vrt2_Qy = -Qy1 + Qy2 + Qy3;

		// Find the value of u and v at the quadrature points
		int v0 = junc->EltoVert[3*elNum];
		int v1 = junc->EltoVert[3*elNum+1];
		int v2 = junc->EltoVert[3*elNum+2];

		double vrt1_z = junc->z[junc->EltoVert[v1]];
		double vrt2_z = junc->z[junc->EltoVert[v2]];

		double vrt1_h = vrt1_zeta + vrt1_z;
		double vrt2_h = vrt2_zeta - vrt2_z;


		double quad1_h = 0.5*(vrt1_h*(1+1/sqrt(3))+vrt2_h*(1-1/sqrt(3)));
		double quad2_h = 0.5*(vrt1_h*(1-1/sqrt(3))+vrt2_h*(1+1/sqrt(3)));

		double quad1_Qx = 0.5*(vrt1_Qx*(1+1/sqrt(3))+vrt2_Qx*(1-1/sqrt(3)));
		double quad2_Qx = 0.5*(vrt1_Qx*(1-1/sqrt(3))+vrt2_Qx*(1+1/sqrt(3)));

		double quad1_Qy = 0.5*(vrt1_Qy*(1+1/sqrt(3))+vrt2_Qy*(1-1/sqrt(3)));
		double quad2_Qy = 0.5*(vrt1_Qy*(1-1/sqrt(3))+vrt2_Qy*(1+1/sqrt(3)));

		double quad1_u = quad1_Qx/quad1_h;
		double quad2_u = quad2_Qx/quad2_h;
		double quad1_v = quad1_Qy/quad1_h;
		double quad2_v = quad2_Qy/quad2_h;

		// Find the normal vector at the edge
		double xv0 = junc->x[v0];
		double xv1 = junc->x[v1];
		double xv2 = junc->x[v2];
		double yv0 = junc->y[v0];
		double yv1 = junc->y[v1];
		double yv2 = junc->y[v2];

		double jac = (xv1-xv0)*(yv2-yv0) - (xv2-xv0)*(yv1-yv0);
		double edgLength = sqrt((xv2-xv1)*(xv2-xv1)+(yv2-yv1)*(yv2-yv1));
		double nx = sign(jac)*(yv2 - yv1)/edgLength;
		double ny = sign(jac)*(xv1 - xv2)/edgLength;
	
		// Calculate the integral of v_L^2 = v_n^2 = (u*nx+v*ny)^2
		double vnSqIntegral = 0.5*edgLength*(vrt2_h - vrt1_h)*(quad1_u*quad1_u*nx*nx + 
			2*quad1_u*quad1_v*nx*ny + quad1_v*quad1_v*ny*ny + quad2_u*quad2_u*nx*nx +
			2*quad2_u*quad2_v*nx*ny + quad2_v*quad2_v*ny*ny);

		if ( fabs(Qhat) > 1e-10 ) 
			beta = vnSqIntegral/Qhat;
		
	}
	Chan->beta[edge] = beta;
	return beta;
}
*/
