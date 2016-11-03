#include <stdio.h>
#include <math.h>
#include "Watershed.h"
#include "Globals.h"
#include "math_functions.h"
#include "constitutive_equations.h"

/**********************************************************************************************//**
* This routine couples the channels with the respective kinematic elements
* *************************************************************************************************/

// This is a one way coupling. Kinematic elements do not get input from channels.
// for now, just assume that there is only one floodplain
void couple_kinEls_with_channels(double time, double dt)
{
	static int Nstep = 0;
	double* totalqL = xcalloc(NumChannels, sizeof(double));
	double totalArea = 0;

// begin coupling
	int NumEl = FloodplainList[0]->NumEl;
	int num = 0;
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl* kinEl = KinematicElList[i];
		if (kinEl->isActive && kinEl->numDownstreamEls == 0)
		{
			//printf("coupling channel with KinEl %d\n", i);
			int connectedChannel = kinEl->connChanNum;
			int connectedChanEl = kinEl->connChanEl;

			int Np = kinEl->Np;
			double A = kinEl->A[Np-1];
			double S0 = kinEl->dz[Np-1];
			double nf = kinEl->NodalnFriction[Np-1];
			double weq = kinEl->weq;
			double H = A/weq;
			double Htil = H - 1e-7;					// current wetting/drying treatment
			if (Htil < 0)
				Htil = 0;

			// Manning's relationship
			double u = sqrt(S0)*pow(Htil, 2.0/3)/nf;
			
			// Chezy's relationship
			//double u = sqrt(S0*Htil)*nf;

			double qL = Htil*u*weq;
			
			totalqL[connectedChannel] += qL;

			
			num++;
		}
	
	}

	for (int i = 0; i < NumChannels; ++i)
	{
		int ChanNumNodes = ChannelList[i]->NumNodes;
		for (int j = 0; j < ChanNumNodes; ++j)
		{
			ChannelList[i]->qL[j] = totalqL[i]/ChannelList[i]->channelLength;
		}
	}
	
	// write out qL every minute
	//if (Nstep % 1 == 0)
	//{
	//	FILE* Qfile;
	//	char fileName[100];
	//	sprintf(fileName, "KinFlowQ.dat");
	//	Qfile = fopen(fileName,"a");
	//	fprintf(Qfile, "%lf \t %lf \n", time/60, totalqL[0]*60);
	//	fclose(Qfile);
	//}

	free(totalqL);
	Nstep++;

}




/**********************************************************************************************//**
* This function couples the channels and junctions with each other by applying the coupling 
* conditions on all of the channels and the junctions
*
* *************************************************************************************************/

// assumes we are using nodal P1 elements in junctions
void couple_channels_with_junctions()
{
	// loop over the channels
	for (int i = 0; i < NumChannels; ++i)
	{

		int NumInflowEdges = ChannelList[i]->NumInflowJunctionEdges;
		double totalInflowQn = 0;
		double totalInflowArea = 0;

		// get channel values at inflow (junction's outflow)
		double A_in = ChannelList[i]->A[1];
		double b_in = ChannelList[i]->NodalB[0];
		double m1_in = ChannelList[i]->Nodalm1[0];
		double m2_in = ChannelList[i]->Nodalm2[0];
		double H_in = getH(A_in, b_in, m1_in, m2_in);
		double Qn_in = ChannelList[i]->Q[1];
		double c_in = sqrt(g*A_in/b_in);
		double u_in = Qn_in/A_in;
		double z_in = ChannelList[i]->NodalZ[0];

		// get channels values at outflow (junction's inflow)
		int  NumNodes = ChannelList[i]->NumNodes;
		double A_out = ChannelList[i]->A[NumNodes];
		double b_out = ChannelList[i]->NodalB[NumNodes-1];
		double m1_out = ChannelList[i]->Nodalm1[NumNodes-1];
		double m2_out = ChannelList[i]->Nodalm2[NumNodes-1];
		double H_out = getH(A_out, b_out, m1_out, m2_out);
		double Qn_out = ChannelList[i]->Q[NumNodes];
		double c_out = sqrt(g*A_out/b_out);
		double u_out = Qn_out/A_out;
		double z_out = ChannelList[i]->NodalZ[NumNodes-1];

		// loop over all the junction edges that are connected to the 
		// beginning of this channel
		for (int j = 0; j < NumInflowEdges; ++j)
		{
			int globalEdgNum = ChannelList[i]->InflowJunctionEdges[2*j];
			int juncNum = ChannelList[i]->InflowJunctionEdges[2*j+1];

			struct TwoDRegion *junc = JunctionList[juncNum];

			int el = junc->EdgtoEls[2*globalEdgNum];	
			int edg = junc->GlobaltoLocalEdg[globalEdgNum*2];
		
			int edgv1 = junc->EdgtoVert[globalEdgNum*2];
			int edgv2 = junc->EdgtoVert[globalEdgNum*2+1];
			double xv1 = junc->Vx[edgv1];
			double xv2 = junc->Vx[edgv2];
			double yv1 = junc->Vy[edgv1];
			double yv2 = junc->Vy[edgv2];
			double edgLength = sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2));

			double nx = junc->nx[el*3+edg];
			double ny = junc->ny[el*3+edg];

			// calculate average value on the edge of the junction
			int n1 = Fmask[0][edg];
			int n2 = Fmask[1][edg];
			double juncZeta = 0.5*(junc->zeta[el][n1] + junc->zeta[el][n2]);
			double juncZ = 0.5*(junc->NodalZ[el][n1] + junc->NodalZ[el][n2]);
			double juncH = juncZeta + juncZ;
			double Qx = 0.5*(junc->Qx[el][n1] + junc->Qx[el][n2]);
			double Qy = 0.5*(junc->Qy[el][n1] + junc->Qy[el][n2]);
		
			double juncQn = Qx*nx+Qy*ny;
	
			totalInflowArea += juncH*edgLength;
			totalInflowQn += juncQn*edgLength;	
		
			JunctionList[juncNum]->bzeta[globalEdgNum] = H_in - juncZ;
			JunctionList[juncNum]->bQn[globalEdgNum] = Qn_in/b_in;

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
			
			struct TwoDRegion *junc = JunctionList[juncNum];

			int el = junc->EdgtoEls[2*globalEdgNum];	
			int edg = junc->GlobaltoLocalEdg[globalEdgNum*2];
		
			int edgv1 = junc->EdgtoVert[globalEdgNum*2];
			int edgv2 = junc->EdgtoVert[globalEdgNum*2+1];
			double xv1 = junc->Vx[edgv1];
			double xv2 = junc->Vx[edgv2];
			double yv1 = junc->Vy[edgv1];
			double yv2 = junc->Vy[edgv2];
			double edgLength = sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2));

			double nx = junc->nx[el*3+edg];
			double ny = junc->ny[el*3+edg];

			// calculate average value on the edge of the junction
			int n1 = Fmask[0][edg];
			int n2 = Fmask[1][edg];
			double juncZeta = 0.5*(junc->zeta[el][n1] + junc->zeta[el][n2]);
			double juncZ = 0.5*(junc->NodalZ[el][n1] + junc->NodalZ[el][n2]);
			double juncH = juncZeta + juncZ;
			double Qx = 0.5*(junc->Qx[el][n1] + junc->Qx[el][n2]);
			double Qy = 0.5*(junc->Qy[el][n1] + junc->Qy[el][n2]);
	
			//if (i == 2)
			//	printf("Qx = %lf Qy = %lf\n", Qx, Qy);

			double juncQn = Qx*nx+Qy*ny;
		
			totalOutflowArea += juncH*edgLength;
			totalOutflowQn += juncQn*edgLength;	

			
			JunctionList[juncNum]->bzeta[globalEdgNum] = H_out - juncZ;
			JunctionList[juncNum]->bQn[globalEdgNum] = -Qn_out/b_out;
				
		} // end outflow junction edges loop

		if (NumOutflowEdges != 0)
		{
		
			int NumNodes = ChannelList[i]->NumNodes;
		
			ChannelList[i]->Q[NumNodes+1] = -totalOutflowQn;
			ChannelList[i]->A[NumNodes+1] = totalOutflowArea;

			//if (i == 2)
			//{
			//	printf("juncA = %lf\n", totalOutflowArea);
			//	printf("juncQ = %lf\n", -totalOutflowQn);
			//}
		}

	} // end channel loop

}

// in the array chval, this function returns the water height and normal flow in the channels connected to each of the edges
// because widths of channels are supposed to be narrower than the resolution of a 2D mesh, we don't allow the possibility that
// one channel might be connected to multiple shallow water element edges
void obtain_channels_values(int *ncheds, int chnums[][2], double chvals[][2])
{
	for (int i = 0; i < (*ncheds); ++i)
	{
		int chanNum = chnums[i][0];
		int connType = chnums[i][1];

		struct channel *myChan = ChannelList[chanNum];

		int chanNumEdges = myChan->NumEdges;
		double b, m1, m2, zval, uval;
		double A;
		double Qn;
		if (connType == 0)
		{
			A = myChan->A[1];
			Qn = myChan->Q[1];
			b = myChan->b[0];
			m1 = myChan->m1[0];
			m2 = myChan->m2[0];
			zval = myChan->z[0];

		}
		else
		{
			int NumNodes = myChan->NumNodes;
			A = myChan->A[NumNodes];
			Qn = myChan->Q[NumNodes];
			b = myChan->b[chanNumEdges-1];
			m1 = myChan->m1[chanNumEdges-1];
			m2 = myChan->m2[chanNumEdges-1];
			Qn = -Qn;
			zval = myChan->z[chanNumEdges-1];
		}
	
		double ht = getH(A, b, m1, m2);
		double zeta = ht - zval;

		chvals[i][0] = zeta;
		chvals[i][1] = Qn;

		//FILE *Sendfile = fopen("QatEnd.dat", "a");
		//fprintf(Sendfile, "%3.16f \t %3.16f\n", -chvals[i][1]);
		//fclose(Sendfile);

		//printf("values sent to dgswem: %3.16f %3.16f \n", chvals[i][0], -chvals[i][1]);
	}
}

void receive_swe_values(int *ncheds, int chnums[][2], double swevals[][2])
{
	for (int i = 0; i < (*ncheds); ++i)
	{
		
		int chanNum = chnums[i][0];
		int connType = chnums[i][1];

		struct channel* myChan = ChannelList[chanNum];
		double zeta = swevals[i][0];
		double Qn = swevals[i][1];

		if (connType == 0)
		{
			double zval = myChan->z[0];
			double bval = myChan->b[0];
			double m1val = myChan->m1[0];
			double m2val = myChan->m2[0];
			double htval = zeta + zval;
			double Area = htval*bval + 0.5*(m1val+m2val)*htval*htval;
			myChan->A[0] = Area;
			myChan->Q[0] = Qn;
		}
		else
		{
			int NumNodes = myChan->NumNodes;
			double zval = myChan->NodalZ[NumNodes-1];
			double bval = myChan->NodalB[NumNodes-1];
			double m1val = myChan->Nodalm1[NumNodes-1];
			double m2val = myChan->Nodalm2[NumNodes-1];
			double htval = zeta+zval;
			double Area = htval*bval + 0.5*(m1val+m2val)*htval*htval;
			myChan->A[NumNodes+1] = Area;
			myChan->Q[NumNodes+1] = -Qn;
		}
		//printf("received values from dgswem: %4.16f %3.16f \n", zeta, -Qn);
	}
}

