
/*******************************************************************************//**
* @file boundary_conditions.c
*
* Specify the boundary conditions on the inflow and outflow channel ends in 
* this file. Currently these conditions need to be compiled. But eventually
* we might have to change the way these conditions are provided so that they
* can be read from a time series data file and won't hvave to be compiled.
* 
* *********************************************************************************/	

#include <math.h>
#include "ChannelsAndJunctions.h"
#include "MeshAttributes.h"

void boundary_conditions()
{
/*	double Q = 1;
	ChannelList[0]->Q[0] = Q;
	double zeta0 = 5;
	double b0 = ChannelList[0]->b[0];
	double z0 = ChannelList[0]->z[0];
	double h0 = zeta0 + z0;
//	ChannelList[0]->A[0] = h0*b0;
	ChannelList[0]->A[0] = ChannelList[0]->A[1];

	ChannelList[1]->Q[0] = Q;
	double zeta1 = 5;
	double b1 = ChannelList[1]->b[0];
	double z1 = ChannelList[1]->z[0];
	double h1 = zeta1 + z1;
//	ChannelList[1]->A[0] = h1*b1;
	ChannelList[1]->A[0] = ChannelList[1]->A[1];	


	int Chan1NumNodes = ChannelList[1]->NumNodes;
	int Chan1LastNode = 2*Chan1NumNodes-1;
	double zeta1 = 5;
	double b1 = ChannelList[1]->b[Chan1NumNodes-1];
	double z1 = ChannelList[1]->z[Chan1NumNodes-1];
	double h1 = zeta1 + z1;
//	ChannelList[1]->Q[Chan1LastNode] = 2*Q;
	ChannelList[1]->Q[Chan1LastNode] = ChannelList[1]->Q[Chan1LastNode-1];
//	ChannelList[1]->A[Chan1LastNode] = ChannelList[1]->A[Chan1LastNode-1];
	ChannelList[1]->A[Chan1LastNode] = h1*b1;	

	int Chan2NumNodes = ChannelList[2]->NumNodes;
	int Chan2LastNode = 2*Chan2NumNodes-1;
	double zeta2 = 5;
	double b2 = ChannelList[2]->b[Chan2NumNodes-1];
	double z2 = ChannelList[2]->z[Chan2NumNodes-1];
	double h2 = zeta2 + z2;
//	ChannelList[2]->Q[Chan2LastNode] = 2*Q;
	ChannelList[2]->Q[Chan2LastNode] = ChannelList[2]->Q[Chan2LastNode-1];
//	ChannelList[2]->A[Chan2LastNode] = ChannelList[2]->A[Chan2LastNode-1];
	ChannelList[2]->A[Chan2LastNode] = h2*b2;

	ChannelList[2]->Q[0] = Q;
	double zeta2 = 5;
	double b2 = ChannelList[2]->b[0];
	double z2 = ChannelList[2]->z[0];
	double h2 = zeta2 + z2;
	ChannelList[2]->A[0] = h2*b2;
*/

	// T Mesh boundary conditions

/*	ChannelList[0]->Q[0] = 0.0425;
	ChannelList[0]->A[0] = ChannelList[0]->A[1];

	int Chan1LastNode = 2*ChannelList[1]->NumNodes-1;
	ChannelList[1]->Q[Chan1LastNode] = ChannelList[1]->Q[Chan1LastNode-1];
	double Chan1width = ChannelList[1]->b[ChannelList[1]->NumNodes-1];
	ChannelList[1]->A[Chan1LastNode] = 0.296*Chan1width;

	ChannelList[2]->Q[0] = 0.1275;
	ChannelList[2]->A[0] = ChannelList[2]->A[1];
*/	

/*	// Angled channels boundary conditions

	ChannelList[0]->Q[0] = 30;
	ChannelList[0]->A[0] = ChannelList[0]->A[1];

	ChannelList[1]->Q[0] = 20;
	ChannelList[1]->A[0] = ChannelList[1]->A[1];

	int Chan2LastNode = 2*ChannelList[2]->NumNodes-1;
	ChannelList[2]->Q[Chan2LastNode] = ChannelList[2]->Q[Chan2LastNode-1];
	double Chan2width = ChannelList[2]->b[ChannelList[2]->NumNodes-1];
	ChannelList[2]->A[Chan2LastNode] = 1.69*Chan2width;
*/

	// Conservation Example
	ChannelList[0]->Q[0] = 10;
	ChannelList[0]->A[0] = ChannelList[0]->A[1];
	//double height = 2.0;
	//ChannelList[0]->A[0] = ChannelList[0]->b[0]*height;

	ChannelList[1]->Q[0] = 10;
	ChannelList[1]->A[0] = ChannelList[1]->A[1];
	//ChannelList[1]->A[0] = ChannelList[1]->b[0]*height;

	int Chan2LastNode = 2*ChannelList[2]->NumNodes-1;
	ChannelList[2]->Q[Chan2LastNode] = ChannelList[2]->Q[Chan2LastNode-1];
	ChannelList[2]->A[Chan2LastNode] = ChannelList[2]->A[Chan2LastNode-1];



}


