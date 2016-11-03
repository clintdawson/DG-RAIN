/************************************************************************************************//**
* @file initialize_junctions.c
*
* This file contains code to initialize the junction structures created in create_channel_network.c
* and their member fields.
*
* ***************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathfunctions.h"
#include "ChannelsAndJunctions.h"
#include "MeshAttributes.h"

/*********************************************************************************************//**
* This function allocates necessary space for all the member fields of the junction structure.
* Then it assigns initial values for the quantities zeta, Qx and Qy. It also assigns an intiial
* wet/dry state for each element of the %junction.
*
* ************************************************************************************************/ 

void initialize_junctions()
{

	int i;
	for (i=0; i<NumJunctions; ++i)
	{
		int NumEl = JunctionList[i]->NumEl;
		int NumEdges = JunctionList[i]->TotalNumEdges;

		JunctionList[i]->zeta = malloc(3*NumEl*sizeof(double));
		JunctionList[i]->Qx = malloc(3*NumEl*sizeof(double));
		JunctionList[i]->Qy = malloc(3*NumEl*sizeof(double));
		
		JunctionList[i]->bzeta = malloc(NumEdges*sizeof(double));
		JunctionList[i]->bQn = malloc(NumEdges*sizeof(double));

		#ifdef WDON
		JunctionList[i]->WD = malloc(NumEl*sizeof(int));
		#endif
	
		for (int k=0; k<NumEl; ++k)
		{
			int j;
			for (j=0; j <3; ++j)
			{
				JunctionList[i]->zeta[index(k,j,3)] = 2;
				JunctionList[i]->Qx[index(k,j,3)] = 0;
				JunctionList[i]->Qy[index(k,j,3)] =0;
			}
			
			#ifdef WDON
			JunctionList[i]->WD[k] = -1;
			#endif
		}

	}
}
