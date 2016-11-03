/************************************************************************************************//**
* @file initialize_channels.c
*
* This file contains code to initialize the channel structures created in create_channel_network.c
* and their member fields.
*
* *************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathfunctions.h"
#include "ChannelsAndJunctions.h"
#include "MeshAttributes.h"

/***********************************************************************************************//**
* This function allocates necessary space for all the member fields of the channel structure. Then
* it assigns initial values for the quantities A and Q. It also assigns an initial wet/dry state
* for each element of the %channel.
* 
* **************************************************************************************************/ 

void initialize_channels()
{	
	for (int i=0; i<NumChannels; ++i)
	{
		int NumNodes = ChannelList[i]->NumNodes;
		
		ChannelList[i]->A = malloc(2*NumNodes*sizeof(double));
		ChannelList[i]->Q = malloc(2*NumNodes*sizeof(double));

		ChannelList[i]->beta = malloc(NumNodes*sizeof(double));
	
		#ifdef WDON
		ChannelList[i]->WD = malloc(ChannelList[i]->NumEl*sizeof(int));
		#endif

		for (int k=0; k<NumNodes; ++k)
		{
			for (int j=0; j<2; ++j)
			{
				double z_node = ChannelList[i]->z[k];
				double b_node = ChannelList[i]->b[k];
				double zeta = 2;
				double height = zeta + z_node;
				if (b_node == 0)
				{
					printf("width zero at node %d \n", k);
					exit(EXIT_FAILURE);
				}
				ChannelList[i]->A[index(k,j,2)] = height*b_node;
				ChannelList[i]->Q[index(k,j,2)] = 0;
			}

			#ifdef WDON
			if (k < NumNodes-1)
				ChannelList[i]->WD[k] = -1;
			#endif
		}

		
	}

}

