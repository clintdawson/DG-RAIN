/********************************************************************************//**
* @file main.c
*
* This file is the driver for ChanNet. When wetting and drying is on, the 
* minimum water threshold and the momentum threshold are also allocated in this file. 
* This needs to be changed later so that it will be a parameter than can be read 
* from a file.
************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ChannelsAndJunctions.h"
#include "SimulationSteps.h"

const double g = 9.810000000000;	/**< gravitational acceleration constant */

/* @cond GLOBAL_VARS_DECLARATION */
int NumChannels;			// Total Number of Channels in the mesh
int NumJunctions;			// Total Number of Junctions in the mesh
struct channel **ChannelList;
struct junction **JunctionList;

/* @endcond */

#ifdef WDON
const double H0 = 1e-3;     		/**< Minimum water height threshold that is maintained in all 
										 elements. Nodes whose water height is below this threshold 
										 are considered dry. */
const double VELZERO = 1e-2;		/**< a momentum threshold. If the magnitude of the momentum 
								  		at a dry node is below this threshold value, then this
								  		momentum will not be transferred to another wet node in
								  		the element.*/

#endif

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		printf("Usage: ./ChanNet 'ChannelNodes.in' 'JunctionNodes.in'\n");
		exit(EXIT_FAILURE);

	}
	
	char *ChannelNodes = argv[1];
	char *JunctionNodes = argv[2];

	//create junctions and channels from the mesh
	create_channel_network(ChannelNodes, JunctionNodes);

	double FinalTime;
	printf("Enter the time you would like to run the simulation till:\n");
	scanf("%lf", &FinalTime);
	
	// initialize
	initialize_channels();
	initialize_junctions();	

	/* Step through time */

	time_evolution(FinalTime);

	return(0);

}


 
	
	
