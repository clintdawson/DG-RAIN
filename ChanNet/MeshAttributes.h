/*************************************************************************//**
* @mainpage ChanNet
*
* This code simulates flow through a network of open channels.
*
* COMPILING THE CODE
* ==================
* 
* To compile the code, type the following in the directory where the Makefile is.
* -# make
*
* For *wetting and drying purposes*, the code has to be recompiled with wetting
* and drying turned on. So in the directory where the Makefile and the source
* files are, do the following:
* 
* -# make clean
*
* -# make WD=-DWDON
* 
* RUNNING THE CODE:
* ==================
* 
* ./ChanNet pathToChannelNodes.in pathToJunctionMesh.in
* 
* DESCRIPTION OF THE INPUT FILES
* =================================
* 
* - **ChannelNodes.in**: This is a file that contains the grid information associated
* with the channels, including the coordinates, the width at the nodes and the
* Manning's n values. 
*
* - **JunctionMesh.in**: This is a file that contains the gird and connectivity information
* associated with the junctions and their discretizations. It should also have 
* information about the connectivity of the channels.
*
* DESCRIPTION OF FILES THAT NEED TO BE CHANGED FOR EACH RUN
* ==========================================================
* 
* - **main.c**: The parameters H0 and VEL_ZERO are defined here. These might need to be 
* tweaked for different wetting/drying cases.
*
* - **initialize_channels.c**: File where the initial conditions on the channels are 
* specified.
*
* - **initialize_junctions.c**: File where the initial conditions on the junctions
* are specified.
*
* - **boundary_conditions.c**: File where the boundary conditions on the channels are
* specified.
*
* @author Prapti Neupane
* ***************************************************************************/

/**********************************************************************//**
*
* @file MeshAttributes.h
*
* This file contains global variables that are used throughout the 
* software for storing the channel and junction structures.
*
* ***********************************************************************/

#ifndef MESH	

#define MESH

#include "ChannelsAndJunctions.h"

extern int NumChannels;							/**< Total number of channels in the network */
extern int NumJunctions;						/**< Total number of junctions in the network */
extern struct channel **ChannelList;			/**< An array that holds the pointers to the channel structures created
													 for each channel in the network */
extern struct junction **JunctionList;			/**< An array that holds the pointers to the junction structures created
													 for each channel in the network */

#endif


