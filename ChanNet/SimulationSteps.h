/**************************************************************************//**
*
* @file SimulationSteps.h
*
* This file contains the function prototypes for functions that represent 
* different stages of the simulation, i.e. function to create the channels
* and junction structure, functions to initialize those structures and the
* function to evolve them in time.
*
* ***************************************************************************/

#ifndef SIMULATION_STEPS

#define SIMULATION_STEPS

#include "ChannelsAndJunctions.h"

/* @cond FUNCTION_PROTOTYPES */
extern void create_channel_network(char *ChannelNodes, char *JunctionNodes);
extern void initialize_channels();
extern void initialize_junctions();
extern void time_evolution(double FinalTime);
/* @endcond */

#endif
