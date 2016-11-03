/************************************************************************************//**
*
* @file oneTimeStep.h
* 
* This file contains function prototypes for the functions that are exectued in order
* to step forward in time.
* ***************************************************************************************/ 


#ifndef ONE_TIME_STEP

#define ONE_TIME_STEP

#include "ChannelsAndJunctions.h"

/* @cond FUNCTION_PROTOTPYES */
extern void minmod(struct channel *Chan);
extern void computeL(struct channel *Chan, double time, int channelNumber, double *RHSA, double *RHSQ);
extern void wetDryStatus1D(struct channel *Chan);
extern void PDop1D(struct channel *Chan);
extern void compute2DL(struct junction *junc, double time, double *RHSZeta, double *RHSQx, double *RHSQy);
extern void SlopeLimiter(struct junction *junc);
extern void wetDryStatus2D(struct junction *junc);
extern void PDop2D(struct junction *junc);
extern void boundary_conditions();
extern void internal_BC();

/* @uncond */
#endif
