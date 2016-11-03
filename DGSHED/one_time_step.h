/************************************************************************************//**
*
* @file oneTimeStep.h
* 
* This file contains function prototypes for the functions that are exectued in order
* to step forward in time.
* ***************************************************************************************/ 


#ifndef ONE_TIME_STEP

#define ONE_TIME_STEP

#include "Watershed.h"

/* @cond FUNCTION_PROTOTPYES */
extern void computeLKinematicEls();
extern void minmodKinematicEls();

extern void computeChanL(struct channel* Chan, double time, double* dt, int channelNumber, int stage, double* RHSA, double* RHSQ);
extern void minmodChan(struct channel* Chan);
extern void minmodHeightChan(struct channel* Chan);
#ifdef WDON
extern void wetDryStatusChan(struct channel *Chan);
extern void PDopChan(struct channel *Chan);
#endif

extern void compute2DL(struct TwoDRegion* region, double time, double* RHSZeta, double* RHSQx, double* RHSQy, double dt);
extern void SlopeLimiter2D(struct TwoDRegion* region);
extern void SlopeLimiterHeight2D(struct TwoDRegion* region);
#ifdef WDON
extern void wetDryStatus2D(struct TwoDRegion* region);
extern void PDop2D(struct TwoDRegion* region, double time);
#endif

extern void couple_kinEls_with_channels(double time, double dt);
extern void couple_channels_with_junctions();

extern void boundary_conditions(double time);


extern void outputFloodplainDataToFile(double time, double finalTime, double recordTimeIntervals);
extern void outputFloodplainNodalDataToFile(double time, double finalTime, double recordTimeIntervals);
extern void outputDataToFile(double time, double finalTime, double recordTimeIntervals);
extern double calculateTotalWater(struct channel *Chan);
extern void append_to_file(char* fileName, double time, double data);
extern void output2DNodalError(double time);
extern void calculate_l2_error_2d(double time);


/* @uncond */
#endif
