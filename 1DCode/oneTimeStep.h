#ifndef ONE_TIME_STEP

#define ONE_TIME_STEP

extern void wetDryStatus();
extern void PDop();
extern void minmod();
extern void SlopeLimiter1D();
extern void SurfaceSlopeLimiter1D();
extern void computeL();
extern void boundary_conditions();

#endif
