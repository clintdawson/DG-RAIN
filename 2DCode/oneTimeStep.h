#ifndef ONE_TIME_STEP

#define ONE_TIME_STEP

extern void compute2DInnerProducts(double SI[][9], double VI[][21], double time);
extern void compute2DL(double time);
extern void SlopeLimiter();
extern void boundary_conditions();

#endif
