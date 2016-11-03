#ifndef GLOBALS

#define GLOBALS

#include <gsl/gsl_matrix.h>

extern int NumEl;
extern int NumEdges;
extern int NumNodes;
extern int P;
extern int Np;
extern gsl_matrix *LIFT;
extern gsl_matrix *VolMat;
extern gsl_matrix *MassMatrix;

extern double *x;
extern double *y;
extern double *z;
extern double *b;
extern double *m1;
extern double *m2;
extern double *db;
extern double *dm1;
extern double *dm2;
extern double *dz;
extern double *dh;
extern double *nFriction;

extern double *NodalX;
extern double *NodalY;
extern double *NodalB;
extern double *Nodalz;
extern double *Nodalm1;
extern double *Nodalm2;
extern double *NodalnFriction;


extern double *A;
extern double *Q;

//extern double *Fhat1L;
//extern double *Fhat1R;
//extern double *Fhat2L;
//extern double *Fhat2R;
extern double *RHSA;
extern double *RHSQ;

extern double dt;
extern double max_lambda;
extern double mindh;

extern const double g;
extern const double H0;
extern const double VELZERO;

extern int *WD;		// keeps track of the wet/dry status of an element

#endif 
