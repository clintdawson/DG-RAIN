#ifndef MESH	

#define MESH

extern int NumEl;
extern int NumNodes;
extern int NumEdges;
extern double *x_coord;
extern double *y_coord;
extern double *z;
extern int *EltoVert;
extern int *EdgtoEls;
extern int *EdgtoVert;
extern int *BdryPrescribed;		// 1 for inflow boundary. 2 for outflow boundary. 0 for no flow.

//extern int *InflowBdryNodes1;
//extern int *InflowBdryNodes2;
//extern int *OutflowBdryNodes;

extern double minEdgLength;
extern double max_lambda;

extern double *zeta;
extern double *Qx;
extern double *Qy;

extern double *Fhat1dotn;
extern double *Fhat2dotn;
extern double *Fhat3dotn;
extern double *RHSZeta;
extern double *RHSQx;
extern double *RHSQy;

extern double *bzeta;
extern double *bQn;

extern const double g;
extern double *nFriction;

extern const double H0;
extern const double VELZERO;
extern int *WD;


#endif
