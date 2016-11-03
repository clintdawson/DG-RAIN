#include <stdio.h>
#include <stdlib.h>
#include <math.h>	
#include "mathfunctions.h"
#include "globals.h"

int NumEl;
int NumNodes;
int NumEdges;
double *x;
double *NodalX;
double *y;
double *NodalY;
double *z;
double *Nodalz;
double *b; 
double *NodalB;
double *m1;
double *Nodalm1;
double *m2;
double *Nodalm2;
double *db;
double *dm1;
double *dm2;
double *dz;
double *dh;
double *S0;
double *nFriction;
double *NodalnFriction;

double mindh;
double max_lambda;

extern void* xcalloc(int items, int size);
extern void calculateNodalCoordinates(int N, int NumEl, double* X, double *Xnodes);
extern void InterpolateWithSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double*x,
	double* NodalX, double* func, double* NodalFunc, double* NodalDeriv);
extern void InterpolateWithGSLSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double*x,
	double* NodalX, double* func, double* NodalFunc, double* NodalDeriv);

void store_mesh(char *FullMesh)
{
	// read in the co-ordinates
	FILE *Coordinates = fopen(FullMesh,"r");
	
	// throw away the first line of the file
	char firstline[100];
	fgets(firstline, sizeof(firstline), Coordinates);

	char secondline[100];
	fgets(secondline, sizeof(secondline), Coordinates);
	sscanf(secondline, "%d", &NumEdges);	

	NumEl = NumEdges-1;
	NumNodes = NumEl*Np+2;
	
	// array to store the x, y and z values at the nodes
	x = malloc(NumEdges*sizeof(double));
	y = malloc(NumEdges*sizeof(double));
	z = calloc(NumEdges,sizeof(double));
	b = malloc(NumEdges*sizeof(double));	
	m1 = malloc(NumEdges*sizeof(double));	
	m2 = malloc(NumEdges*sizeof(double));	
	nFriction = malloc(NumEdges*sizeof(double));
	//double *PtS0 = malloc(NumEdges*sizeof(double));

	// variable to temporarily store the node numbers while reading the file
	int nodes;
	
	int j;
	for (j=0; j<NumEdges; ++j)
	{	
		char grid[10000];
		fgets(grid, sizeof(grid), Coordinates);
		sscanf(grid, "%d %lf %lf %lf %lf %lf %lf %lf", &nodes, &x[j], &y[j], &z[j], &b[j], &m1[j], &m2[j], &nFriction[j]);
	//sscanf(grid, "%d %lf %lf %lf %lf %lf %lf %lf", &nodes, &x[j], &y[j], &PtS0[j], &b[j], &m1[j], &m2[j], &nFriction[j]);

	}

	fclose(Coordinates);

	dm1 = xcalloc(NumNodes-2,sizeof(double));
	dm2 = xcalloc(NumNodes-2,sizeof(double));
	dh = malloc(NumEl*sizeof(double));
	NodalX = xcalloc(NumNodes-2, sizeof(double));
	NodalY = xcalloc(NumNodes-2, sizeof(double));
	NodalB = xcalloc(NumNodes-2, sizeof(double));
	Nodalm1 = xcalloc(NumNodes-2, sizeof(double));
	Nodalm2 = xcalloc(NumNodes-2, sizeof(double));
	Nodalz = xcalloc(NumNodes-2, sizeof(double));
	NodalnFriction = xcalloc(NumNodes-2,sizeof(double));
	dz = xcalloc(NumNodes-2, sizeof(double));
	db = xcalloc(NumNodes-2, sizeof(double));
 	
	calculateNodalCoordinates(P, NumEl, x, NodalX);
 	calculateNodalCoordinates(P, NumEl, y, NodalY);
 	InterpolateWithGSLSplines(NumEdges, NumNodes-2, x, NodalX, b, NodalB, db);
	InterpolateWithGSLSplines(NumEdges, NumNodes-2, x, NodalX, z, Nodalz, dz); 
	InterpolateWithGSLSplines(NumEdges, NumNodes-2, x, NodalX, m1, Nodalm1, dm1);
	InterpolateWithGSLSplines(NumEdges, NumNodes-2, x, NodalX, m2, Nodalm2, dm2);

	//calculateNodalCoordinates(P, NumEl, z, Nodalz);
 	//calculateNodalCoordinates(P, NumEl, b, NodalB);
	//calculateNodalCoordinates(P, NumEl, m1, Nodalm1);
 	//calculateNodalCoordinates(P, NumEl, m2, Nodalm2);
	//calculateNodalCoordinates(P, NumEl, PtS0, dz);
	//InterpolateWithGSLSplines(NumEdges, NumNodes-2, x, NodalX, PtS0, dz, NULL);
	//InterpolateWithSplines(NumEdges, NumNodes-2, x, NodalX, z, Nodalz, dz); 

	//calculateNodalCoordinates(P, NumEl, nFriction, NodalnFriction);
	InterpolateWithSplines(NumEdges, NumNodes-2, x, NodalX, nFriction, NodalnFriction, NULL); 
	
	mindh = 100000;

	for (int i =0; i < NumEl; ++i)
	{
		dh[i] = sqrt((x[i+1] - x[i])*(x[i+1]-x[i])+(y[i+1] - y[i])*(y[i+1]-y[i]));
		mindh = fmin(mindh, dh[i]);

		//double S0val = (z[i+1]-z[i])/dh[i];
	 	//double dbval = (b[i+1] - b[i])/dh[i];
		for (int j = 0; j < Np; j++)
		{
			//dz[i*Np+j] = S0val;
			//db[i*Np+j] = dbval;
		}
	
	}

}

