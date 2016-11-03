#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nubspline.h>
#include "SimulationSteps.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

void eval_NUBspline_1d_d(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val);

void eval_NUBspline_1d_d_vg(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val, double* restrict grad);

void eval_NUBspline_1d_d_vgl(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val, double* restrict grad, double* restrict lapl);


void InterpolateWithSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double* x, double* NodalX, double *func, double *NodalFunc, double* NodalDeriv)
{

	NUgrid* myGrid = create_general_grid(x, NumPointsToInterpolate);
	BCtype_d myBC;
	myBC.lCode = NATURAL;
	myBC.rCode = NATURAL;
	NUBspline_1d_d *mySpline = create_NUBspline_1d_d(myGrid, myBC, func);
	for(int i =0; i < NumPointsToEvaluate; i++)
	{
		if (NodalDeriv != NULL)
			eval_NUBspline_1d_d_vg(mySpline, NodalX[i], &NodalFunc[i], &NodalDeriv[i]);
		else
			eval_NUBspline_1d_d(mySpline, NodalX[i], &NodalFunc[i]);
	}

}

void InterpolateWithGSLSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double* x, double* NodalX, double* func, double* NodalFunc, double* NodalDeriv)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, NumPointsToInterpolate);
	gsl_spline_init(spline, x, func, NumPointsToInterpolate);

	for(int j = 0; j < NumPointsToEvaluate; j++)
	{
		NodalFunc[j] = gsl_spline_eval(spline, NodalX[j], acc);
		if (NodalDeriv != NULL)
			NodalDeriv[j] = gsl_spline_eval_deriv(spline, NodalX[j], acc);
	}

}
