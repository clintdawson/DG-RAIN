#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <nubspline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

extern void* xcalloc(int NumPoints, int size);

/*void eval_NUBspline_1d_d(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val);

void eval_NUBspline_1d_d_vg(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val, double* restrict grad);

void eval_NUBspline_1d_d_vgl(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val, double* restrict grad, double* restrict lapl);


void InterpolateWithSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double* x, double* y, double* NodalX, double* NodalY, double *func, double *NodalFunc, double* NodalDeriv)
{

	double *s = xcalloc(NumPointsToInterpolate, sizeof(double));
	double *NodalS = xcalloc(NumPointsToEvaluate, sizeof(double));
	for (int j = 0; j < NumPointsToInterpolate; j++)
	{
		s[i] = sqrt(pow(x[j]-x[0],2)+pow(y[j]-y[0],2));

	}

	for(int j = 0; j < NumPointsToEvaluate, j++)
	{
		NodalS[i] = sqrt(pow(NodalX[j] -x[0],2) + pow(NodalY[j]-y[0],2));
	}

	NUgrid* myGrid = create_general_grid(s, NumPointsToInterpolate);
	BCtype_d myBC;
	myBC.lCode = NATURAL;
	myBC.rCode = NATURAL;
	NUBspline_1d_d *mySpline = create_NUBspline_1d_d(myGrid, myBC, func);
	for(int i =0; i < NumPointsToEvaluate; i++)
	{
		if (NodalDeriv != NULL)
			eval_NUBspline_1d_d_vg(mySpline, NodalS[i], &NodalFunc[i], &NodalDeriv[i]);
		else
			eval_NUBspline_1d_d(mySpline, NodalS[i], &NodalFunc[i]);
	}

	free(s);
	free(NodalS);

}
*/

void InterpolateWithGSLSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double* x, double* y, double* NodalX, double* NodalY, double* func, double* NodalFunc, double* NodalDeriv)
{
	double *s = xcalloc(NumPointsToInterpolate, sizeof(double));
	double *NodalS = xcalloc(NumPointsToEvaluate, sizeof(double));
	for (int j = 0; j < NumPointsToInterpolate; j++)
	{
		s[j] = sqrt(pow(x[j]-x[0],2)+pow(y[j]-y[0],2));

	}

	for(int j = 0; j < NumPointsToEvaluate; j++)
	{
		NodalS[j] = sqrt(pow(NodalX[j] -x[0],2) + pow(NodalY[j]-y[0],2));
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, NumPointsToInterpolate);
	gsl_spline_init(spline, s, func, NumPointsToInterpolate);

	for(int j = 0; j < NumPointsToEvaluate; j++)
	{
		NodalFunc[j] = gsl_spline_eval(spline, NodalS[j], acc);
		if (NodalDeriv != NULL)
			NodalDeriv[j] = gsl_spline_eval_deriv(spline, NodalS[j], acc);
	}

	gsl_interp_accel_free(acc);
	gsl_spline_free(spline);
	free(s);
	free(NodalS);

}
