/*********** Functions that are definitions and thus should not be changed ***************/

#include <math.h>
#include "mathfunctions.h"

double getH(double A, double b, double m1, double m2)
{
	/** A  = wet cross-sectional area
	 *  b = width 
	 *  m1 = reciprocal of the slope of one side of the trapezoid (delta x / delta y) 
	 *  m2 = reciprocal of the slope of the other side of the trapezoid
	 *  Therefore, if m1 = 0 and m2 = 0, cross-section is rectangular
	 *  */
	
	double H;
	if (m1 > 0 || m2 > 0)
	{
		H = (-b+sqrt(b*b+ 2*A*m1 + 2*A*m2))/(m1+m2);
	}
	else 
		H = A/b;
	
	return H;
}

double getI1(double A, double b, double m1, double m2)
{
	/** A  = wet cross-sectional area
	 *  b = width 
	 *  m1 = reciprocal of the slope of one side of the trapezoid (delta x / delta y) 
	 *  m2 = reciprocal of the slope of the other side of the trapezoid
	 *  Therefore, if m1 = 0 and m2 = 0, cross-section is rectangular
	 *  */

	double I1;
	double h = getH(A, b, m1, m2);
	I1 = h*h*b/2 + h*h*h*m1/6 + h*h*h*m2/6;
	return(I1);
}
double getI2(double A, double b, double db, double m1, double dm1, double m2, double dm2)
{
	double I2;
	double h = getH(A, b, m1, m2);
	I2 = db*h*h/2 + h*h*h*dm1/6 + h*h*h*dm2/6;
	return (I2);
}

double getS_f(double A, double Q, double b, double m1, double m2, double n)
{
	double h = getH(A, b, m1, m2);
	double B = b + h*m1 + h*m2;
	double R = (0.5*h*(b+B))/(b+2*sqrt(pow((0.5*(B-b)),2) + h*h));
	double S_f = (pow(n,2)*Q*fabs(Q))/(pow(A,2)*pow(R,1.3333));
	//double S_f = (pow(n,2)*Q*fabs(Q))/(pow(A,2)*pow(h,4./3));
	return (S_f);
} 
