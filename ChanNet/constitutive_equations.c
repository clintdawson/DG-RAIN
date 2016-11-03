/*********************************************************************************************//**
* @file constitutive_equations.c
* 
* This file contains definitions for calculating the pressure terms (I1 and I2) along with the
* friction slope (\f$S_f\f$ using Manning's formula) that appear in the 1-D St. Venant equations.
* Currently, calculation of I1 and I2 assumes that the channels have rectangular cross-sections.
* The code can be easily extended for trapezoidal cross-sections.
*
**************************************************************************************************/

#include <math.h>
#include "mathfunctions.h"

/**********************************************************************************************//**
* Function that calculates and returns a double value (I1) that represents the hydrostatic 
* pressure term at any point in a channel with rectangular cross-sections
* @param [in] A cross-sectional area at the quadrature point
* @param [in] b width of the channel at the quadrature point
*
* *************************************************************************************************/ 
double getI1(double A, double b)
{
	double I1;
	double h = A/b;
	I1 = pow(h,2)*b/2;
	return(I1);
}


/********************************************************************************************//**
* Function that calculates and returns a double value (I2) that represents the wall pressure
* term at any point in a channel with rectangular cross-sections
* @param [in] A cross-sectional area at the quadrature point
* @param [in] b width at the quadrature point
* @param [in] db value of the derivative of width in the element we are currently working on
*
************************************************************************************************/
double getI2(double A, double b, double db)
{
	double I2;
	double h = A/b;
	I2 = db*h*h/2;
	return (I2);
}


/*****************************************************************************************//**
* Function that calculates and returns a double value (\f$S_f\f$) that represents the friction
* slope term at any point in a channel. The friction slope is calculated using Manning's 
* formula.
* @param [in] A cross-sectional area at the quadrature point
* @param[in] Q volumetric discharge at the quadrature point
* @param[in] b width at the quadrature point
* @param[in] n value of Manning's coefficient at the quadrature point
*
* ********************************************************************************************/
double getS_f(double A, double Q, double b, double n)
{
	double h = A/b;
	double S_f = (pow(n,2)*Q*fabs(Q))/(pow(A,2)*pow(h,2./3));
	return (S_f);
} 
