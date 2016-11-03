
/**********************************************************************************//**
* @file mathfunctions.c
*
* This file contains some simple math functions that I wasn't sure is included in 
* the C math library.
*
* *********************************************************************************/

#include <stdio.h>

/****************************************************************//**
 * a function that returns -1 if a <0, 0 if a = 0 and 1 if a >1 
 *
 * @param [in] a a double value whose sign we want to determine
 * ******************************************************************/
int sign(double a)
{
	if (a < 0)
		return (-1);
	if (a == 0)
		return (0);
	else
		return (1);
}

/***********************************************************************//**
* a function to compute s = a*vector1 + b*vector2
*
* @param [in] vector1 an array containing doubles
* @param [in] vector2 an array containing doubles
* @param [in] size size of the two vectors being added
* @param [in] a a double number by which the first vector is scaled
* @param [in] b a double number by which the second vector is scaled
* @param [out] sum a double array where the result will be stored
*
* ***************************************************************************/
void addVectors (double vector1[], double vector2[], double sum[], int size, double a, double b)
{		
	int i;
	for (i =0; i < size; ++i)
			sum[i] = a*vector1[i] + b*vector2[i];

}

/********************************************************************************//**
* a function to compute perform a matrix-vector multiplication: a = M*b
* @param [in] mat a pointer to the matrix M stored as an array in row-major order
* @param [in] vec a pointer to the vector b 
* @param [in] row integer number of rows of matrix M
* @param [in] col integer number of columns of matrix M
* @param [out] prodvec an array where the result will be stored
* ********************************************************************************/
void MatrixVectorMultiply(double *mat, double vec[], double prodvec[], int row, int col){

int i,j;

for (i=0;i<row;++i){
	double sum = 0;
	for(j=0;j<col;++j){
		sum +=mat[col*i+j]*vec[j] ;
	}
	prodvec[i] = sum;
}
}
