#include <stdio.h>

// a function that returns -1 if a <0, 0 if a = 0 and 1 if a >1
int sign(double a)
{
	if (a < 0)
		return (-1);
	if (a == 0)
		return (0);
	else
		return (1);
}

// A function to compute s = a*v1 + b*v2, where a and b are scalars and v1 and v2 are vectors

void addVectors (double vector1[], double vector2[], double sum[], int size, double a, double b)
{		
	int i;
	for (i =0; i < size; ++i)
			sum[i] = a*vector1[i] + b*vector2[i];

}

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
