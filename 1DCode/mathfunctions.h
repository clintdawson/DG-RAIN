// Ensure that this header file is only included once

#ifndef ADDITIONAL_MATH_FUNCTIONS

	#define ADDITIONAL_MATH_FUNCTIONS

	// Define a macro for max, min, abs and sign functions and PI 
	#ifndef max
		#define max(a,b) (((a) > (b)) ? (a):(b))
	#endif

	#ifndef min
		#define min(a,b) (((a) < (b)) ? (a):(b))
	#endif

//	#ifndef abs
//		#define abs(a) (((a) >= 0) ? (a):(-a))
//	#endif
	
	#ifndef index
		#define index(a,b,col)  col*a+b 
	#endif

	#ifndef PI
		#define PI 3.14159265
	#endif

	extern int sign(double a);
	extern int mod (int a, int b);
	extern void addVectors (double vector1[], double vector2[], double sum[], int size,
		 double a, double b);
	extern void MatrixVectorMultiply(double *mat, double vec[], double prodvec[], int row, int col); 

#endif
