#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

void *xcalloc(int items, int size)
{
	void *ptr = calloc(items, size);
	if (ptr == NULL)
	{
		printf("Unable to allocate memory\n");
		exit(EXIT_FAILURE);
	}
	return ptr;
}

void GetLGLPoints(int N, double *x)
{
	if (N == 0)
	{
		printf("Support for constants doesn't exist\n");
		exit(1);
	}
	else if (N == 1)
	{
		x[0] = -1.0;
		x[1] = 1.0;
	}
	else if (N==2)
	{
		x[0] = -1.0;
		x[1] = 0.0;
		x[2] = 1.0;
	}
	else if (N==3)
	{
		x[0] = -1.0;
		x[1] = -0.4472;
		x[2] = 0.4472;
		x[3] = 1.0;
	}
	else if (N==4)
	{
		x[0] = -1.0;
		x[1] = -0.6547;
		x[2] = 0;
		x[3] = 0.6547;
		x[4] = 1.0;
	}
	else
	{
		printf("Currently only have up to fourth order polynomial encoded for LGL points\n");
		exit(1);
	}
}

void GetLGLWeights(int N, double *w)
{
	if (N == 0)
	{
		printf("Support for constants doesn't exist\n");
		exit(1);
	}
	else if (N == 1)
	{
		w[0] = 1.0;
		w[1] = 1.0;
	}
	else if (N==2)
	{
		w[0] = 0.3333333333333;
		w[1] = 1.3333333333333;
		w[2] = 0.3333333333333;
	}
	else if (N==3)
	{
		w[0] = 0.16666666666667;
		w[1] = 0.83333333333333;
		w[2] = 0.83333333333333;
		w[3] = 0.16666666666667;
	}
	else if (N==4)
	{
		w[0] = 0.1;
		w[1] = 0.5444444;
		w[2] = 0.7111111;
		w[3] = 0.5444444;
		w[4] = 0.1;
	}
	else
	{
		printf("Currently only have up to fourth order polynomial encoded for LGL points\n");
		exit(1);
	}

}

void calculateNodalCoordinates(int N, int NumEl, double* X, double* Xnodes)
{
	double r[N+1];
	GetLGLPoints(N, r);
	for(int i = 0; i < NumEl; i++)
	{
		for (int j = 0; j < N+1; j++)
		{
			Xnodes[i*(N+1)+j] = X[i] + 0.5*(r[j]+1)*(X[i+1]-X[i]);
		}
	}

}

void LegendrePoly(double *x, int Np, int N, double* P)
{
	// Purpose: Evaluates the normalized Legendre polynomials 
	// 					at points x for order N and returns the values in the
	// 					array P that is Np long
	
	double aold = 0.0, anew=0.0;

	double **PL = xcalloc(N+1, sizeof(double*));
	for (int i =0; i < N+1; i++)
	{
		PL[i] =  xcalloc(Np, sizeof(double));
	}

	// Initial values P_0(x) and P_1(x)
	if (N==0)
	{
		double constVal = 1.0/sqrt(2);
		for (int j = 0; j < Np; j++)
		{
			P[j] = constVal;
		}
		
		free(PL[0]);
		free(PL); 
		return;
	}
	else
	{
		for (int j = 0; j < Np; j++)
		{
			PL[0][j] = 1.0/sqrt(2);
		}
	}
	
	double *prow = xcalloc(Np, sizeof(double));
	double coeff = sqrt(1.5);
	for (int i = 0; i < Np; i++)
	{
		prow[i] = coeff*x[i];
	}

	if(N==1)
	{
		for (int i = 0; i < Np; i++)
		{
			P[i] = prow[i];
		}
		free(prow);
		free(PL[0]);
		free(PL[1]);
		free(PL);
		return;
	}
	else
	{
		for (int j =0; j < Np; j++)
		{
			PL[1][j] = coeff*x[j];
		}

	}

	// Repeat value in recurrence
	aold = sqrt(1.0/3);

	// Forward recurrence
	for (int i = 2; i <= N; i++)
	{
		anew = sqrt(i*i/((2.0*i+1)*(2*i-1)));
		for (int j = 0; j < Np; j++)
		{
			PL[i][j] = (x[j]*PL[i-1][j]-aold*PL[i-2][j])/anew;
		}
		aold = anew;
	}

	for (int j = 0; j < Np; j++)
	{
		P[j] = PL[N][j];
	}

	for (int i =0; i < N; i++)
	{
		free(PL[i]);
	}
	free(PL);

	return;
}

void Vandermonde1D(double *x, int Np, int N, gsl_matrix** V1D)
{
	// Purpose: Evaluates the Vandermonde matrix V_{ij} = phi_i(r_j);
	
	*V1D = gsl_matrix_alloc(Np, N+1);
	double V1DT[Np];
	for(int i = 0; i < N+1; i++)
	{
			LegendrePoly(x, Np, i, V1DT);
			for (int j = 0; j < Np; j++)
			{
				gsl_matrix_set(*V1D, i, j, V1DT[j]);
				//gsl_matrix_set(*V1D, j, i, V1DT[j]);
			}
	}
}

void GradLegendrePoly(double *x, int Np, int N, double* DP)
{
	// Purpose: Evaluates the derivatives of normalized Legendre polynomial  
	// 					at points x for order N and returns the values in the
	// 					array P that is Np long
	//
	
	double aold, anew;

	double **DPL = xcalloc(N+1, sizeof(double*));
	for (int i =0; i < N+1; i++)
	{
		DPL[i] =  xcalloc(Np, sizeof(double));
	}

	// Initial values P_0(x) and P_1(x)
	if (N==0)
	{
		for (int j = 0; j < Np; j++)
		{
			DP[0] = 0;
		}
		free(DPL[0]);
		free(DPL); 
		return;
	}
	else
	{
		for (int j = 0; j < Np; j++)
		{
			DPL[0][j] = 0;
		}
	}
	
	double *dprow = xcalloc(Np, sizeof(double));
	double deriv = sqrt(1.5);
	for (int i = 0; i < Np; i++)
	{
		dprow[i] = deriv;
	}

	if(N==1)
	{
		for (int i = 0; i < Np; i++)
		{
			DP[i] = deriv;
		}
		free(dprow);
		free(DPL[0]);
		free(DPL[1]);
		free(DPL);
		return;
	}
	else
	{
		for (int j =0; j < Np; j++)
		{
			DPL[1][j] = deriv;
		}

	}

	// Repeat value in recurrence
	aold = sqrt(1.0/3);

	// Forward recurrence
	for (int i = 2; i <= N; i++)
	{
		anew = sqrt(i*i/((2.0*i+1)*(2*i-1)));
		double P[Np];
		LegendrePoly(x, Np, i-1, P);
		for (int j = 0; j < Np; j++)
		{
			DPL[i][j] = 1.0/(anew)*(P[j] + x[j]*DPL[i-1][j] - aold*DPL[i-2][j]);
		}
		aold = anew;
	}

	for (int j = 0; j < Np; j++)
	{
		DP[j] = DPL[N][j];
	}

	for (int i =0; i < N; i++)
	{
		free(DPL[i]);
	}
	free(DPL);

	return;

}

void GradVandermonde1D(double *r, int Np, int N, gsl_matrix** GV)
{
	// Purpose: Evaluates the gradient of modal basis (i) at (r)
	// 					at order N and return them in GV whose size is
	// 					ixNp
	//
	*GV = gsl_matrix_alloc(Np, N+1);
	double GVT[Np];
	for (int i =0; i < N+1; i++)
	{
		GradLegendrePoly(r, Np, i, GVT);
		for (int j = 0; j < Np; j++)
			gsl_matrix_set(*GV,i,j, GVT[j]);
			//gsl_matrix_set(*GV,j,i,GVT[j]);
	}

}


void VolIntMat1D(double *x, int N, gsl_matrix* V, gsl_matrix *Vr, gsl_matrix** VolMat)
{
	
	gsl_permutation *p = gsl_permutation_alloc(N+1);
	gsl_matrix *Vcopy = gsl_matrix_alloc(N+1, N+1);
	gsl_matrix_memcpy(Vcopy, V);

	gsl_matrix *Vinv = gsl_matrix_alloc(N+1, N+1);
	
	int s;
	gsl_linalg_LU_decomp(Vcopy, p, &s);
	gsl_linalg_LU_invert(Vcopy, p, Vinv); 

	gsl_matrix *tmp1 = gsl_matrix_alloc(N+1, N+1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Vinv, Vinv, 0.0, tmp1);

	gsl_matrix *tmp2 = gsl_matrix_alloc(N+1, N+1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vr, tmp1, 0.0, tmp2);

	*VolMat = gsl_matrix_alloc(N+1, N+1);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, tmp2, 0.0, *
VolMat);
	//gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vr, Vinv, 0.0, *Dr);

	gsl_permutation_free(p);
	gsl_matrix_free(Vcopy);
	gsl_matrix_free(Vinv);
	gsl_matrix_free(tmp1);
	gsl_matrix_free(tmp2);

}

void Lift1D(int Np, gsl_matrix *V, gsl_matrix **LIFT)
{
	*LIFT = gsl_matrix_alloc(Np,2);
	
	gsl_matrix *EMAT = gsl_matrix_calloc(Np, 2);
	gsl_matrix_set(EMAT, 0, 0, 1.0);
	gsl_matrix_set(EMAT, Np-1, 1, -1.0);
	
	gsl_matrix *tmp = gsl_matrix_alloc(Np,2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, EMAT, 0.0, tmp);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, tmp, 0.0, *LIFT);

	gsl_matrix_free(EMAT);
	gsl_matrix_free(tmp);

}

void calculateMassMatrix(int Np, gsl_matrix* V, gsl_matrix** MassMatrix)
{
	*MassMatrix = gsl_matrix_alloc(Np,Np);
	gsl_matrix *MassMatrixInv = gsl_matrix_alloc(Np,Np);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, V, 0.0, MassMatrixInv);
	
	int s = 0;	
	gsl_permutation *p = gsl_permutation_alloc(Np);
	gsl_linalg_LU_decomp(MassMatrixInv, p, &s);
	gsl_linalg_LU_invert(MassMatrixInv, p, *MassMatrix);

	gsl_matrix_free(MassMatrixInv);
	gsl_permutation_free(p);

}

void calculateLIFTVolMat(int N, int Np, gsl_matrix **LIFT, gsl_matrix **VolMat,gsl_matrix **MassMatrix)
{
	double r[Np];
	GetLGLPoints(N,r);
	
	gsl_matrix *V; 
	gsl_matrix *GV;

	Vandermonde1D(r,Np,N,&V);
	GradVandermonde1D(r, Np, N, &GV);
	VolIntMat1D(r, N, V, GV, VolMat);
	Lift1D(Np, V, LIFT);
	calculateMassMatrix(Np, V,MassMatrix);

	gsl_matrix_free(V);
	gsl_matrix_free(GV);

}




