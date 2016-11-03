#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "globals.h"
#include "SimulationSteps.h"


//#define WDON

int P;
int Np;
gsl_matrix *LIFT;
gsl_matrix *VolMat;
gsl_matrix *MassMatrix;

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		printf("Usage: ./Simulation 'Fort.14'\n");
		printf("Fort.14: Mesh for the channel\n");
		exit(EXIT_FAILURE);

	}
	
	char *Mesh = argv[1];
	
	//printf("Enter the polynomial approximation order:\n");
	//scanf("%d", &P);

	P = 1;
	Np = P+1;


	store_mesh(Mesh);

/*	for (int i = 0; i < NumEl; i++)
	{
		for (int j = 0; j < Np; j++)
			printf("x = %lf \t z = %lf \t S0 = %3.18lf\n", NodalX[i*Np+j] , Nodalz[i*Np+j], dz[i*Np+j]);
	}

	exit(1);
*/	// Create VolMat, LIFT and MassMatrix
	calculateLIFTVolMat(P, Np, &LIFT, &VolMat, &MassMatrix);

	double FinalTime;

	printf("Enter the time you would like to run the simulation till:\n");
	scanf("%lf", &FinalTime);
	

	initialize();

	/* Impose boundary conditions and couple the channels to the junction */
	boundary_conditions();

	/* Step through time */

	time_evolution(FinalTime);

	return(0);

}


 
	
	
