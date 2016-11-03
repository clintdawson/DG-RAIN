#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MeshAttributes.h"
#include "SimulationSteps.h"

extern void extractPointstoPlot(int NumElx);

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		printf("Usage: ./Simulation 'Fort.14'\n");
		printf("Fort.14: Mesh that includes the entire domain\n");
		exit(EXIT_FAILURE);

	}
	
	char *FullMesh = argv[1];

	//create EltoVert, EdgtoEls and EdgtoVert
	store_mesh(FullMesh);
	extract_mesh_attributes();
	
	int NumElx = 96;
	//extractPointstoPlot(NumElx);
	

	double FinalTime;

	printf("Enter the time you would like to run the simulation till:\n");
	scanf("%lf", &FinalTime);
	

	// initialize
	initialize();	

	// Impose boundary conditions and couple the channels to the junction 
	boundary_conditions();

	// Step through time 

	time_evolution(FinalTime);

	return(0);


}


 
	
	
