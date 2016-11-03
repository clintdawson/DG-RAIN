#ifndef SIMULATION_STEPS

#define SIMULATION_STEPS

#include <gsl/gsl_matrix.h>

extern void calculateLIFTVolMat(int N, int Np, gsl_matrix **LIFT, gsl_matrix **VolMat, gsl_matrix** MassMatrix);

extern void store_mesh(char *Mesh);
extern void initialize();
extern void time_evolution(double FinalTime);
extern void boundary_conditions();
#endif
