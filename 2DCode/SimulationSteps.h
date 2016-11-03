#ifndef SIMULATION_STEPS

#define SIMULATION_STEPS

extern void store_mesh(char *Mesh);
extern void extract_mesh_attributes();
extern void initialize();
extern void time_evolution(double FinalTime);
extern void boundary_conditions();
#endif
