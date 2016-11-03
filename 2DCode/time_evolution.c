#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "MeshAttributes.h"
#include "mathfunctions.h"
#include "oneTimeStep.h"

extern const double g;

void oneTimeStep2D (double time, double dt)
{

	// store H, Qx, Qy 
	gsl_vector *gzeta = gsl_vector_alloc(3*NumEl);
	gsl_vector *gQx = gsl_vector_alloc(3*NumEl);
	gsl_vector *gQy = gsl_vector_alloc(3*NumEl);

	int j = 0;
	for (int i=0; i < NumEl; ++i) {
		for (int k=0; k < 3; ++k){
			gsl_vector_set(gzeta,j,zeta[index(i,k,3)]);
			gsl_vector_set(gQx,j, Qx[index(i,k,3)]);
			gsl_vector_set(gQy,j, Qy[index(i,k,3)]);
			++j;
		}
	}

	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("Initial values:\n zeta \n");
	gsl_vector_fprintf(stdout,gzeta,"%e");
	printf("\n Qx \n");
	gsl_vector_fprintf(stdout,gQx,"%e");
	printf("\n Qy \n");
	gsl_vector_fprintf(stdout,gQy,"%e");
	printf("\n");
	#endif
	/***************************************************************************/

	// Calculate the wet-dry status of the elements
//	wetDryStatus2D();

	// Intermediate RK step
	compute2DL(time);

	gsl_vector *gRHSZeta = gsl_vector_alloc(3*NumEl);
	gsl_vector *gRHSQx = gsl_vector_alloc(3*NumEl);
	gsl_vector *gRHSQy = gsl_vector_alloc(3*NumEl);	

	j = 0;
	for (int i=0; i < NumEl; ++i) {
		for (int k=0; k < 3; ++k){
			gsl_vector_set(gRHSZeta,j,RHSZeta[index(i,k,3)]);
			gsl_vector_set(gRHSQx,j,RHSQx[index(i,k,3)]);
			gsl_vector_set(gRHSQy,j,RHSQy[index(i,k,3)]);
			++j;
		}
	}

	/************ print out RHSZeta, RHSQx and RHSQy for debugging ********************/
	#ifdef DEBUG
	printf("After first computeL:\n RHSZeta \n");
	gsl_vector_fprintf(stdout,gRHSZeta,"%e");
	printf("\n RHSQx \n");
	gsl_vector_fprintf(stdout,gRHSQx,"%e");
	printf("\n RHSQy \n");
	gsl_vector_fprintf(stdout,gRHSQy,"%e");
	printf("\n");
	#endif
	/***************************************************************************/

	// w_1 = w + dt*L; 
	gsl_vector *gzeta_1 = gsl_vector_alloc(3*NumEl);
	gsl_vector *gQx_1 = gsl_vector_alloc(3*NumEl);
	gsl_vector *gQy_1 = gsl_vector_alloc(3*NumEl);
	for (int i=0; i<3*NumEl; ++i)
	{
		gsl_vector_set(gzeta_1,i,gsl_vector_get(gzeta,i));
		gsl_vector_set(gQx_1,i,gsl_vector_get(gQx,i));
		gsl_vector_set(gQy_1,i,gsl_vector_get(gQy,i));
	}

	gsl_vector_scale(gRHSZeta,dt);
	gsl_vector_scale(gRHSQx,dt);
	gsl_vector_scale(gRHSQy,dt);
	gsl_vector_add(gzeta_1,gRHSZeta);
	gsl_vector_add(gQx_1, gRHSQx);
	gsl_vector_add(gQy_1, gRHSQy);
	
	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("After first computeL:\n Zeta \n");
	gsl_vector_fprintf(stdout,gzeta_1,"%e");
	printf("\n Qx \n");
	gsl_vector_fprintf(stdout,gQx_1,"%e");
	printf("\n Qy \n");
	gsl_vector_fprintf(stdout,gQy_1,"%e");
	printf("\n");
	#endif
	/***************************************************************************/
	
	j = 0;
	for (int i=0; i<NumEl; ++i){
		for (int k =0; k < 3; ++k){
			zeta[index(i,k,3)] = gsl_vector_get(gzeta_1,j);
			Qx[index(i,k,3)] = gsl_vector_get(gQx_1,j);
			Qy[index(i,k,3)] = gsl_vector_get(gQy_1,j);
			++j;
		}
	}

/*	printf("zeta\n");
	for (int i=0; i < NumEl; ++i) 
	{
		for (int k =0; k <3; ++k)
		{
			printf("%1.15f\t",junc->zeta[index(i,k,3)]);
		}
		printf("\n");

	}

	printf("Qx\n");
	for (int i=0; i < NumEl; ++i) 
	{
		for (int k =0; k <3; ++k)
		{
			printf("%1.15f\t",junc->Qx[index(i,k,3)]);
		}
		printf("\n");

	}

	printf("Qy\n");
	for (int i=0; i < NumEl; ++i) 
	{
		for (int k =0; k <3; ++k)
		{
			printf("%1.15f\t",junc->Qy[index(i,k,3)]);
		}
		printf("\n");

	}
*/
	// Apply minmod slope limiter on w_1
	SlopeLimiter();

	// Apply the positive-depth operator on w_1
//	PDop2D();

	j = 0;			
	for (int i=0; i<NumEl; ++i)
	{
		for(int k=0; k <3; ++k)
		{
			gsl_vector_set(gzeta_1,j,zeta[index(i,k,3)]);
			gsl_vector_set(gQx_1,j,Qx[index(i,k,3)]);
			gsl_vector_set(gQy_1,j,Qy[index(i,k,3)]);
			++j;
		}
	}
	
	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("After first slope limiting:\n zeta \n");
	gsl_vector_fprintf(stdout,gzeta_1,"%e");
	printf("\n Qx \n");
	gsl_vector_fprintf(stdout,gQx_1,"%e");
	printf("\n Qy \n");
	gsl_vector_fprintf(stdout,gQy_1,"%e");
	printf("\n");
	#endif
	/***************************************************************************/
	
	// w_1 = w_1 + w;
	gsl_vector_add(gzeta_1,gzeta);
	gsl_vector_add(gQx_1,gQx);
	gsl_vector_add(gQy_1,gQy);

	// Calcualte the wet-dry status of the elements again
//	wetDryStatus2D();

	// Compute w
	compute2DL(time);

	j = 0;
	for (int i=0; i < NumEl; ++i) {
		for (int k=0; k < 3; ++k){
			gsl_vector_set(gRHSZeta,j,RHSZeta[index(i,k,3)]);
			gsl_vector_set(gRHSQx,j,RHSQx[index(i,k,3)]);
			gsl_vector_set(gRHSQy,j,RHSQy[index(i,k,3)]);
			++j;
		}
	}

	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("After second computeL:\n RHSZeta \n");
	gsl_vector_fprintf(stdout,gRHSZeta,"%e");
	printf("\n RHSQx \n");
	gsl_vector_fprintf(stdout,gRHSQx,"%e");
	printf("\n RHSQy \n");
	gsl_vector_fprintf(stdout,gRHSQy,"%e");
	printf("\n");
	#endif
	/***************************************************************************/

	gsl_vector_scale(gzeta_1,0.5);
	gsl_vector_scale(gQx_1,0.5);
	gsl_vector_scale(gQy_1,0.5);
	gsl_vector_scale(gRHSZeta,dt/2);
	gsl_vector_scale(gRHSQx,dt/2);
	gsl_vector_scale(gRHSQy,dt/2);

	gsl_vector_add(gzeta_1,gRHSZeta);
	gsl_vector_add(gQx_1,gRHSQx);
	gsl_vector_add(gQy_1,gRHSQy);

	j = 0;
	for (int i=0; i<NumEl; ++i){
		for (int k =0; k < 3; ++k){
			zeta[index(i,k,3)] = gsl_vector_get(gzeta_1,j);
			Qx[index(i,k,3)] = gsl_vector_get(gQx_1,j);
			Qy[index(i,k,3)] = gsl_vector_get(gQy_1,j);
			++j;
		}
	}	
	
	// Apply minmod slope limiter to w
	SlopeLimiter();

	// Apply the positive depth operator on w
//	PDop2D();

	gsl_vector_free(gzeta);
	gsl_vector_free(gQx);
	gsl_vector_free(gQy);
	gsl_vector_free(gzeta_1);
	gsl_vector_free(gQx_1);
	gsl_vector_free(gQy_1);
	gsl_vector_free(gRHSZeta);
	gsl_vector_free(gRHSQx);
	gsl_vector_free(gRHSQy);


}

void time_evolution(double FinalTime)
{

	// output the file with x y and z information for the elements
	boundary_conditions();
	
	double time = 0;
	double dt = 1e-7;
	int Nstep = 0;
	int fileNumber = 1;

	// write out initial conditions
	FILE* Initfile1;
	FILE* Initfile2;
	FILE* Initfile3;
//	FILE* Initfile4;

	char flname1[100];
	char flname2[100];
	char flname3[100];
//	char flname4[100];
	sprintf(flname1, "000_zeta.dat");
	sprintf(flname2, "000_Qx.dat");
	sprintf(flname3, "000_Qy.dat");
//	sprintf(flname4, "000_Qdotn.dat");
	Initfile1 = fopen(flname1, "w");
	Initfile2 = fopen(flname2, "w");
	Initfile3 = fopen(flname3, "w");
//	Initfile4 = fopen(flname4, "w");
	
/*	fprintf(Initfile1, "H at time 0\n"); 
	fprintf(Initfile2, "Qx at time 0 \n"); 
	fprintf(Initfile3, "Qy at time 0 \n");
	fprintf(Initfile4, "Qdotn at time 0 \n");
*/
	int k;
	for (int j=0; j<NumEl; ++j)
	{
		for(k=0; k<3; ++k)	
		{
			double zeta_val = zeta[index(j,k,3)];
			double zval = (z[EltoVert[j*3+k]]+z[EltoVert[j*3+((k+1)%3)]])/2;
			double H = zeta_val + zval;
			double Qx_val = Qx[index(j,k,3)];
			double Qy_val = Qy[index(j,k,3)];
		//	double Qdotn = Qx_val*normvecx[index(j,k,3)] + Qy_val*normvecy[index(j,k,3)];
			fprintf(Initfile1, "%1.15f\t\t",zeta_val);
			fprintf(Initfile2, "%1.15f\t\t",Qx_val);
			fprintf(Initfile3, "%1.15f\t\t",Qy_val);
		//	fprintf(Initfile4, "%1.15f\t\t",Qdotn);
		}
		fprintf(Initfile1, "\n");
		fprintf(Initfile2, "\n");
		fprintf(Initfile3, "\n");
//		fprintf(Initfile4, "\n");
	}
	fclose(Initfile1);
	fclose(Initfile2);
	fclose(Initfile3);
//	fclose(Initfile4);


	//wetDryStatus2D();
	while (time<FinalTime)
	{
		time = time + dt;
		printf("dt = %e\n", dt);
		printf("time = %e\n", time);	
		
		oneTimeStep2D(time,dt);
		
		boundary_conditions();

		/********* Output data every twentienth step *************/		
		Nstep = Nstep + 1;
		int NstepinSec = floor(1./dt);
		if (time == FinalTime)
	//	if (Nstep%500 == 0 || time == FinalTime)
//		if ((Nstep%(10*NstepinSec) == 0) || time == FinalTime)
		{
		
			FILE* file3;
			FILE* file4;
			FILE* file5;
//			FILE* file6;

			char filename1[100];
			char filename2[100];
			char filename3[100];
//			char filename4[100];
			
/*			sprintf(filename1, "Zeta%d", NumEl);
			sprintf(filename2, "Qx%d", NumEl);
			sprintf(filename3, "Qy%d", NumEl);
*/

			sprintf(filename1, "%3.3d_H.dat", fileNumber);
			sprintf(filename2, "%3.3d_Qx.dat", fileNumber);
			sprintf(filename3, "%3.3d_Qy.dat", fileNumber);
	
		//	sprintf(filename1, "%3.3d_H.dat", fileNumber);
		//	sprintf(filename2, "%3.3d_Qx.dat", fileNumber);
		//	sprintf(filename3, "%3.3d_Qy.dat", fileNumber);
//			sprintf(filename4, "%3.3d_Q.n.dat", fileNumber);
			file3 = fopen(filename1, "w");
			file4 = fopen(filename2, "w");
			file5 = fopen(filename3, "w");
//			file6 = fopen(filename4, "w");
			
			//fprintf(file3, "H at time %e\n", time); 
			//fprintf(file4, "Qx at time %e \n", time); 
			//fprintf(file5, "Qy at time %e \n", time);
//			fprintf(file6, "Qdotn at time %e \n", time);
			
			int k;
			for (int j=0; j<NumEl; ++j)
			{
				for(k=0; k<3; ++k)	
				{
					double zeta_val = zeta[index(j,k,3)];
					double Qx_val = Qx[index(j,k,3)];
					double Qy_val = Qy[index(j,k,3)];
					double zval = (z[EltoVert[j*3+k]]+z[EltoVert[j*3+((k+1)%3)]])/2;
					double hval = zeta_val + zval;

//					double Qdotn = Qx_val*normvecx[index(j,k,3)] + Qy_val*normvecy[index(j,k,3)];
					fprintf(file3, "%1.15f\t\t",zeta_val);
					fprintf(file4, "%1.15f\t\t",Qx_val);
					fprintf(file5, "%1.15f\t\t",Qy_val);
//					fprintf(file6, "%1.15f\t\t",Qdotn);
				}
				fprintf(file3, "\n");
				fprintf(file4, "\n");
				fprintf(file5, "\n");
//				fprintf(file6, "\n");
			}
			fclose(file3);
			fclose(file4);
			fclose(file5);
//			fclose(file6);
				
			++fileNumber;
		}
			
		
		/**************** Set next time step ****************/

		        
		// set next time step
	//	printf("max_lambda = %e\n", max_lambda);
	//	printf("minEdgLength = %e\n", minEdgLength);
		dt = 0.1*minEdgLength/max_lambda;
		if ((time + dt) > FinalTime)
			dt = FinalTime - time;

	}
}
