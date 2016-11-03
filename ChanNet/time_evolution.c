/************************************************************************//**
* @file time_evolution.c
*
* This file contains code to advance the solutions one step forward in time
*
* *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ChannelsAndJunctions.h"
#include "MeshAttributes.h"
#include "mathfunctions.h"
#include "oneTimeStep.h"

/* @cond GLOBAL_VARIABLES */
extern const double g;
/* @endcond */

/***************************************************************************//**
*
* Function to advance a channel structure one step forward in time using 2-stage
* RKDG scheme
* @param[in] Chan channel structure that we are currently working on
* @param[in] time current time of the simulation
* @param[in] dt size of the time step
* @param[in] channelNumber identity (number) of the channel structure that we
* are currently working on
*
* *********************************************************************************/
void oneTimeStep (struct channel *Chan, double time, double dt, int channelNumber)
{
	int TotalNumNodes = 2*Chan->NumNodes;
	double RHSA[TotalNumNodes];
	double RHSQ[TotalNumNodes];

	// Calculate wet/dry status of elements
	#ifdef WDON
	wetDryStatus1D(Chan);
	#endif

	// store A and Q 
	gsl_vector *A = gsl_vector_alloc(TotalNumNodes);
	gsl_vector *Q = gsl_vector_alloc(TotalNumNodes);
	for (int i=0; i < TotalNumNodes; ++i) 
	{
		gsl_vector_set(A,i,(*Chan).A[i]);
		gsl_vector_set(Q,i,(*Chan).Q[i]);
	}

	/******** For debugging print out RHSA and RHSQ)************/
	#ifdef DEBUG
	printf("Intial A and Q:\n");
	gsl_vector_fprintf(stdout, A, "%1.16f");
	printf("\n");
	gsl_vector_fprintf(stdout, Q, "%1.16f");
	printf("\n");
	#endif
	/***********************************************************/
	
	// Intermediate RK step
	computeL(Chan, time, channelNumber, RHSA, RHSQ);
	
	gsl_vector *gRHSA = gsl_vector_alloc(TotalNumNodes);
	gsl_vector *gRHSQ = gsl_vector_alloc(TotalNumNodes);
	for (int i=0; i<TotalNumNodes; ++i)
	{
		gsl_vector_set(gRHSA,i,RHSA[i]);
		gsl_vector_set(gRHSQ,i,RHSQ[i]);
	}
	
	/******** For debugging print out RHSA and RHSQ)************/
	#ifdef DEBUG
	printf("RHSA and RHSQ after first computeL:\n");
	gsl_vector_fprintf(stdout, RHSA, "%1.16f");
	printf("\n");
	gsl_vector_fprintf(stdout, RHSQ, "%1.16f");
	printf("\n");
	#endif
	/***********************************************************/


	// w_1 = w + dt*L 
	gsl_vector *A_1 = gsl_vector_alloc(TotalNumNodes);
	gsl_vector *Q_1 = gsl_vector_alloc(TotalNumNodes);
	for (int i=0; i<TotalNumNodes; ++i)
	{
		gsl_vector_set(A_1,i,gsl_vector_get(A,i));
		gsl_vector_set(Q_1,i,gsl_vector_get(Q,i));
	}
	gsl_vector_scale(gRHSA,dt);
	gsl_vector_scale(gRHSQ,dt);
	gsl_vector_add(A_1,gRHSA);
	gsl_vector_add(Q_1,gRHSQ);

	/******** For debugging print out A_1 and Q_1************/
	#ifdef DEBUG
	printf("A_1 and Q_1\n");
	gsl_vector_fprintf(stdout, A_1, "%1.16f");
	printf("\n");
	gsl_vector_fprintf(stdout, Q_1, "%1.16f");
	printf("\n");
	#endif
	/***********************************************************/

	for (int i=0; i <TotalNumNodes; ++i)
	{
		Chan->A[i] = gsl_vector_get(A_1,i);
		Chan->Q[i] = gsl_vector_get(Q_1,i);
	}

	// Apply minmod slope limiter on w_1
	minmod(Chan);

	#ifdef WDON
	// apply the positive-depth operaton on w_1
	PDop1D(Chan);

	// Caculate wet/dry status again
	wetDryStatus1D(Chan);
	
	#endif 
	
//	internal_BC();
	boundary_conditions();

	for (int i = 0; i<TotalNumNodes; ++i)
	{
		gsl_vector_set(A_1,i,Chan->A[i]);
		gsl_vector_set(Q_1,i,Chan->Q[i]);
	}
	
	/******** For debugging print out A and Q************/
/*	#ifdef DEBUG
	printf("A and Q after first minmod\n");
	gsl_vector_fprintf(stdout, A_1, "%e");
	printf("\n");
	gsl_vector_fprintf(stdout, Q_1, "%e");
	printf("\n");
	#endif
*/	/***********************************************************/

	// Compute w
	computeL(Chan, time, channelNumber, RHSA, RHSQ);

	for (int i=0; i<TotalNumNodes; ++i)
	{
		gsl_vector_set(gRHSA, i, RHSA[i]);
		gsl_vector_set(gRHSQ, i, RHSQ[i]);
	}

	/******** For debugging print out RHSA and RHSQ)************/
	#ifdef DEBUG
	printf("RHSA and RHSQ after second computeL:\n");
	gsl_vector_fprintf(stdout, RHSA, "%1.16f");
	printf("\n");
	gsl_vector_fprintf(stdout, RHSQ, "%1.16f");
	printf("\n");
	#endif
	/***********************************************************/

	// w_1 = w + w_1 
	gsl_vector_add(A_1,A);
	gsl_vector_add(Q_1,Q);

	// w = (w1 + dt*L)/2
	gsl_vector_scale(A_1,0.5);
	gsl_vector_scale(Q_1,0.5);
	gsl_vector_scale(gRHSA,dt/2);
	gsl_vector_scale(gRHSQ,dt/2);
	
	gsl_vector_add(A_1,gRHSA);
	gsl_vector_add(Q_1,gRHSQ);
	
	for (int i=0; i < TotalNumNodes; ++i) 
	{
		Chan->A[i] = gsl_vector_get(A_1,i);
		Chan->Q[i] = gsl_vector_get(Q_1,i);
	}

	/******** For debugging print out A and Q)************/
	#ifdef DEBUG
	printf("A and Q after second computeL:\n");
	gsl_vector_fprintf(stdout, A_1, "%1.16f");
	printf("\n");
	gsl_vector_fprintf(stdout, Q_1, "%1.16f");
	printf("\n");
	#endif
	/***********************************************************/

	// Apply minmod slope limiter to w
	minmod(Chan);
	
	#ifdef WDON
	PDop1D(Chan);
	#endif

	gsl_vector_free(A);
	gsl_vector_free(Q);
	gsl_vector_free(A_1);
	gsl_vector_free(Q_1);
	gsl_vector_free(gRHSA);
	gsl_vector_free(gRHSQ);

	boundary_conditions();
}

/***************************************************************************//**
*
* Function to advance a junction structure one step forward in time using 2-stage
* RKDG scheme
* @param[in] junc the junction tructure that we are currently working on
* @param[in] time current time of the simulation
* @param[in] dt size of the time step
*
* *********************************************************************************/

void oneTimeStep2D (struct junction *junc, double time, double dt)
{

	int NumEl = junc->NumEl;
	int NumEdges = 3*NumEl;
	double RHSZeta[NumEdges], RHSQx[NumEdges], RHSQy[NumEdges];

	#ifdef WDON
	wetDryStatus2D(junc);
	#endif

	// store H, Qx, Qy 
	gsl_vector *zeta = gsl_vector_alloc(NumEdges);
	gsl_vector *Qx = gsl_vector_alloc(NumEdges);
	gsl_vector *Qy = gsl_vector_alloc(NumEdges);

	
	int j = 0;
	for (int i=0; i < NumEl; ++i) {
		for (int k=0; k < 3; ++k){
			gsl_vector_set(zeta,j,junc->zeta[index(i,k,3)]);
			gsl_vector_set(Qx,j, junc->Qx[index(i,k,3)]);
			gsl_vector_set(Qy,j, junc->Qy[index(i,k,3)]);
			++j;
		}
	}

	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("Initial values:\n zeta \n");
	gsl_vector_fprintf(stdout,zeta,"%e");
	printf("\n Qx \n");
	gsl_vector_fprintf(stdout,Qx,"%e");
	printf("\n Qy \n");
	gsl_vector_fprintf(stdout,Qy,"%e");
	printf("\n");
	#endif
	/***************************************************************************/

	// Intermediate RK step
	compute2DL(junc, time, RHSZeta, RHSQx, RHSQy);

	gsl_vector *gRHSZeta = gsl_vector_alloc(NumEdges);
	gsl_vector *gRHSQx = gsl_vector_alloc(NumEdges);
	gsl_vector *gRHSQy = gsl_vector_alloc(NumEdges);	

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
	gsl_vector_fprintf(stdout,RHSZeta,"%e");
	printf("\n RHSQx \n");
	gsl_vector_fprintf(stdout,RHSQx,"%e");
	printf("\n RHSQy \n");
	gsl_vector_fprintf(stdout,RHSQy,"%e");
	printf("\n");
	#endif
	/***************************************************************************/

	// w_1 = w + dt*L; 
	gsl_vector *zeta_1 = gsl_vector_alloc(NumEdges);
	gsl_vector *Qx_1 = gsl_vector_alloc(NumEdges);
	gsl_vector *Qy_1 = gsl_vector_alloc(NumEdges);
	for (int i=0; i<NumEdges; ++i)
	{
		gsl_vector_set(zeta_1,i,gsl_vector_get(zeta,i));
		gsl_vector_set(Qx_1,i,gsl_vector_get(Qx,i));
		gsl_vector_set(Qy_1,i,gsl_vector_get(Qy,i));
	}

	gsl_vector_scale(gRHSZeta,dt);
	gsl_vector_scale(gRHSQx,dt);
	gsl_vector_scale(gRHSQy,dt);
	gsl_vector_add(zeta_1,gRHSZeta);
	gsl_vector_add(Qx_1, gRHSQx);
	gsl_vector_add(Qy_1, gRHSQy);
	
	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("After first computeL:\n Zeta \n");
	gsl_vector_fprintf(stdout,zeta_1,"%e");
	printf("\n Qx \n");
	gsl_vector_fprintf(stdout,Qx_1,"%e");
	printf("\n Qy \n");
	gsl_vector_fprintf(stdout,Qy_1,"%e");
	printf("\n");
	#endif
	/***************************************************************************/
	
	j = 0;
	for (int i=0; i<NumEl; ++i){
		for (int k =0; k < 3; ++k){
			junc->zeta[index(i,k,3)] = gsl_vector_get(zeta_1,j);
			junc->Qx[index(i,k,3)] = gsl_vector_get(Qx_1,j);
			junc->Qy[index(i,k,3)] = gsl_vector_get(Qy_1,j);
			++j;
		}
	}

	// Apply minmod slope limiter on w_1
	SlopeLimiter(junc);
//	internal_BC();

	#ifdef WDON
	PDop2D(junc);
	wetDryStatus2D(junc);
	#endif

	j = 0;			
	for (int i=0; i<NumEl; ++i)
	{
		for(int k=0; k <3; ++k)
		{
			gsl_vector_set(zeta_1,j,junc->zeta[index(i,k,3)]);
			gsl_vector_set(Qx_1,j,junc->Qx[index(i,k,3)]);
			gsl_vector_set(Qy_1,j,junc->Qy[index(i,k,3)]);
			++j;
		}
	}
	
	/************ print out zeta, Qx and  Qy for debugging ********************/
	#ifdef DEBUG
	printf("After first slope limiting:\n zeta \n");
	gsl_vector_fprintf(stdout,zeta_1,"%e");
	printf("\n Qx \n");
	gsl_vector_fprintf(stdout,Qx_1,"%e");
	printf("\n Qy \n");
	gsl_vector_fprintf(stdout,Qy_1,"%e");
	printf("\n");
	#endif
	/***************************************************************************/
	
	// w_1 = w_1 + w;
	gsl_vector_add(zeta_1,zeta);
	gsl_vector_add(Qx_1, Qx);
	gsl_vector_add(Qy_1,Qy);

	// Compute w
	compute2DL(junc, time, RHSZeta, RHSQx, RHSQy);

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
	gsl_vector_fprintf(stdout,RHSZeta,"%e");
	printf("\n RHSQx \n");
	gsl_vector_fprintf(stdout,RHSQx,"%e");
	printf("\n RHSQy \n");
	gsl_vector_fprintf(stdout,RHSQy,"%e");
	printf("\n");
	#endif
	/***************************************************************************/

	gsl_vector_scale(zeta_1,0.5);
	gsl_vector_scale(Qx_1,0.5);
	gsl_vector_scale(Qy_1,0.5);
	gsl_vector_scale(gRHSZeta,dt/2);
	gsl_vector_scale(gRHSQx,dt/2);
	gsl_vector_scale(gRHSQy,dt/2);

	gsl_vector_add(zeta_1,gRHSZeta);
	gsl_vector_add(Qx_1,gRHSQx);
	gsl_vector_add(Qy_1,gRHSQy);

	j = 0;
	for (int i=0; i<NumEl; ++i){
		for (int k =0; k < 3; ++k){
			junc->zeta[index(i,k,3)] = gsl_vector_get(zeta_1,j);
			junc->Qx[index(i,k,3)] = gsl_vector_get(Qx_1,j);
			junc->Qy[index(i,k,3)] = gsl_vector_get(Qy_1,j);
			++j;
		}
	}	
	
	// Apply minmod slope limiter to w
	SlopeLimiter(junc);

	#ifdef WDON
	PDop2D(junc);
	#endif

	gsl_vector_free(zeta);
	gsl_vector_free(Qx);
	gsl_vector_free(Qy);
	gsl_vector_free(zeta_1);
	gsl_vector_free(Qx_1);
	gsl_vector_free(Qy_1);
	gsl_vector_free(gRHSZeta);
	gsl_vector_free(gRHSQx);
	gsl_vector_free(gRHSQy);

}

/*******************************************************************//**
*
* Function to evolve the channels as well as the junctions through time
* until the final simulation time. It also outputs the data at a certain
* time
*
* *********************************************************************/
void time_evolution(double FinalTime)
{
		
	//Output initial conditions 
	FILE* file1;
	FILE* file2;
	
	for (int i=0; i<NumChannels; ++i)
	{
		char filename1[100];
		char filename2[100];
		sprintf(filename1, "000_Zeta%d.dat", i);
		sprintf(filename2, "000_Q%d.dat", i);
		file1 = fopen(filename1, "w");
		file2 = fopen(filename2, "w");
			
		fprintf(file1, "Height at time 0 for Channel %d\n", i); 
		fprintf(file2, "Q at time 0 for Channel %d\n", i);
			
		int NumNodes = ChannelList[i]->NumNodes;	
		
		for (int j = 0; j<NumNodes; ++j)
		{	
			double B;
			double z;
			
			B = ChannelList[i]->b[j];
			z = ChannelList[i]->z[j];
			double x_val = ChannelList[i]->x[j];
			double y_val = ChannelList[i]->y[j];
			for (int k=0; k<2; ++k)
			{
		
				double A = ChannelList[i]->A[index(j,k,2)]; 
				double Q = ChannelList[i]->Q[index(j,k,2)];
				double H = A/B;
				double zeta = H - z;
				fprintf(file1, "%e\t%e\t%e\t%e\n", x_val, y_val, -z, H);
				fprintf(file2, "%e\t%e\t%e\t%e\n", x_val, y_val, -z, Q);
			}
		}
		fclose(file1);
		fclose(file2);
	}

	// output the file with x y and z information for the elements
/*	FILE *juncMesh;
	for (int i=0; i <NumJunctions; ++i)
	{
		char filename1[100];
		sprintf(filename1, "Junction%dMesh.dat",i);
		juncMesh = fopen(filename1, "w");
		int NumEl = JunctionList[i]->NumEl;
		int NumNodes = JunctionList[i]->NumNodes;
		int NumEdges = JunctionList[i]->TotalNumEdges;
		
		fprintf(juncMesh, "NodeNum \t x \t\t y \t\t z\n");
		for (int j = 0; j < NumNodes; ++j)
		{
			double xcor = JunctionList[i]->x[j];
			double ycor = JunctionList[i]->y[j];
			double zcor = JunctionList[i]->z[j];
			fprintf(juncMesh, "%d\t%e\t%e\t%e\n",j+1, xcor, ycor, zcor);
		}

	
		fprintf(juncMesh, "\n EltoVert \n");
		for (int j = 0; j < NumEl; ++j)
		{
			int n1 = JunctionList[i]->EltoVert[j*3]+1;
			int n2 = JunctionList[i]->EltoVert[j*3+1]+1;
			int n3 = JunctionList[i]->EltoVert[j*3+2]+1;
			fprintf(juncMesh, "%d\t3\t%d\t%d\t%d\n", j+1, n1, n2,n3);

		} 

		fprintf(juncMesh, "\n EdgtoVert \n");
		for (int j =0; j < NumEdges; ++j)
		{
			int n1 = JunctionList[i]->EdgtoVert[j*2]+1;
			int n2 = JunctionList[i]->EdgtoVert[j*2+1]+1;
			fprintf(juncMesh, "%d\t%d\n",n1,n2);

		}

		fprintf(juncMesh, "\n EdgtoEls \n");
		for (int j =0; j < NumEdges; ++j)
		{
			int n1 = JunctionList[i]->EdgtoEls[j*2];
			int n2 = JunctionList[i]->EdgtoEls[j*2+1];
			fprintf(juncMesh, "%d\t%d\n",n1,n2);

		}

		fclose(juncMesh);
		
	}
*/
	boundary_conditions();
	
	double time = 0;
	double dt = 1e-5;
	int Nstep = 0;
	int fileNumber = 1;
	
	double minLength = 99999999;	
	for (int i = 0; i < NumChannels; ++i)
	{
		minLength = fmin(minLength, ChannelList[i]->mindh);
	}

	for (int i = 0; i < NumJunctions; ++i)
	{
		minLength = fmin(minLength, JunctionList[i]->minEdgLength);
	}

	while (time < FinalTime)
	{
		time = time + dt;
		printf("dt = %e\n", dt);
		printf("time = %e\n", time);	
	
		double maxLambda = 0;

		// step forward in time
		for(int i=0; i<NumChannels; ++i)
		{
			oneTimeStep(ChannelList[i],time,dt,i);
			maxLambda = fmax(maxLambda, ChannelList[i]->max_lambda);
		
		}


		internal_BC();
	
		for(int i=0; i<NumJunctions; ++i)
		{
			oneTimeStep2D(JunctionList[i],time,dt);
			maxLambda = fmax(maxLambda, JunctionList[i]->max_lambda);
		}
		internal_BC();

		boundary_conditions();

		/********* Output data every twentienth step *************/		
		Nstep = Nstep + 1;
		int NstepinOneSec = floor(1./dt);
		if (Nstep%5000== 0 || time == FinalTime)
		//if (Nstep% NstepinOneSec == 0)
		{
		
			//Files for the channels 
			FILE* file1;
			FILE* file2;
			
			for (int i=0; i<NumChannels; ++i)
			{
				char filename1[100];
				char filename2[100];
				sprintf(filename1, "%3.3d_Zeta%d.dat", fileNumber,i);
				sprintf(filename2, "%3.3d_Q%d.dat", fileNumber,i);
				file1 = fopen(filename1, "w");
				file2 = fopen(filename2, "w");
			
				fprintf(file1, "Height at time %e for Channel %d\n", time, i); 
				fprintf(file2, "Q at time %e for Channel %d\n", time,i);
				
				int NumNodes = ChannelList[i]->NumNodes;			

				for (int j = 0; j<NumNodes; ++j)
				{
					double B;
					double z;
			
					B = ChannelList[i]->b[j];
					z = ChannelList[i]->z[j];
					double x_val = ChannelList[i]->x[j];
					double y_val = ChannelList[i]->y[j];
	
					int k;
					for (k=0; k<2; ++k)
					{
			
						double A = ChannelList[i]->A[index(j,k,2)]; 
						double Q = ChannelList[i]->Q[index(j,k,2)];
						double H = A/B;
						double zeta = H - z;
						fprintf(file1, "%e\t%e\t%e\t%e\n", x_val, y_val, -z, H);
						fprintf(file2, "%e\t%e\t%e\t%e\n", x_val, y_val, -z, Q);
					}

				
				}
				fclose(file1);
				fclose(file2);
			}
	
			FILE* file3;
			FILE* file4;
			FILE* file5;

			for (int i=0; i<NumJunctions;++i)
			{
				char filename1[100];
				char filename2[100];
				char filename3[100];
				sprintf(filename1, "%3.3d_Height_junc%d.dat", fileNumber,i);
				sprintf(filename2, "%3.3d_Qx_junc%d.dat", fileNumber,i);
				sprintf(filename3, "%3.3d_Qy_junc%d.dat", fileNumber,i);
				file3 = fopen(filename1, "w");
				file4 = fopen(filename2, "w");
				file5 = fopen(filename3, "w");
			
				fprintf(file3, "Zeta at time %e for Junction %d\n", time, i); 
				fprintf(file4, "Qx at time %e for Junction %d\n", time, i); 
				fprintf(file5, "Qy at time %e for Junction %d\n", time,i);
			
				int j;
				int k;
				for (j=0; j<JunctionList[i]->NumEl; ++j)
				{
					for(k=0; k<3; ++k)	
					{
						double zeta = JunctionList[i]->zeta[index(j,k,3)];
						double Qx = JunctionList[i]->Qx[index(j,k,3)];
						double Qy = JunctionList[i]->Qy[index(j,k,3)];
						double z1 = JunctionList[i]->z[JunctionList[i]->EltoVert[j*3+k]];
						double z2 = JunctionList[i]->z[JunctionList[i]->EltoVert[j*3+(k+1)%3]];
						double z = 0.5*(z1+z2);
						double height = zeta + z;
						fprintf(file3, "%1.15f\t\t",height);
						fprintf(file4, "%1.15f\t\t",Qx);
						fprintf(file5, "%1.15f\t\t",Qy);
					}
					fprintf(file3, "\n");
					fprintf(file4, "\n");
					fprintf(file5, "\n");
				}
				fclose(file3);
				fclose(file4);
				fclose(file5);
			}
				
			++fileNumber;
		}
			
		
		/**************** Set next time step ****************/

		        
		// set next time step
		//printf("max_lambda = %e\n", maxlambda);
		//printf("mindx = %e\n", mindx);
		dt = 0.4*minLength/maxLambda;
		if ((time + dt) > FinalTime)
			dt = FinalTime - time;

	}
}
