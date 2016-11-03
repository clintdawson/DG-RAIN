#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"
#include "globals.h"
#include "oneTimeStep.h"
#include "constitutive_equations.h"

double time;
void oneTimeStep (double time, double dt)
{

	// Calculate the wet-dry status of the elements
	#ifdef WDON
	wetDryStatus();
	#endif

	// Intermediate RK step
	computeL();
	
	double oldA[NumNodes];
	double oldQ[NumNodes];
	for (int i =0; i < NumNodes; ++i)
	{
		oldA[i] = A[i];
		oldQ[i] = Q[i];
		A[i] = A[i] + dt*RHSA[i];
		Q[i] = Q[i] + dt*RHSQ[i];
		//printf("RHSA = %lf \n" ,RHSA[i]);
		//printf("RHSQ = %lf\n", RHSQ[i]);
	}


	#ifdef DEBUG
	printf("After first computeL\n");
	printf("A\n");	
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",A[i]);
	printf("Q\n");
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",Q[i]);
	#endif

	// Apply minmod slope limiter on w_1
	minmod();
	
	#ifdef DEBUG
	printf("After first minmod\n");
	printf("A\n");	
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",A[i]);
	printf("Q\n");
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",Q[i]);
	#endif

	#ifdef WDON
	// apply the positive-depth operator on w_1
	PDop();
	#endif

	#ifdef DEBUG
	printf("After second PDop \n");
	printf("A\n");	
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",A[i]);
	printf("Q\n");
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",Q[i]);
	#endif


	#ifdef WDON
	// Calculate the wet-dry status of the elements
	wetDryStatus();
	#ifdef DEBUG
	for (int i =0; i < NumEl; ++i)
	{
		printf("%d   ", WD[i]);
	}
	printf("\n");
	#endif
	#endif

	boundary_conditions();

	// Compute w
	computeL();

	for (int i =0; i < NumNodes; ++i)
	{
		A[i] = (oldA[i]+A[i]+dt*RHSA[i])/2;
		Q[i] = (oldQ[i]+Q[i]+dt*RHSQ[i])/2;
	}

	#ifdef DEBUG
	printf("After second computeL\n");
	printf("A\n");	
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",A[i]);
	printf("Q\n");
	for(int i =0; i < NumNodes; ++i)
		printf("%1.15f\n",Q[i]);

	#endif 

	// Apply minmod slope limiter to w
	minmod();


	#ifdef DEBUG
		printf("After second minmod\n");
		printf("A\n");	
		for(int i =0; i < NumNodes; ++i)
			printf("%1.15f\n",A[i]);
		printf("Q\n");
		for(int i =0; i < NumNodes; ++i)
			printf("%1.15f\n",Q[i]);
		#endif

		#ifdef WDON
		// Apply the positive depth operator on w
		PDop();

		#ifdef DEBUG
		printf("After third PDop\n");
		printf("A\n");	
		for(int i =0; i < NumNodes; ++i)
			printf("%1.15f\n",A[i]);
		printf("Q\n");
		for(int i =0; i < NumNodes; ++i)
			printf("%1.15f\n",Q[i]);
		#endif
		#endif

		boundary_conditions();
	}

	void time_evolution(double FinalTime)
	{
			
/*		//Output initial conditions 
		FILE* file1;
		FILE* file2;
		
		char filename1[100];
		char filename2[100];
		sprintf(filename1, "000_Zeta.dat");
		sprintf(filename2, "000_Q.dat");
		file1 = fopen(filename1, "w");
		file2 = fopen(filename2, "w");
			
		fprintf(file1, "Zeta at time 0\n"); 
		fprintf(file2, "Q at time 0 \n");
			
		for (int k = 0; k < NumEl; k++)
		{
			for (int j = 0; j < Np; j++)
			{	
				int nodeNum = k*Np+j;
				double bval = NodalB[nodeNum];
				double xval = NodalX[nodeNum];
				double yval = NodalY[nodeNum];
				double zval = Nodalz[nodeNum];
				double m1val = Nodalm1[nodeNum];
				double m2val = Nodalm2[nodeNum];
				double Aval = A[nodeNum+1];
				double Qval = Q[nodeNum+1];
				double Hval = getH(Aval, bval, m1val, m2val);
				double zetaVal = Hval - zval;
				
				fprintf(file1, "%e\t%e\t%e\t%e\t%e\n", xval, yval, -zval, Hval, zetaVal);
				fprintf(file2, "%e\t%e\t%e\t%e\n", xval, yval, -zval, Qval);

			}

		}

		fclose(file1);
		fclose(file2);
	 
*/	 
		time = 0;
		dt = 1e-5;
		int Nstep = 0;
		int fileNumber = 1;


		while (time<FinalTime)
		{
			time = time + dt;
			printf("dt = %e\n", dt);
			printf("time = %e\n", time);	
			
			oneTimeStep(time,dt);

			boundary_conditions();

			/********* Output data every twentienth step *************/		
			Nstep = Nstep + 1;
			int NstepinOneSec = floor(1./dt);
			if ( time == FinalTime)
			//if (Nstep%1000 == 0 || time == FinalTime)
			//if (Nstep%(NstepinOneSec/10) == 0 || time == FinalTime)
			{
			
				//Files for the channels 
				FILE* file1;
				FILE* file2;
				
				char filename1[100];
				char filename2[100];
				sprintf(filename1, "%3.3d_Zeta.dat", fileNumber);
				sprintf(filename2, "%3.3d_Q.dat", fileNumber);
				file1 = fopen(filename1, "w");
				file2 = fopen(filename2, "w");
			
				fprintf(file1, "Zeta at time %e for Channel \n", time);
				fprintf(file1, "x\t\t y \t\t -z \t\t height \t zeta\n"); 
				fprintf(file2, "Q at time %e for Channel \n", time);
		
				for (int k = 0; k < NumEl; k++)
				{
					for (int j = 0; j < Np; j++)
					{	
						int nodeNum = k*Np+j;
						double bval = NodalB[nodeNum];
						double xval = NodalX[nodeNum];
						double yval = NodalY[nodeNum];
						double zval = Nodalz[nodeNum];
						double Aval = A[nodeNum+1];
						double Qval = Q[nodeNum+1];
						double m1val = Nodalm1[nodeNum]; 
						double m2val = Nodalm2[nodeNum];
						double Hval = getH(Aval, bval, m1val, m2val);
						//double Hval = Aval/bval;
						double zetaVal = Hval - zval;
						
						fprintf(file1, "%e\t%e\t%e\t%e\t%3.14f\n", xval, yval, -zval, Hval, zetaVal);
						fprintf(file2, "%e\t%e\t%e\t%3.14f\n", xval, yval, -zval, Qval);

					}

				}
			
	/*	
				for (int j = 0; j<NumEl; ++j)
				{
					double Bval  = (b[j]+b[j+1])/2;
					double x_val = (x[j]+x[j+1])/2;
					double y_val = (y[j]+y[j+1])/2;
					double zval = (z[j]+z[j+1])/2;	

					double Aval = (A[j*2+1] + A[j*2+2])/2;
					double Qval = (Q[j*2+1] + Q[j*2+2])/2;
					double H = Aval/Bval;
					double zeta = H - zval;

					fprintf(file1, "%e\t%e\t%e\t%e\t%e\n", x_val, y_val, -zval,H, zeta);
					fprintf(file2, "%e\t%e\t%e\t%e\n", x_val, y_val, -zval, Qval);
						
				}

	*/
			fclose(file1);
			fclose(file2);
			
				++fileNumber;
			}
				
			
			/**************** Set next time step ****************/

							
			// set next time step
			//printf("time evolution max_lambda = %e\n", max_lambda);
			//printf("mindh = %e\n", mindh);
			dt = 0.9*mindh/max_lambda;
			if ((time + dt) > FinalTime)
				dt = FinalTime - time;
			//printf("dt = %lf\n",dt);

	}
	
		
}
