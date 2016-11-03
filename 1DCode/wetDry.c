/******************** This file contains all the functions associated with the wet/dry treatment *****************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "globals.h"

#ifdef WDON
int BigTheta(double a)
{
	if (a > 0)
		return(1);
	else
		return(0);
}

void wetDryStatus()
{
	for (int i =0; i < NumEl; ++i)
	{
		double A0 = A[2*i+1];
		double A1 = A[2*i+2];
		double b0 = b[i];
		double b1 = b[i+1];
		double h0 = A0/b0;
		double h1 = A1/b1;
		double avgH = (h0+h1)/2;

		int maxHeightIndex;
		double maxHeight;
		if (h0 > h1)
		{
			maxHeightIndex = i;
			maxHeight = h0;
		}
		else
		{
			maxHeightIndex = i+1;
			maxHeight = h1;
		}
		
		double ZetaMax = maxHeight - z[maxHeightIndex];
		double z0 = z[i];
		double z1 = z[i+1];
		double zmin = z0*(z0 < z1) + z1*(z1 < z0); 

		if (WD[i] == -1) 	// if uninitialized
			WD[i] = BigTheta(avgH-H0)*BigTheta(ZetaMax - (H0-zmin));
		else
			WD[i] = WD[i]*BigTheta(avgH - H0)+(1-WD[i])*BigTheta(avgH-H0)*BigTheta(ZetaMax - (H0 - zmin));

	}


}

// post processing to ensure that the water depth remains positive
void PDop()
{
/*	printf("Initial A and Q \n");
	printf("A\n");
	for (int i =0; i < NumEl; ++i)
	{
		for (int k =0; k <2; ++k)
		{
			printf("%e\n",A[i*2+k]);
		}
		
	}
	printf("\n");
	printf("Q\n");
	for (int i =0; i < NumEl; ++i)
	{
		for (int k =0; k <2; ++k)
		{
			printf("%e\n",Q[i*2+k]);
		}
	}
*/
	// redistribute mass as necessary
	for (int i =0; i < NumEl; ++i)
	{
		double A0 = A[2*i+1];
		double A1 = A[2*i+2];
		double h0 = A0/b[i];
		double h1 = A1/b[i+1];
		double avgH = (h0+h1)/2;
		
//		printf("avgH = %e \n", avgH);

		if (h0 >= H0 && h1 >= H0)
		{	
			h0 = h0;
			h1 = h1;
		}		
	
		else if (avgH <= H0)
		{
			h0 = avgH;
			h1 = avgH;
//			printf("avgH is less than H0\n");
		}

		else
		{
//			printf("h0 = %e \t h1 = %e\n", h0, h1);
//			printf("avgH is greater than H0\n");
			if (h0 <= h1)
			{	
				double oldh0 = h0;
				h0 = H0;
				h1 = h1 - (h0 - oldh0);
			}

			else
			{
				double oldh1 = h1;
				h1 = H0;
				h0 = h0 - (h1-oldh1);
			}

		}

		A0 = h0*b[i];
		A1 = h1*b[i+1];
		A[2*i+1] = A0;
		A[2*i+2] = A1;
		if (isnan(A0) || isnan(A1))
		{
			printf("water height negative during wetting and drying treatment for el %d with avgH = %e\n", i, avgH);
			exit(EXIT_FAILURE);
		}
	
	}			

	// redistribute momentum as necessary
	for (int el = 0; el < NumEl; ++el)
	{
		int wetDryFlagNode[2];
		double u[2];
		for (int i =0; i < 2; ++i)
		{
			double Anode = A[2*el+i+1];
			double bnode = b[el+i];
			double hnode = Anode/bnode;
			wetDryFlagNode[i] = BigTheta(hnode - H0);
			
			double Qnode = Q[2*el+i+1];
			u[i] = Qnode/Anode;
	
		}
		if (wetDryFlagNode[0] == 0 && wetDryFlagNode[1] == 0)
		{
			u[0] = 0;
			u[1] = 0;
		//	Q[2*el+1] = 0;
		//	Q[2*el+2] = 0;


		}

		else
		{
			int npos = wetDryFlagNode[0] + wetDryFlagNode[1];
			//double delu = u[0]*(1-wetDryFlagNode[0]) + u[1]*(1-wetDryFlagNode[1]);
			double delu = u[0]*(1-wetDryFlagNode[0])*(fabs(Q[2*el+1]>VELZERO)) + u[1]*(1-wetDryFlagNode[1])*(fabs(Q[2*el+2])>VELZERO);
		//	double delQ = (1-wetDryFlagNode[0])*Q[2*el+1] + (1-wetDryFlagNode[1])*Q[2*el+2];
			for (int i =0; i < 2; ++i)
			{
				u[i] = wetDryFlagNode[i]*(u[i] + delu/npos);
				//Q[2*el+i] = wetDryFlagNode[i]*(Q[2*el+i]+delQ/npos);
			}
			

		}
		Q[2*el+1] = A[2*el+1]*u[0];
		Q[2*el+2] = A[2*el+2]*u[1];


		if (isnan(u[0]) || isnan(u[1]))
			printf("WD-modified velocity NAN for element %d \n",el);


	}
/*
	printf("A\n");
	for (int i =0; i < NumEl; ++i)
	{
		for (int k =0; k <2; ++k)
		{
			printf("%e\n",A[i*2+k]);
		}
		
	}
	printf("\n");
	printf("Q\n");
	for (int i =0; i < NumEl; ++i)
	{
		for (int k =0; k <2; ++k)
		{
			printf("%e\n",Q[i*2+k]);
		}
	}

	printf("Next\n");
*/
}

#endif









