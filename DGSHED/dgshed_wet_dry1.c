/*****************************************************************************************//**
* @file wetDry.c
*
* This file contains all the functions associated with the wetting and drying treatment.
*
* ****************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Watershed.h"
#include "Globals.h"
#include "constitutive_equations.h"

#ifdef WDON

/***********************************************************************************//**
* 
* A function that checks to see if a value is positive.
* @param a a double value to be tested
* @return 1 if positive, 0 if negative
* **************************************************************************************/
int BigTheta(double a)
{
	if (a > 0)
		return(1);
	else
		return(0);
}

/*********************************************************************************//**
*
* A function that assigns a wet or dry status to each element of a channel and stores
* it in channel::WD. 0 represents a dry element and 1 represents a wet element.
* WARNING: FILE HAS NOT BEEN MODIFIED FOR P !=1. ONLY TURN ON WHEN P = 1
* @param[in] Chan pointer to the channel structure that we are currently working on
* **********************************************************************************/
void wetDryStatusChan(struct channel *Chan)
{
	for (int i =0; i < Chan->NumEl; ++i)
	{
		double A0 = Chan->A[2*i+1];
		double A1 = Chan->A[2*i+2];
		double b0 = Chan->b[i];
		double b1 = Chan->b[i+1];
		double m10 = Chan->m1[i];
		double m11 = Chan->m1[i+1];
		double m20 = Chan->m2[i];
		double m21 = Chan->m2[i+1];
		double h0 = getH(A0, b0, m10, m20);
		double h1 = getH(A1, b1, m11, m21);
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
		
		double ZetaMax = maxHeight - Chan->z[maxHeightIndex];
		double z0 = Chan->z[i];
		double z1 = Chan->z[i+1];
		double zmin = z0*(z0 < z1) + z1*(z1 < z0); 

		if (Chan->WD[i] == -1) 	// if uninitialized
			Chan->WD[i] = BigTheta(avgH-H0)*BigTheta(ZetaMax - (H0-zmin));
		else
			Chan->WD[i] = Chan->WD[i]*BigTheta(avgH - H0)+(1-Chan->WD[i])*BigTheta(avgH-H0)*BigTheta(ZetaMax - (H0 - zmin));

	}


}

/**************************************************************************************//**
*
* This function does some post processing to ensure that the water depth remains positive
* in the 1-D channels
* @param[in] Chan pointer to the channel structure that we are currently working on
*
* *****************************************************************************************/
void PDopChan(struct channel *Chan)
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
	for (int i =0; i < Chan->NumEl; ++i)
	{
		double A0 = Chan->A[2*i+1];
		double A1 = Chan->A[2*i+2];
		double b0 = Chan->b[i];
		double b1 = Chan->b[i+1];
		double m10 = Chan->m1[i];
		double m11 = Chan->m1[i+1];
		double m20 = Chan->m2[i];
		double m21 = Chan->m2[i+1];
		double h0 = getH(A0, b0, m10, m20);
		double h1 = getH(A1, b1, m11, m21);
		double avgH = (h0+h1)/2;
		double ht[2] = {h0, h1};
		
		if (h0 >= H0 && h1 >= H0)
		{	
			ht[0] = h0;
			ht[1] = h1;
		}		
	
		else if (avgH < H0)
		{
			ht[0] = avgH;
			ht[1] = avgH;
		}

		else
		{
			int minHeightIndex;
			double minHeight;
			int remInd;
			double remHt;
			if (h0 < h1)
			{
				minHeightIndex = 0;
				minHeight = h0;
				remInd = 1;
				remHt = h1;
			}
			else
			{
				minHeightIndex = 1;
				minHeight = h1;
				remInd = 0;
				remHt = h0;
			}
		
			ht[minHeightIndex] = H0;
			ht[remInd] = fmax(H0, remHt - (H0 - minHeight));

		}

		A0 = ht[0]*b0 + 0.5*m10*ht[0]*ht[0] + 0.5*m20*ht[0]*ht[0];
		A1 = ht[1]*b1 + 0.5*m11*ht[1]*ht[1] + 0.5*m21*ht[1]*ht[1];
		Chan->A[2*i+1] = A0;
		Chan->A[2*i+2] = A1;
		
		if (A0 < 0 || A1 < 0)
		{
			printf("water height negative during wetting and drying treatment for el %d with avgH = %e\n", i, avgH);
			exit(EXIT_FAILURE);
		}
	
	}			

	// redistribute momentum as necessary
	for (int el = 0; el < Chan->NumEl; ++el)
	{
		int wetDryFlagNode[2];
		double Q[2];
		for (int i =0; i < 2; ++i)
		{
			double Anode = Chan->A[2*el+ i+1];
			double bnode = Chan->b[el+i];
			double m1node = Chan->m1[el+i];
			double m2node = Chan->m2[el+i];
			double hnode = getH(Anode, bnode, m1node, m2node);

			wetDryFlagNode[i] = BigTheta(hnode - H0);

			Q[i] = Chan->Q[2*el+i+1];
	
		}
		
		if (wetDryFlagNode[0] == 0 && wetDryFlagNode[1] == 0)
		{
			Chan->Q[2*el+1] = 0;
			Chan->Q[2*el+2] = 0;

		}

		else
		{
			int npos = wetDryFlagNode[0] + wetDryFlagNode[1];
			double delQ = Q[0]*(1-wetDryFlagNode[0]) + Q[1]*(1-wetDryFlagNode[1]);
			//double delu = u[0]*(1-wetDryFlagNode[0])*(fabs(Chan->Q[2*el+1]>VELZERO)) + u[1]*(1-wetDryFlagNode[1])*(fabs(Chan->Q[2*el+2])>VELZERO);
			for (int i =0; i < 2; ++i)
			{
				Q[i] = wetDryFlagNode[i]*(Q[i] + delQ/npos);
			}
		
			if (isnan(Q[0]) || isnan(Q[1]))
			{
				printf("npos = %d\n", npos);
				printf("wetDryFlagNode[0] = %d wetDryFlagNode[1] = %d\n", wetDryFlagNode[0], wetDryFlagNode[1]);
				printf("Q[0] = %lf Q[1] = %lf\n", Chan->Q[2*el+1], Chan->Q[2*el+2]);
				printf("WD-modified velocity NAN for element %d \n",el);

			}		

		}
		Chan->Q[2*el+1] = Q[0];
		Chan->Q[2*el+2] = Q[1];
		

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

/*********************************************************************************//**
*
* A function that assigns a wet or dry status to each element of a TwoDRegion and stores
* it in TwoDRegion::WD. 0 represents a dry element and 1 represents a wet element.
* @param[in] currRegion pointer to the the TwoDRegion structure that we are currently working on
* **********************************************************************************/
// only works for P=1 elements
void wetDryStatus2D(struct TwoDRegion *junc)
{
	for (int el = 0; el < junc->NumEl; ++el)
	{
		double zetaEl[3] = {junc->zeta[el][0], junc->zeta[el][1], junc->zeta[el][2]};
		double zEl[3] = {junc->NodalZ[el][0], junc->NodalZ[el][1], junc->NodalZ[el][2]};
	
		double height[3] = {zetaEl[0] + zEl[0], zetaEl[1] + zEl[1], zetaEl[2] + zEl[2]};

		double avgH = (height[0]+height[1]+height[2])/3;

		int maxHeightIndex = 0;
		double maxHeight = height[0];

		for (int i =1; i < 3; ++i)
		{
			if (height[i] > maxHeight)
			{
				maxHeightIndex = i;
				maxHeight = height[i];
			} 
		}

		double ZetaMax = zetaEl[maxHeightIndex];
		
		double zmin = zEl[0];
		for (int i = 1; i < 3; ++i)
		{
			if (zEl[i] < zmin)
			{
				zmin = zEl[i];
			} 

		}

		if (junc->WD[el] == -1)		// if uninitialized
			junc->WD[el] = BigTheta(avgH - H0)*BigTheta(ZetaMax - (H0 - zmin));
		else
			junc->WD[el] = junc->WD[el]*BigTheta(avgH - H0) + (1-junc->WD[el])*BigTheta(avgH - H0)*BigTheta(ZetaMax - (H0 - zmin));

	}

}

/**************************************************************************************//**
*
* This function does some post processing to ensure that the water depth remains positive
* in the 2-D junctions
* @param[in] junc pointer to the junction structure that we are currently working on
*
* *****************************************************************************************/
void PDop2D(struct TwoDRegion *junc, double time)
{
	// redistribute mass and momentum as necessary
	for (int el = 0; el < junc->NumEl; ++el)
	{
		
		double zeta[3] = {junc->zeta[el][0], junc->zeta[el][1], junc->zeta[el][2]};
		double z[3] = {junc->NodalZ[el][0], junc->NodalZ[el][1], junc->NodalZ[el][2]};
		double avgH = 0;
		double Ht[3];
		for (int s = 0; s < 3; s++)
		{
			Ht[s] = zeta[s] + z[s];
			avgH += Ht[s];
		}
		avgH /= 3;

		// redistribute mass
		if (Ht[0] < H0 || Ht[1] < H0 || Ht[2] < H0)
		{
			if (avgH < H0)
			{
				for (int i =0; i < 3; ++i)
					Ht[i] = avgH;	
			}

			else
			{
				int minHtIndex = 0;
				double minHt = Ht[0];
				for (int i = 1; i < 3; ++i)
				{
					if (Ht[i] < minHt)
					{
						minHtIndex = i;
						minHt = Ht[i];
					} 

				}
				
				int nextMinHtIndex = (minHtIndex+1)%3;
				double nextMinHt = Ht[nextMinHtIndex];
				for (int i = 0; i < 3; ++i)
				{
					if (Ht[i] < nextMinHt && i != minHtIndex)
					{
						nextMinHtIndex = i;
						nextMinHt = Ht[i];
					}
				}

				int remIndex = 2;
				for (int i =0; i <3; ++i)
				{
					if (i != minHtIndex && i != nextMinHtIndex)
					{
						remIndex = i;
						break;
					}
				}
				
				Ht[minHtIndex] = H0;
				Ht[nextMinHtIndex] = fmax(H0, nextMinHt - 0.5*(Ht[minHtIndex] - minHt));
				Ht[remIndex] = Ht[remIndex] - (Ht[minHtIndex] - minHt) - (Ht[nextMinHtIndex] - nextMinHt); 

			}

			if (Ht[0] < 0 || Ht[1] < 0 || Ht[2] < 0)
			{
				printf("water height negative during wetting and drying treatment for el %d with avgH = %e at time %e for domain type %d\n", el, avgH, time, junc->type);
				printf("zeta1 = %lf, zeta2 = %lf, zeta3 = %lf\n", zeta[0], zeta[1], zeta[2]);
				printf("z1 = %lf, z2 = %lf, z3 = %lf\n", z[0], z[1], z[2]);
				exit(EXIT_FAILURE);
			}
		
			// Assign the nodal values
			junc->zeta[el][0] = Ht[0] - z[0];
			junc->zeta[el][1] = Ht[1] - z[1];
			junc->zeta[el][2] = Ht[2] - z[2];
		}

		int NodalWDFlag[3];
		
		// Store the wet and dry nodes
		for (int i = 0; i < 3; ++i)
		{
			NodalWDFlag[i] = (Ht[i] > H0);

		}		
		
		// redistribute momentum if necessary

		
		if (NodalWDFlag[0] == 0 && NodalWDFlag[1] == 0 && NodalWDFlag[2] == 0)
		{
			for (int i = 0; i < 3; ++i)
			{
				junc->Qx[el][i] = 0;
				junc->Qy[el][i] = 0;
	
			}
		}

		else if (NodalWDFlag[0] == 0 || NodalWDFlag[1] == 0 || NodalWDFlag[2] == 0)
		{
			// Store the nodal values of the velocities in each direction
			double u[3];
			double v[3];

			double Qx[3] = {junc->Qx[el][0], junc->Qx[el][1], junc->Qx[el][2]};
			double Qy[3] = {junc->Qy[el][0], junc->Qy[el][1], junc->Qy[el][2]};


			int npos = 0;
			double delQx = 0;
			double delQy = 0;
			for (int i = 0; i < 3; ++i)
			{
				npos += NodalWDFlag[i];
				delQx += Qx[i]*(1-NodalWDFlag[i]) * (fabs(Qx[i]) > VELZERO);
				delQy += Qy[i]*(1-NodalWDFlag[i]) * (fabs(Qy[i]) > VELZERO);

			}			

			for (int i = 0; i < 3; ++i)
			{
				if (NodalWDFlag[i] > 0)
				{
					if (Qx[i]*delQx > 0) 
						Qx[i] = Qx[i]+delQx/npos;

					if (Qy[i]*delQy > 0)
						Qy[i] = Qy[i]+delQy/npos;
					
				}
				if (isnan(Qx[i]) || isnan(Qy[i]))
				{
					printf("WD-modified flux NAN for element %d \n",el);
					printf("Ht[0] = %lf, Ht[1] = %lf, Ht[2] = %lf\n", Ht[0], Ht[1], Ht[2]);
					printf("Qx[0] = %lf, Qx[1] = %lf, Qx[2] = %lf\n", Qx[0], Qx[1], Qx[2]);
					printf("Qy[0] = %lf, Qy[1] = %lf, Qy[2] = %lf\n", Qy[0], Qy[1], Qx[2]);
					exit(EXIT_FAILURE);
				}	
			}
	
			for (int i = 0; i < 3; ++i)
			{
				junc->Qx[el][i] = Qx[i];
				junc->Qy[el][i] = Qy[i];

			}


			//for (int i = 0; i < 3; ++i)
			//{
			//	u[i] = Qx[i]/Ht[i];
			//	v[i] = Qy[i]/Ht[i];
			//}
			//int npos = 0;
			//double delu = 0;
			//double delv = 0;
			//for (int i = 0; i < 3; ++i)
			//{
			//	npos += NodalWDFlag[i];
			//	delu += u[i]*(1-NodalWDFlag[i])*(fabs(Qx[i]) > VELZERO);
			//	delv += v[i]*(1-NodalWDFlag[i])*(fabs(Qy[i]) > VELZERO);

			//}			

			//for (int i = 0; i < 3; ++i)
			//{
			//	u[i] = NodalWDFlag[i]*(u[i]+delu/npos);
			//	v[i] = NodalWDFlag[i]*(v[i]+delv/npos);

			//	if (isnan(u[i]) || isnan(v[i]))
			//	{
			//		printf("WD-modified velocity NAN for element %d \n",el);
			//		printf("Ht[0] = %lf, Ht[1] = %lf, Ht[2] = %lf\n", Ht[0], Ht[1], Ht[2]);
			//		printf("Qx[0] = %lf, Qx[1] = %lf, Qx[2] = %lf\n", Qx[0], Qx[1], Qx[2]);
			//		printf("Qy[0] = %lf, Qy[1] = %lf, Qy[2] = %lf\n", Qy[0], Qy[1], Qx[2]);
			//		exit(EXIT_FAILURE);
			//	}	
			//}
			//for (int i = 0; i < 3; ++i)
			//{
			//	Qx[i] = Ht[i]*u[i];
			//	Qy[i] = Ht[i]*v[i];

			//}
	
			//for (int i = 0; i < 3; ++i)
			//{
			//	junc->Qx[el][i] = Qx[i];
			//	junc->Qy[el][i] = Qy[i];

			//}

		}

	
	}

}

#endif

