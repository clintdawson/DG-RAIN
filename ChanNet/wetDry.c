/*****************************************************************************************//**
* @file wetDry.c
*
* This file contains all the functions associated with the wetting and drying treatment.
*
* ****************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ChannelsAndJunctions.h"

#ifdef WDON

/** @cond EXTERN_VARIABLES */
extern const double H0;       
extern const double VELZERO;
/** @endcond */ 

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
* @param[in] Chan pointer to the channel structure that we are currently working on
* **********************************************************************************/
void wetDryStatus1D(struct channel *Chan)
{
	for (int i =0; i < Chan->NumEl; ++i)
	{
		double A0 = Chan->A[2*i+1];
		double A1 = Chan->A[2*i+2];
		double b0 = Chan->b[i];
		double b1 = Chan->b[i+1];
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
void PDop1D(struct channel *Chan)
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
		double h0 = A0/Chan->b[i];
		double h1 = A1/Chan->b[i+1];
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

		A0 = h0*Chan->b[i];
		A1 = h1*Chan->b[i+1];
		Chan->A[2*i+1] = A0;
		Chan->A[2*i+2] = A1;
		if (isnan(A0) || isnan(A1))
		{
			printf("water height negative during wetting and drying treatment for el %d with avgH = %e\n", i, avgH);
			exit(EXIT_FAILURE);
		}
	
	}			

	// redistribute momentum as necessary
	for (int el = 0; el < Chan->NumEl; ++el)
	{
		int wetDryFlagNode[2];
		double u[2];
		for (int i =0; i < 2; ++i)
		{
			double Anode = Chan->A[2*el+i+1];
			double bnode = Chan->b[el+i];
			double hnode = Anode/bnode;
			wetDryFlagNode[i] = BigTheta(hnode - H0);
			
			double Qnode = Chan->Q[2*el+i+1];
			u[i] = Qnode/Anode;
	
		}
		if (wetDryFlagNode[0] == 0 && wetDryFlagNode[1] == 0)
		{
			u[0] = 0;
			u[1] = 0;


		}

		else
		{
			int npos = wetDryFlagNode[0] + wetDryFlagNode[1];
			//double delu = u[0]*(1-wetDryFlagNode[0]) + u[1]*(1-wetDryFlagNode[1]);
			double delu = u[0]*(1-wetDryFlagNode[0])*(fabs(Chan->Q[2*el+1]>VELZERO)) + u[1]*(1-wetDryFlagNode[1])*(fabs(Chan->Q[2*el+2])>VELZERO);
			for (int i =0; i < 2; ++i)
			{
				u[i] = wetDryFlagNode[i]*(u[i] + delu/npos);
			}
			

		}
		Chan->Q[2*el+1] = Chan->A[2*el+1]*u[0];
		Chan->Q[2*el+2] = Chan->A[2*el+2]*u[1];


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

/*********************************************************************************//**
*
* A function that assigns a wet or dry status to each element of a junction and stores
* it in junction::WD. 0 represents a dry element and 1 represents a wet element.
* @param[in] junc pointer to the the junction structure that we are currently working on
* **********************************************************************************/
void wetDryStatus2D(struct junction *junc)
{
	for (int el = 0; el < junc->NumEl; ++el)
	{
		double midZetaEl[3] = {junc->zeta[el*3], junc->zeta[el*3+1], junc->zeta[el*3+2]};
		double zetaEl[3] = {midZetaEl[0] - midZetaEl[1] + midZetaEl[2], midZetaEl[0] + midZetaEl[1] - midZetaEl[2], -midZetaEl[0] + midZetaEl[1] + midZetaEl[2]};

		double zEl[3] = {junc->z[junc->EltoVert[el*3]], junc->z[junc->EltoVert[el*3+1]], junc->z[junc->EltoVert[el*3+2]]};

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
void PDop2D(struct junction *junc)
{
	// redistribute mass and momentum as necessary
	for (int el = 0; el < junc->NumEl; ++el)
	{
		double midZeta0 = junc->zeta[3*el];
		double midZeta1 = junc->zeta[3*el+1];
		double midZeta2 = junc->zeta[3*el+2];
	
		double z0 = junc->z[junc->EltoVert[el*3]]; 
		double z1 = junc->z[junc->EltoVert[el*3+1]];
		double z2 = junc->z[junc->EltoVert[el*3+2]];
	
		double midz0 = 0.5*(z0+z1);
		double midz1 = 0.5*(z1+z2);
		double midz2 = 0.5*(z2+z0);

		double midHt0 = midZeta0 + midz0;
		double midHt1 = midZeta1 + midz1;
		double midHt2 = midZeta2 + midz2;
	
		double Ht0 = midHt0 - midHt1 + midHt2;
		double Ht1 = midHt0 + midHt1 - midHt2;
		double Ht2 = -midHt0 + midHt1 + midHt2;

		double Ht[3] = {Ht0, Ht1, Ht2};

		double avgH = (Ht[0]+Ht[1]+Ht[2])/3;

		// redistribute mass
		if (Ht[0] <= H0 || Ht[1] <= H0 || Ht[2] <= H0)
		{
			if (avgH <= H0)
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
				
				int nextMinHtIndex = 1;
				double nextMinHt = Ht[1];
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

			if (isnan(Ht[0]) || isnan(Ht[1]) || isnan(Ht[2]))
			{
				printf("water height negative during wetting and drying treatment for el %d with avgH = %e\n", el, avgH);
				exit(EXIT_FAILURE);
			}

			midHt0 = 0.5*(Ht[0]+Ht[1]);
			midHt1 = 0.5*(Ht[1]+Ht[2]);
			midHt2 = 0.5*(Ht[2]+Ht[0]);
 
			// Calculate the values at the midpoints
			junc->zeta[3*el] = midHt0 - midz0;
			junc->zeta[3*el+1] = midHt1 - midz1;
			junc->zeta[3*el+2] = midHt2 - midz2;
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
				junc->Qx[3*el+i] = 0;
				junc->Qy[3*el+i] = 0;
	
			}
		}

		else if (NodalWDFlag[0] == 0 || NodalWDFlag[1] == 0 || NodalWDFlag[2] == 0)
		{
			// Store the nodal values of the velocities in each direction
			double u[3];
			double v[3];


			// Calculate new u and v at the nodes
			double midQx0 = junc->Qx[3*el];
			double midQx1 = junc->Qx[3*el+1];
			double midQx2 = junc->Qx[3*el+2];
			
			double midQy0 = junc->Qy[3*el];
			double midQy1 = junc->Qy[3*el+1];
			double midQy2 = junc->Qy[3*el+2];
		
/*			midQx0 = midQx0*(fabs(midQx0) > VELZERO);
			midQx1 = midQx1*(fabs(midQx1) > VELZERO);
			midQx2 = midQx2*(fabs(midQx2) > VELZERO);

			midQy0 = midQy0*(fabs(midQy0) > VELZERO);
			midQy1 = midQy1*(fabs(midQy1) > VELZERO);
			midQy2 = midQy2*(fabs(midQy2) > VELZERO);
*/		

			double Qx_el[3] = {midQx0 - midQx1 + midQx2, midQx0 + midQx1 - midQx2, -midQx0 + midQx1 + midQx2};
			double Qy_el[3] = {midQy0 - midQy1 + midQy2, midQy0 + midQy1 - midQy2, -midQy0 + midQy1 + midQy2};

			for (int i = 0; i < 3; ++i)
			{
				u[i] = Qx_el[i]/Ht[i];
				v[i] = Qy_el[i]/Ht[i];

			}	


			int npos = 0;
			double delu = 0;
			double delv = 0;
			for (int i = 0; i < 3; ++i)
			{
				npos += NodalWDFlag[i];
				delu += u[i]*(1-NodalWDFlag[i])*(fabs(Qx_el[i]) > VELZERO);
				delv += v[i]*(1-NodalWDFlag[i])*(fabs(Qy_el[i]) > VELZERO);
				

			}			

			for (int i = 0; i < 3; ++i)
			{
				u[i] = NodalWDFlag[i]*(u[i]+delu/npos);
				v[i] = NodalWDFlag[i]*(v[i]+delv/npos);

				if (isnan(u[i]) || isnan(v[i]))
				{
					printf("WD-modified velocity NAN for element %d \n",el);
					exit(EXIT_FAILURE);
				}	
			}
			for (int i = 0; i < 3; ++i)
			{
				Qx_el[i] = Ht[i]*u[i];
				Qy_el[i] = Ht[i]*v[i];

			}
	
			for (int i = 0; i < 3; ++i)
			{
				int ind1 = i % 3;
				int ind2 = (i+1) % 3;
				junc->Qx[3*el+i] = 0.5*(Qx_el[ind1] + Qx_el[ind2]);
				junc->Qy[3*el+i] = 0.5*(Qy_el[ind1] + Qy_el[ind2]);

			}

		}

	
	}

}

#endif









