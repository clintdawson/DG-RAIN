# include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "MeshAttributes.h"


int BigTheta(double a)
{
	if  (a > 0)
		return(1);
	else
		return(0);

}

void wetDryStatus2D()
{
	for (int el = 0; el < NumEl; ++el)
	{
		double midZetaEl[3] = {zeta[el*3], zeta[el*3+1], zeta[el*3+2]};
		double zetaEl[3] = {midZetaEl[0] - midZetaEl[1] + midZetaEl[2], midZetaEl[0] + midZetaEl[1] - midZetaEl[2], -midZetaEl[0] + midZetaEl[1] + midZetaEl[2]};

		double zEl[3] = {z[EltoVert[el*3]], z[EltoVert[el*3+1]], z[EltoVert[el*3+2]]};

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

		if (WD[el] == -1)		// if uninitialized
			WD[el] = BigTheta(avgH - H0)*BigTheta(ZetaMax - (H0 - zmin));
		else
			WD[el] = WD[el]*BigTheta(avgH - H0) + (1-WD[el])*BigTheta(avgH - H0)*BigTheta(ZetaMax - (H0 - zmin));


	}

}

// post processing to ensure that the water depth remains positive
void PDop2D()
{
	// redistribute mass and momentum as necessary
	for (int el = 0; el < NumEl; ++el)
	{
		double midZeta0 = zeta[3*el];
		double midZeta1 = zeta[3*el+1];
		double midZeta2 = zeta[3*el+2];
	
		double z0 = z[EltoVert[el*3]]; 
		double z1 = z[EltoVert[el*3+1]];
		double z2 = z[EltoVert[el*3+2]];
	
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
			zeta[3*el] = midHt0 - midz0;
			zeta[3*el+1] = midHt1 - midz1;
			zeta[3*el+2] = midHt2 - midz2;
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
				Qx[3*el+i] = 0;
				Qy[3*el+i] = 0;
	
			}
		}

		else if (NodalWDFlag[0] == 0 || NodalWDFlag[1] == 0 || NodalWDFlag[2] == 0)
		{
			// Store the nodal values of the velocities in each direction
			double u[3];
			double v[3];


			// Calculate new u and v at the nodes
			double midQx0 = Qx[3*el];
			double midQx1 = Qx[3*el+1];
			double midQx2 = Qx[3*el+2];
			
			double midQy0 = Qy[3*el];
			double midQy1 = Qy[3*el+1];
			double midQy2 = Qy[3*el+2];
		
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
				Qx[3*el+i] = 0.5*(Qx_el[ind1] + Qx_el[ind2]);
				Qy[3*el+i] = 0.5*(Qy_el[ind1] + Qy_el[ind2]);

			}

		}

	
	}

}

