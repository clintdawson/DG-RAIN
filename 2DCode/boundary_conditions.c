
/****************** Boundary Conditions ************************************/	

#include <math.h>
#include "MeshAttributes.h"

//double *bzeta;
//double *bQn;

void boundary_conditions()
{
	for (int edg = 0; edg<NumEdges; ++edg)
	{
		double qn;
		double zeta_b;

		// if inflow boundary
		if (BdryPrescribed[edg] == 1)
		{
			double x1 = x_coord[EdgtoVert[edg*2]];
			double x2 = x_coord[EdgtoVert[edg*2+1]];
			double y1 = y_coord[EdgtoVert[edg*2]];
			double y2 = y_coord[EdgtoVert[edg*2+1]];
			double edglen = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	
			qn = -0.042/0.914;
			zeta_b = -1000;
			// qn = -140;
			//zeta_b = 12;


			// MacDonald boundary condition
		/*	zeta_b = -1000;
			qn = -2;
		*/

	
		//	zeta_b = -10000;
			//qn = -0.0425/(8*edglen);
//		qn = -30*(edglen);
			//zeta_b = -1000;

		}

		// if outflow boundary
		if (BdryPrescribed[edg] == 2)
		{
		//	qn = 100000;
		//	zeta_b = 0.66;
		
			// MacDonald boundary condition
		/*	zeta_b = 2.87844 + 0.0008;	
			qn = 1000000;
		*/

			//zeta_b = 0.296;
			//qn = 100000000;
			//double zval = 0.5*(z[EdgtoVert[edg*2]] + z[EdgtoVert[edg*2+1]]);
			//zeta_b = 1.69 - zval;
			zeta_b = 0.296;
			qn = 1000000;

		}
		
		if (BdryPrescribed[edg] == 3)
		{	double x1 = x_coord[EdgtoVert[edg*2]];
			double x2 = x_coord[EdgtoVert[edg*2+1]];
			double y1 = y_coord[EdgtoVert[edg*2]];
			double y2 = y_coord[EdgtoVert[edg*2+1]];
			double edglen = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

			//zeta_b = -10000;
			//qn = -0.1275/(8*edglen);
			zeta_b = -10000;
			qn = -0.1275/0.914;
		//qn = -20/(4*edglen);
		}

		bzeta[edg] = zeta_b;
		bQn[edg] = qn;
	}

}



