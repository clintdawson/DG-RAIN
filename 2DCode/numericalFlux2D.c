#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"
#include "MeshAttributes.h"

double RoeFlux(double *Fhatdotn, double zeta_in, double zeta_ex, double Qx_in, double Qx_ex,
	 double Qy_in, double Qy_ex, double z_edge, double nx, double ny, double localG)
{
	double H_in = zeta_in + z_edge;
	double H_ex = zeta_ex + z_edge;

	double u_in = Qx_in/H_in;
	double u_ex = Qx_ex/H_ex;

	double v_in = Qy_in/H_in;
	double v_ex = Qy_ex/H_ex;

	double jump[3] = {zeta_in - zeta_ex, Qx_in - Qx_ex, Qy_in - Qy_ex};

	double F1x_in = H_in*u_in;
	double F1x_ex = H_ex*u_ex;

	double F1y_in = H_in*v_in;
	double F1y_ex = H_ex*v_ex;
	 
	double F2x_in = H_in*u_in*u_in + 0.5*localG*(H_in*H_in-z_edge*z_edge);
	double F2x_ex = H_ex*u_ex*u_ex + 0.5*localG*(H_ex*H_ex-z_edge*z_edge);

	double F2y_in = H_in*u_in*v_in;
	double F2y_ex = H_ex*u_ex*v_ex;

	double F3x_in = H_in*u_in*v_in;
	double F3x_ex = H_ex*u_ex*v_ex;

	double F3y_in = H_in*v_in*v_in + 0.5*localG*(H_in*H_in-z_edge*z_edge);
 	double F3y_ex = H_ex*v_ex*v_ex + 0.5*localG*(H_ex*H_ex-z_edge*z_edge);

	double Fn_in[3] = {F1x_in*nx+F1y_in*ny, F2x_in*nx+F2y_in*ny, F3x_in*nx+F3y_in*ny};
		
	double Fn_ex[3] = {F1x_ex*nx+F1y_ex*ny, F2x_ex*nx+F2y_ex*ny, F3x_ex*nx+F3y_ex*ny};

		

	// Roe's averages
	double Hhat = 0.5*(H_in+H_ex);
	double uhat = (u_in*sqrt(H_in)+u_ex*sqrt(H_ex))/(sqrt(H_in)+sqrt(H_ex));
	double vhat = (v_in*sqrt(H_in)+v_ex*sqrt(H_ex))/(sqrt(H_in)+sqrt(H_ex));

	double Rhat[3][3] = {{1, 0, 1},{uhat-sqrt(g*Hhat)*nx, -ny, uhat+sqrt(g*Hhat)*nx},{vhat-sqrt(g*Hhat)*ny, nx,vhat+sqrt(g*Hhat)*ny}};
	double Rhatinv[3][3];
	double denom = -2*sqrt(g*Hhat);
	Rhatinv[0][0] = (-sqrt(g*Hhat)-nx*uhat-ny*vhat)/denom;
	Rhatinv[0][1] = nx/denom;
	Rhatinv[0][2] = ny/denom;
	Rhatinv[1][0] = (-2*sqrt(g*Hhat)*ny*uhat+2*sqrt(g*Hhat)*nx*vhat)/denom;
	Rhatinv[1][1] = 2*sqrt(g*Hhat)*ny/denom;
	Rhatinv[1][2] = -2*sqrt(g*Hhat)*nx/denom;
	Rhatinv[2][0] = (-sqrt(g*Hhat)+nx*uhat+ny*vhat)/denom;
	Rhatinv[2][1] = -nx/denom;
	Rhatinv[2][2] = -ny/denom;

	double lambda1_in = u_in*nx+v_in*ny-sqrt(g*H_in);
	double lambda2_in = u_in*nx+v_in*ny;
	double lambda3_in = u_in*nx+v_in*ny+sqrt(g*H_in);

	double lambda1_ex = u_ex*nx+v_ex*ny-sqrt(g*H_ex);
	double lambda2_ex = u_ex*nx+v_ex*ny;
	double lambda3_ex = u_ex*nx+v_ex*ny+sqrt(g*H_ex);

	double lambdahat1 = uhat*nx+vhat*ny-sqrt(g*Hhat);
	double lambdahat2 = uhat*nx+vhat*ny;
	double lambdahat3 = uhat*nx+vhat*ny+sqrt(g*Hhat);

	double epsilon1 = fmax(0,(lambdahat1-lambda1_in));
	epsilon1 = fmax(epsilon1, (lambda1_ex-lambdahat1));
	lambdahat1 = fmax(epsilon1, fabs(lambdahat1));

	double epsilon2 = fmax(0,(lambdahat2-lambda2_in));
	epsilon2 = fmax(epsilon2, (lambda2_ex-lambdahat2));
	lambdahat2 = fmax(epsilon2, fabs(lambdahat2));

	double epsilon3 = fmax(0,(lambdahat3-lambda3_in));
	epsilon3 = fmax(epsilon3, (lambda3_ex-lambdahat3));
	lambdahat3 = fmax(epsilon3, fabs(lambdahat3));

	double AbsLambdaHat[3][3] = {{lambdahat1,0,0},{0,lambdahat2,0},{0,0,lambdahat3}};
	
	double prod1[3];
	double prod2[3];
	double R_Lam_Rinv_jump[3]; 

	MatrixVectorMultiply(&Rhatinv[0][0],jump,prod1,3,3);
	MatrixVectorMultiply(&AbsLambdaHat[0][0],prod1,prod2,3,3);
	MatrixVectorMultiply(&Rhat[0][0],prod2,R_Lam_Rinv_jump,3,3);

	for (int k=0; k < 3; ++k)
	{
		Fhatdotn[k] = 0.5*(Fn_in[k]+Fn_ex[k]) + 0.5*R_Lam_Rinv_jump[k];
	} 		

	double current_max_lam = max(fabs(lambda1_in), fabs(lambda1_ex));
	current_max_lam = max(current_max_lam, fabs(lambda2_in));
	current_max_lam = max(current_max_lam, fabs(lambda2_ex));
	current_max_lam = max(current_max_lam, fabs(lambda3_in));
	current_max_lam = max(current_max_lam, fabs(lambda3_ex));

	return current_max_lam;
}
