#include <stdio.h>
#include <math.h>
#include "mathfunctions.h"
#include "globals.h"
#include "constitutive_equations.h"

extern double getI1(double A, double b, double m1, double m2);

void RoeFlux(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, 
	double m1, double m2, double localG)
{
	double I1_L, I1_R;
	double F1_L, F1_R;
	double F2_L, F2_R;
	double c_L, c_R;
	double u_L, u_R;
	double lambda1_L, lambda1_R, lambda2_L, lambda2_R;
	double uhat, chat, alpha1, alpha2;
	double lambdahat1, lambdahat2;
	double epsilon1, epsilon2;

	double h_L = getH(A_L,b,m1,m2);
	double h_R = getH(A_R, b, m1, m2);
	
	// width at the free surface
	double B_L = b+m1*h_L+m2*h_L;
	double B_R = b+m1*h_R+m2*h_R;

	I1_L = getI1(A_L, b, m1, m2);
	I1_R = getI1(A_R, b, m1, m2);

	F1_L = Q_L;
	F1_R = Q_R;

	F2_L = pow(Q_L,2)/A_L + localG*I1_L;
	F2_R = pow(Q_R,2)/A_R + localG*I1_R;
		
	c_L = sqrt(g*A_L/B_L);
	c_R = sqrt(g*A_R/B_R);
		
	u_L = Q_L/A_L;
	u_R = Q_R/A_R;
	
	lambda1_L = u_L + c_L;
	lambda1_R = u_R + c_R;
	lambda2_L = u_L - c_L;
	lambda2_R = u_R - c_R;
		
	uhat = (Q_L/sqrt(A_L) + Q_R/sqrt(A_R))/(sqrt(A_L) + sqrt(A_R));
	chat = sqrt(g/2*(A_L/B_L + A_R/B_R));
		
	lambdahat1 = uhat + chat;
	lambdahat2 = uhat - chat;

	epsilon1 = max(0, (lambdahat1 - lambda1_L));
	epsilon1 = max(epsilon1, (lambda1_R - lambdahat1));
	lambdahat1 = max(epsilon1,fabs(lambdahat1));
		
	epsilon2 = max(0, (lambdahat2 - lambda2_L));
	epsilon2 = max(epsilon2, (lambda2_R - lambdahat2));
	lambdahat2 = max(epsilon2,fabs(lambdahat2));

	// Eigenvectors
	double R1[2] = {1, uhat + chat};
	double R2[2] = {1, uhat - chat};
	double Rhat[2][2] = {{1,1},{uhat+chat,uhat-chat}};
	double dtm = 1/(Rhat[0][0]*Rhat[1][1]-Rhat[0][1]*Rhat[1][0]);
	double Rhatinv[2][2] = {{dtm*Rhat[1][1], -dtm*Rhat[0][1]},{-dtm*Rhat[1][0], dtm*Rhat[0][0]}};

	double RoeJac[2][2];
	RoeJac[0][0] = lambdahat1*Rhat[0][0]*Rhatinv[0][0]+lambdahat2*Rhat[0][1]*Rhatinv[1][0]; 
	RoeJac[0][1] = lambdahat1*Rhat[0][0]*Rhatinv[0][1]+lambdahat2*Rhat[0][1]*Rhatinv[1][1]; 
	RoeJac[1][0] = lambdahat1*Rhat[1][0]*Rhatinv[0][0]+lambdahat2*Rhat[1][1]*Rhatinv[1][0]; 
	RoeJac[1][1] = lambdahat1*Rhat[1][0]*Rhatinv[0][1]+lambdahat2*Rhat[1][1]*Rhatinv[1][1]; 

	double jump[2] = {A_L - A_R, Q_L - Q_R};

	Fhat[0] = 0.5*(F1_L + F1_R) + 0.5*(RoeJac[0][0]*jump[0]+RoeJac[0][1]*jump[1]);
	Fhat[1] = 0.5*(F2_L + F2_R) + 0.5*(RoeJac[1][0]*jump[0]+RoeJac[1][1]*jump[1]);

}

void LF(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double m1, double m2, double localG)
{
	double I1_L, I1_R;
	double F1_L, F1_R;
	double F2_L, F2_R;
	double c_L, c_R;
	double u_L, u_R;

	double h_L = getH(A_L,b, m1, m2);
	double h_R = getH(A_R,b, m1, m2);
	I1_L = getI1(A_L, b, m1, m2);
	I1_R = getI1(A_R, b, m1, m2);

	F1_L = Q_L;
	F1_R = Q_R;

	F2_L = pow(Q_L,2)/A_L + I1_L;
	F2_R = pow(Q_R,2)/A_R + I1_R;
	//F2_L = pow(Q_L,2)/A_L + localG*I1_L;
	//F2_R = pow(Q_R,2)/A_R + localG*I1_R;

	// free surface width
	double B_L = b+m1*h_L+m2*h_L;
	double B_R = b+m1*h_R + m2*h_R;

	c_L = sqrt(g*A_L/B_L);
	c_R = sqrt(g*A_R/B_R);
		
	u_L = Q_L/A_L;
	u_R = Q_R/A_R;

	double C = max(fabs(u_L + c_L), fabs(u_R+c_R));
	C = max(C, fabs(u_L - c_L));
	C = max(C, fabs(u_R - c_R));
	double jump[2] = {A_L - A_R, Q_L - Q_R};
	Fhat[0] = 0.5*(F1_L+F1_R) + 0.5*C*jump[0];
	Fhat[1] = 0.5*(F2_L+F2_R) + 0.5*C*jump[1];
	

}


