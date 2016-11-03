#ifndef CONSTITUTIVE_RELATIONS

#define CONSTITUTIVE_RELATIONS

extern double getH(double A, double b, double m1, double m2);
extern double getI1(double A, double b, double m1, double m2);
extern double getI2(double A, double b, double db, double m1, double dm1, double m2, double dm2);
extern double getS_f(double A, double Q, double b, double m1, double m2, double n);

#endif
