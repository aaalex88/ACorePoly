#ifndef ACORE_LIB3_H
#define ACORE_LIB3_H

#include "Acore_lib2.h"


void GenerMultSignalStatic(int n, ACore * cores, int N, double dt, double * signal);

void GetAmplMultOld(int n, double *baseFr, int N, double *signal, double *re, double *im, double dt, ACore *cores, int maxAB);

bool SolveSLU(int N, double* A,double *B, double *X);
int  GetAmplMult2Std(double *baseFr, int N, double *signal, double *re, double *im, double dt, ACore *cores, int maxAB);
int Get2Cores(double fr1_st, double fr1_en, double fr2_st, double fr2_en,
			   int N, double *signal, double *re, double *im, double dt, ACore *cores, int maxAB);

double alpha(double w0, double w1, double l)
{
	if(w0*w0==w1*w1)
		return 0;
	return (1.0/(PI*l))*(w0/(w0*w0-w1*w1))*(sin(PI_2*w0*l));
}

double beta(double w0, double w1, double l)
{
	if(w0*w0==w1*w1)
		return 0;
	return (1.0/(PI*l))*(w1/(w0*w0-w1*w1))*(cos(PI_2*w0*l)-1);
}

double gamma(double w0, double w1, double l)
{
	if(w0*w0==w1*w1)
		return 0;
	return (1.0/(PI*l))*(w0/(w0*w0-w1*w1))*(1-cos(PI_2*w0*l));
}

double delta(double w0, double w1, double l)
{
	if(w0*w0==w1*w1)
		return 0;
	return (1.0/(PI*l))*(w1/(w0*w0-w1*w1))*(sin(PI_2*w0*l));
}


bool GetAmplOPT(int N, double * signal, int K, double * Ampl, double * basis);
void GetAmplitudesOPT(int nCores, ACore * cores, int N, double *signal, double dt);
void GetAmplitudesOPT(int N, double *signal, double dt, 
					  double baseFr, double * A, double * B, int maxAB=AMPL_MAX_AB);




#endif