#ifndef ACORE_LIB3_CPP
#define ACORE_LIB3_CPP

#include "ACore_lib3.h"

void GenerMultSignalStatic(int n, ACore * cores, int N, double dt, double * signal)
{
	int i,j,k;
	long double ld;

	for(i=0;i<N;i++)
	{
		ld=0;
		//signal[i]=0;
		for(j=0;j<n;j++)
		{
			for(k=0;k<cores[j].numAmpl;k++)
			{
				ld+=cores[j].A[k]*cos(PI_2*double(k+1)*dt*i*cores[j].baseFr)+cores[j].B[k]*sin(PI_2*double(k+1)*dt*i*cores[j].baseFr);
				//signal[i]+=core->A[k]*cos(PI_2*double(k+1)*dt*i*core->baseFr)+core->B[k]*sin(PI_2*double(k+1)*dt*i*core->baseFr);
			}
		}
		signal[i]=ld;
	}
}
//

void GetAmplMultOld(int n, double *baseFr, int N, 
					double *signal, double *re, double *im, double dt, ACore *cores, int maxAB)
{
	int i,j,k;

	double *sig_gen,*re_gen,*im_gen;
	sig_gen=(double*)malloc(N*sizeof(double));
	re_gen=(double*)malloc(N*sizeof(double));
	im_gen=(double*)malloc(N*sizeof(double));


	for(i=0;i<n;i++)
	{
		cores[i].numAmpl=maxAB;
		cores[i].baseFr=baseFr[i];
		for(j=0;j<maxAB;j++)
		{
			cores[i].A[j]=re[int((j+1)*baseFr[i]*(dt*N)+0.5)];
			cores[i].B[j]=im[int((j+1)*baseFr[i]*(dt*N)+0.5)];
		}
	}

	for(int ii=0;ii<10;ii++)
	{
		GenerMultSignalStatic(n,cores,N,dt,sig_gen);
		fft(N,sig_gen,re_gen,im_gen);
		for(i=0;i<n;i++)
			for(j=0;j<maxAB;j++)
			{
				double z1_re,z1_im,z2_re,z2_im,z2_mod2,c_re,c_im;
				double tmp_re,tmp_im;
				z1_re=re[int((j+1)*baseFr[i]*(dt*N)+0.5)];
				z1_im=im[int((j+1)*baseFr[i]*(dt*N)+0.5)];
				z2_re=re_gen[int((j+1)*baseFr[i]*(dt*N)+0.5)];
				z2_im=im_gen[int((j+1)*baseFr[i]*(dt*N)+0.5)];
				z2_mod2=z2_re*z2_re+z2_im*z2_im;
				c_re=z1_re*z2_re+z1_im*z2_im;	
				c_im=z1_im*z2_re-z1_re*z2_im;
				c_re/=z2_mod2;
				c_im/=z2_mod2;

				tmp_re=cores[i].A[j]*c_re-cores[i].B[j]*c_im;
				tmp_im=cores[i].A[j]*c_im+cores[i].B[j]*c_re;
				cores[i].A[j]=tmp_re;
				cores[i].B[j]=tmp_im;
			}
	}


	free(sig_gen);
	free(re_gen);
	free(im_gen);
}
//

bool SolveSLU(int N, double* A,double *B, double *X)
{
	//A[i][j]
	int i,j,k;
	double coef;

	for(i=0;i<N;i++)
	{
		for(j=i+1;j<N;j++)
		{
			coef=A[i+N*j]/A[i+N*i];
			A[i+N*j]=0;
			for(k=i+1;k<N;k++)
				A[k+N*j]-=coef*A[k+N*i];
			B[j]-=coef*B[i];
		}
	}

	for(i=N-1;i>=0;i--)
	{
		if(A[i+N*i]==0)
			return false;
		X[i]=B[i];
		for(j=i+1;j<N;j++)
			X[i]-=X[j]*A[j+N*i];
		X[i]/=A[i+N*i];
	}

	return true;
}
//

int  GetAmplMult2Std(double *baseFr, int N, double *signal, double *re, double *im, double dt, ACore *cores, int maxAB)
{
	int i,j;

	double l0=int(baseFr[0]*N*dt)/baseFr[0];
	double l1=int(baseFr[1]*N*dt)/baseFr[1];
	cores[0].baseFr=baseFr[0];
	cores[1].baseFr=baseFr[1];
	cores[0].numAmpl=maxAB;
	cores[1].numAmpl=maxAB;
	cores[0].t0=0;
	cores[1].t0=0;

	GetAmplitudesSPEC(N,signal,re,im,dt,baseFr[0],cores[0].A,cores[0].B,maxAB);
	GetAmplitudesSPEC(N,signal,re,im,dt,baseFr[1],cores[1].A,cores[1].B,maxAB);


	double *B=(double*)malloc(4*maxAB*sizeof(double));
	double *X=(double*)malloc(4*maxAB*sizeof(double));
	double *A=(double*)malloc(16*maxAB*maxAB*sizeof(double));

//#define IND_A(i,j,K,L) (A[((i)+(K)*maxAB) + 4*maxAB*((j)+(L)*maxAB)])

	for(i=0;i<maxAB;i++)
	{
		B[i]=cores[0].A[i];
		B[i+maxAB]=cores[0].B[i];
		B[i+2*maxAB]=cores[1].A[i];
		B[i+3*maxAB]=cores[1].B[i];

		for(j=0;j<maxAB;j++)
		{
			A[i+4*maxAB*j]=A[(i+maxAB)+4*maxAB*(j+maxAB)]=A[(i+2*maxAB)+4*maxAB*(j+2*maxAB)]=
				A[(i+3*maxAB)+4*maxAB*(j+3*maxAB)]=(i==j)?1:0;
			A[i+4*maxAB*(j+maxAB)]=A[(i+maxAB)+4*maxAB*j]=
				A[(i+2*maxAB)+4*maxAB*(j+3*maxAB)]=A[(i+3*maxAB)+4*maxAB*(j+2*maxAB)]=0;

			A[(i+2*maxAB)+4*maxAB*j]=alpha((i+1)*baseFr[1],(j+1)*baseFr[0],l0);
			A[(i+3*maxAB)+4*maxAB*j]=gamma((i+1)*baseFr[1],(j+1)*baseFr[0],l0);
			A[(i+2*maxAB)+4*maxAB*(j+maxAB)]=beta((i+1)*baseFr[1],(j+1)*baseFr[0],l0);
			A[(i+3*maxAB)+4*maxAB*(j+maxAB)]=delta((i+1)*baseFr[1],(j+1)*baseFr[0],l0);

			A[i+4*maxAB*(j+2*maxAB)]=alpha((i+1)*baseFr[0],(j+1)*baseFr[1],l1);
			A[(i+maxAB)+4*maxAB*(j+2*maxAB)]=gamma((i+1)*baseFr[0],(j+1)*baseFr[1],l1);
			A[i+4*maxAB*(j+3*maxAB)]=beta((i+1)*baseFr[0],(j+1)*baseFr[1],l1);
			A[(i+maxAB)+4*maxAB*(j+3*maxAB)]=delta((i+1)*baseFr[0],(j+1)*baseFr[1],l1);
		}
	}

/*	FILE * f1=fopen("debug.txt","w");
	for(j=0;j<4*maxAB;j++)
	{
		for(i=0;i<4*maxAB;i++)
			fprintf(f1,"	%.3lf",A[i+j*4*maxAB]);
		fprintf(f1,"	%.3lf",B[j]);
		fprintf(f1,"\n");
	}
	fclose(f1);
//*/

	if(!SolveSLU(4*maxAB,A,B,X))
		return 0;

	for(i=0;i<maxAB;i++)
	{
		cores[0].A[i]=X[i];
		cores[0].B[i]=X[i+maxAB];
		cores[1].A[i]=X[i+2*maxAB];
		cores[1].B[i]=X[i+3*maxAB];
	}

	return 1;
}
//

int Get2Cores(double fr1_st, double fr1_en, double fr2_st, double fr2_en,
			   int N, double *signal, double *re, double *im, double dt, ACore *cores, int maxAB)
{
	int i,j,k;

	double fr[2];
	double fr_d[2];
	
	fr[0]=DichotomyBaseFreq(N,signal,re,im,dt,fr1_st,fr1_en);
//	fr_d[0]=(fr1_en-fr1_st)*0.5;
	fr[1]=DichotomyBaseFreq(N,signal,re,im,dt,fr2_st,fr2_en);
//	fr_d[1]=(fr2_en-fr2_st)*0.5;

	//if(!GetAmplMult2Std(fr,N,signal,re,im,dt,cores,maxAB))
	//	return 0;
	cores[0].baseFr = fr[0];
	cores[1].baseFr = fr[1];
	GetAmplitudesOPT(2,cores,N,signal,dt);

	//*/ ¬ычесть €дра и посчитать как одиночные;
	double * sig  = (double*)malloc(N*sizeof(double));
	double * rest = (double*)malloc(N*sizeof(double));

	fr[0]=0.5*(fr1_st+fr1_en);
	fr_d[0]=0.5*(fr1_en-fr1_st);
	fr[1]=0.5*(fr2_st+fr2_en);
	fr_d[1]=0.5*(fr2_en-fr2_st);

	for(int iter=0; iter<10;iter++)
	{
		for(i=0;i<2;i++)
		{
			for(k=0;k<N;k++)
			{
				rest[k]=signal[k];
			}
			for(j=0;j<2;j++)
			{
				if(j==i)
					continue;
				GenerateSignalStatic(N,sig,dt,&cores[j]);
				for(k=0;k<N;k++)
				{
					rest[k]-=sig[k];
				}
			}
			fr[i]=DichotomyBaseFreq(N,rest,re,im,dt,fr[i]-fr_d[i],fr[i]+fr_d[i]);
			cores[i].baseFr = fr[i];
			cores[i].numAmpl = maxAB;
			//fr_d[i]=...;
		}
	
		/*
		if(!GetAmplMult2Std(fr,N,signal,re,im,dt,cores,maxAB))
		{
			free(rest);
			free(sig);
			return 0;
		}
		//*/
		GetAmplitudesOPT(2,cores,N,signal,dt);
	}

	free(rest);
	free(sig);//*/
	return 1;
}
// 

bool GetAmplOPT(int N, double * signal, int K, double * Ampl, double * basis)
{
	// N - len of signal
	// K - num of basis vectors
	// Ampl - result ampl for vectors
	// basis[N*K] = {{E1},{E2},...,{EK}} basis[n+N*k]

	double * M = (double *) malloc(K*K*sizeof(double));
	if (!M)
	{
		return false;
	}
	double * B = (double *) malloc(K*sizeof(double));
	if (!B)
	{
		free(M);
		return false;
	}

	int i,j,n;
	double tmp;

	for (j = 0; j < K; j++)
	{
		for (i = 0; i < K; i++)
		{
			tmp = 0;
			for ( n = 0; n < N; n++)
			{
				tmp += basis[n+N*i] * basis[n+N*j];
			}
			M[i+K*j] = tmp;
		}
		
		tmp = 0;
		for ( n = 0; n < N; n++)
		{
			tmp += signal[n] * basis[n+N*j];
		}
		B[j] = tmp;
	}

	bool res = SolveSLU(K, M, B, Ampl);

	free(B);
	free(M);

	return res;
}
//

void GetAmplitudesOPT(int nCores, ACore * cores, int N, double *signal, double dt)
{
	int i,j,n;
	int numAmpl=0;
	for (i = 0; i < nCores; i++)
	{
		numAmpl += 2*cores[i].numAmpl;
	}
	double * basis = (double *) malloc(N*numAmpl*sizeof(double));
	if(!basis)
	{
		return;
	}
	double * ampl = (double *) malloc(numAmpl*sizeof(double));
	if (!ampl)
	{
		free(basis);
		return;
	}

	int g_i=0; // global 'i' for abs ampl index in basis & ampl
	for(i = 0; i < nCores; i++)
	{
		//A1,B1,A2,B2,...
		for(j = 0; j < cores[i].numAmpl;j++)
		{
			double fr = (j+1)*cores[i].baseFr;
			for(n=0;n<N;n++)
			{
				basis[n+N*g_i]=cos(PI_2*fr*n*dt);
				basis[n+N*(g_i+1)]=sin(PI_2*fr*n*dt);
			}
			g_i+=2;
		}
	}

	if (!GetAmplOPT(N,signal,numAmpl,ampl,basis))
	{
		return;
	}


	g_i=0;
	for(i = 0; i < nCores; i++)
	{
		for(j = 0; j < cores[i].numAmpl;j++)
		{
			cores[i].A[j] = ampl[g_i];
			cores[i].B[j] = ampl[g_i+1];
			g_i += 2;
		}
	}

	free(ampl);
	free(basis);
}

//
void GetAmplitudesOPT(int N, double *signal, double dt, 
					  double baseFr, double * A, double * B, int maxAB)
{
	int i,j,n;
	int numAmpl=2*maxAB+1;

	double * basis = (double *) malloc(N*numAmpl*sizeof(double));
	if(!basis)
	{
		return;
	}
	double * ampl = (double *) malloc(numAmpl*sizeof(double));
	if (!ampl)
	{
		free(basis);
		return;
	}


	for(j = 0; j < maxAB;j++)
	{
		double fr = (j+1)*baseFr;
		for(n=0;n<N;n++)
		{
			basis[n+N*(2*j)]=cos(PI_2*fr*n*dt);
			basis[n+N*(2*j+1)]=sin(PI_2*fr*n*dt);
		}
	}
	for(n=0;n<N;n++)
	{
		basis[n+N*(numAmpl-1)]=1;
	}

	if (!GetAmplOPT(N,signal,numAmpl,ampl,basis))
	{
		return;
	}


	for(j = 0; j < maxAB;j++)
	{
		A[j] = ampl[2*j];
		B[j] = ampl[2*j+1];
	}

	free(ampl);
	free(basis);
}

#endif