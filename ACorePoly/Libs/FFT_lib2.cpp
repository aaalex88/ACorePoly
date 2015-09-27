#ifndef FFT_LIB2_CPP
#define FFT_LIB2_CPP

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "FFT_lib2.h"


const unsigned char reverse256[]= {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0,
    0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8,
    0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4,
    0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC,
    0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2,
    0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA,
    0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6,
    0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE,
    0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1,
    0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9,
    0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5,
    0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED,
    0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3,
    0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB,
    0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7,
    0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF,
    0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF,
};



void ampl_count(unsigned int N, double * re, double * im, double * x)
{
	unsigned int i;
	for(i=0;i<N;i++)
		x[i]=sqrt(re[i]*re[i] + im[i]*im[i]);
}


void fft_N2_real(unsigned int N, double * x, double * re, double * im) // N = 2^T
{
//	memset((void*)im,0,N*sizeof(double));

	double S;
	unsigned int I,J;
	unsigned char *Ic;
    unsigned char *Jc;

	unsigned int i,id2,k,m,mpNd2;
	long double WX,WY,TempX,TempY;

	Ic = (unsigned char*) &I;
	Jc = (unsigned char*) &J;

	for(i=0;i<N;i++)
	{
		re[i]=x[i];//memcpy
		im[i]=0.0;//memset
	}



	int T=-1;
	int _N=N;
	while(_N)
	{
		_N=_N>>1;
		T++;
	}
	for(I=1;I<N-1;I++)
	{
		Jc[0] = reverse256[Ic[3]];
	    Jc[1] = reverse256[Ic[2]];
	    Jc[2] = reverse256[Ic[1]];
	    Jc[3] = reverse256[Ic[0]];
	    J >>= (32 - T);
	    if (I < J)
	    {
	        S = re[I];
	        re[I] = re[J];
	        re[J] = S;
	    }
	}
	re[0]=x[0];
	re[N-1]=x[N-1];



	for(i=2,id2=1;i<=N;id2=i,i+=i)
	{
		double _i=i;
		for(k=0;k<id2;k++)
		{
			WX=cos(-PI_2*k/_i);
			WY=sin(-PI_2*k/_i);
			for(m=k;m<N;m+=i)
			{
				mpNd2=m+id2;
				TempX=WX*re[mpNd2]-WY*im[mpNd2];
				TempY=WY*re[mpNd2]+WX*im[mpNd2];
				re[mpNd2]=re[m]-TempX;
				im[mpNd2]=im[m]-TempY;
				re[m]=re[m]+TempX;
				im[m]=im[m]+TempY;
			}
		}
	}


}

void fft_N2_compl(unsigned int N, double * x, double * y, double * re, double * im) // N = 2^T
{

	double S;
	unsigned int I,J;
	unsigned char *Ic = (unsigned char*) &I;
    unsigned char *Jc = (unsigned char*) &J;

	unsigned int i,id2,k,m,mpNd2;
	long double WX,WY,TempX,TempY;

	for(i=0;i<N;i++)
	{
		re[i]=x[i];//memcpy
		im[i]=y[i];//memcpy
	}

	int T=-1;
	int _N=N;
	while(_N)
	{
		_N=_N>>1;
		T++;
	}
	for(I=1;I<N-1;I++)
	{
		Jc[0] = reverse256[Ic[3]];
	    Jc[1] = reverse256[Ic[2]];
	    Jc[2] = reverse256[Ic[1]];
	    Jc[3] = reverse256[Ic[0]];
	    J >>= (32 - T);
	    if (I < J)
	    {
	        S = re[I];
	        re[I] = re[J];
	        re[J] = S;
			S = im[I];
	        im[I] = im[J];
	        im[J] = S;
	    }
	}
	re[0]=x[0];
	im[0]=y[0];
	re[N-1]=x[N-1];
	im[N-1]=y[N-1];

	for(i=2,id2=1;i<=N;id2=i,i+=i)
	{
		for(k=0;k<id2;k++)
		{
			WX=cos(-PI_2*k/double(i));
			WY=sin(-PI_2*k/double(i));
			for(m=k;m<N;m+=i)
			{
				mpNd2=m+id2;
				TempX=WX*re[mpNd2]-WY*im[mpNd2];
				TempY=WY*re[mpNd2]+WX*im[mpNd2];
				re[mpNd2]=re[m]-TempX;
				im[mpNd2]=im[m]-TempY;
				re[m]=re[m]+TempX;
				im[m]=im[m]+TempY;
			}
		}
	}



}


void fft_N2_rev(unsigned int N, double * x, double * y, double * re, double * im) // N = 2^T
{

	double S;
	unsigned int I,J;
	unsigned char *Ic = (unsigned char*) &I;
    unsigned char *Jc = (unsigned char*) &J;

	unsigned int i,id2,k,m,mpNd2;
	long double WX,WY,TempX,TempY;

	for(i=0;i<N;i++)
	{
		re[i]=x[i];//memcpy
		im[i]=y[i];//memcpy
	}

	int T=-1;
	int _N=N;
	while(_N)
	{
		_N=_N>>1;
		T++;
	}
	for(I=1;I<N-1;I++)
	{
		Jc[0] = reverse256[Ic[3]];
	    Jc[1] = reverse256[Ic[2]];
	    Jc[2] = reverse256[Ic[1]];
	    Jc[3] = reverse256[Ic[0]];
	    J >>= (32 - T);
	    if (I < J)
	    {
	        S = re[I];
	        re[I] = re[J];
	        re[J] = S;
			S = im[I];
	        im[I] = im[J];
	        im[J] = S;
	    }
	}
	re[0]=x[0];
	im[0]=y[0];
	re[N-1]=x[N-1];
	im[N-1]=y[N-1];

	for(i=2,id2=1;i<=N;id2=i,i+=i)
	{
		for(k=0;k<id2;k++)
		{
			WX=cos(PI_2*k/double(i)); // +
			WY=sin(PI_2*k/double(i)); // +
			for(m=k;m<N;m+=i)
			{
				mpNd2=m+id2;
				TempX=WX*re[mpNd2]-WY*im[mpNd2];
				TempY=WY*re[mpNd2]+WX*im[mpNd2];
				re[mpNd2]=re[m]-TempX;
				im[mpNd2]=im[m]-TempY;
				re[m]=re[m]+TempX;
				im[m]=im[m]+TempY;
			}
		}
	}

	long double _1N=1.0/(long double)N;
	for(i=0;i<N;i++)
	{
		re[i]*=_1N;
		im[i]*=_1N;
	}



}



void fft_any_real(unsigned int N, double * x, double * re, double * im)
{
	unsigned int i;
	double *x2re,*x2im,*w_re,*w_im;
	long double PI_2_N=PI_2/(long double)N;
	long double TempX,TempY;
	int T=0;  //
	unsigned int N_=N; //~T=-1,_N=2N
	while(N_)
	{
		N_=N_>>1;
		T++;
	}
	N_=1<<(T+1);

	x2re=(double*)malloc(sizeof(double)*4*N_);
	x2im=x2re+N_;
	w_re=x2im+N_;
	w_im=w_re+N_;

	for(i=0;i<N;i++)
	{
		x2re[i]=x[i]*cos(PI_2_N*i*i*0.5);
		x2im[i]=x[i]*sin(PI_2_N*i*i*0.5);
	}
	memset(x2re+N,0,sizeof(double)*(N_-N));
	memset(x2im+N,0,sizeof(double)*(N_-N));

	int N22=2*N-2;
	for(i=0;i<N_;i++,N22--)
	{
		w_re[i]=cos(-PI_2_N*N22*N22*0.5);
		w_im[i]=sin(-PI_2_N*N22*N22*0.5);
	}

	fft_N2_compl(N_,x2re,x2im,x2re,x2im);
	fft_N2_compl(N_,w_re,w_im,w_re,w_im);
	for(i=0;i<N_;i++)
	{
		TempX=x2re[i]*w_re[i]-x2im[i]*w_im[i];
		TempY=x2re[i]*w_im[i]+x2im[i]*w_re[i];
		x2re[i]=TempX;
		x2im[i]=TempY;
	}
	fft_N2_rev(N_,x2re,x2im,x2re,x2im);

	N22=2*N-2;
	for(i=0;i<N;i++,N22--)
	{
		re[i]=x2re[N22]*cos(PI_2_N*i*i*0.5)-x2im[N22]*sin(PI_2_N*i*i*0.5);
		im[i]=x2re[N22]*sin(PI_2_N*i*i*0.5)+x2im[N22]*cos(PI_2_N*i*i*0.5);
	}

	free(x2re);
}


void fft(unsigned int N, double * x, double * re, double * im)
{
	int T=-1;
	int _N=N;
	while(_N)
	{
		_N=_N>>1;
		T++;
	}
	if(N==(1<<T))
		fft_N2_real(N,x,re,im);
	else
		fft_any_real(N,x,re,im);
}


void Ogib_N2(int N, double * x,double * out, double * tmp)
{
	int i;

	
	fft_N2_real(N,x,out,tmp);

	memset(out+N/2, 0, (N/2)*sizeof(double));
	memset(tmp+N/2, 0, (N/2)*sizeof(double));

	fft_N2_rev(N,out,tmp,out,tmp);
	
	for(i=0;i<N;i++)
		out[i]=sqrt( 4*tmp[i]*tmp[i] + x[i]*x[i] );


}







#endif