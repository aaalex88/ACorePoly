#ifndef ACORE_LIB_H
#define ACORE_LIB_H


#include "WAV_lib.h"



#define FR_BUF_SIZE 2048 // 4096
#define PEAK_BUF_SIZE 1024


#define MIN_HZ 70.0
#define MAX_HZ 700.0

#define EPS_ADD 0.001

#define NUM_SECT 12
const double section[NUM_SECT]={0.0/6.0 , 1.0/6.0 , 1.2/6.0 , 1.5/6.0 , 2.0/6.0 , 2.4/6.0,
3.0/6.0 , 3.6/6.0 , 4.0/6.0 , 4.5/6.0 , 4.8/6.0 , 5.0/6.0 };

double Relevant(double fr, double * ampl, int N, double step);
// fr - проверяемая частота, ampl - массив отсчетов частот, N - кол-во отсчетов, step - шаг частоты 
// между отсчетами

double MainACore(double * ampl, int N, double step);



double fr_buf[2*FR_BUF_SIZE]; // 0..FR_BUF_SIZE-1 - re; FR_BUF_SIZE..2*FR_BUF_SIZE-1 - im;


int Statistic(WavSound w, char * output, int gdi_tmp,  double gdi_sc, bool gdi_draw);
int Statistic(WavSound w, double **array, int *N);


int peak_buf[PEAK_BUF_SIZE];
int FindPeaks(double *frArr, int N, int *peakArr, int maxP);	// один хрен, просто разными способами
void FindFormants(int N, double * x, int maxP, int *peakArr);	//


double Max(double a, double b);
double Min(double a, double b);

double FindACore(int *peakArr, int maxP, double *frArr, int frN);


void compl_div(double z1re, double z1im, double z2re, double z2im, double *re, double *im);
double ACoreRandomGenerate(double* gen, int N);
void ACoreGenerate(double freq, double * gen, double *re, double *im, int N);
void ACoreGenerateEX(double freq, double * gen, double *ar, int N);
int ACoreSynthes(WavSound w, char * output, int gdi_tmp, double gdi_sc);


void ACoreFileMake(WavSound w);


#endif