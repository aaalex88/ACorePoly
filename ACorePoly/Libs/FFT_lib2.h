#ifndef FFT_LIB2_H
#define FFT_LIB2_H


#define PI   (3.1415926535897932384626433832795)
#define PI_2 (6.283185307179596476925286766559)


void ampl_count(unsigned int N, double * re, double * im, double * x);
void fft_N2_real(unsigned int N, double * x, double * re, double * im);
void fft_N2_compl(unsigned int N, double * x, double * y, double * re, double * im);
void fft_N2_rev(unsigned int N, double * x, double * y, double * re, double * im);
void fft_any_real(unsigned int N, double * x, double * re, double * im);

void fft(unsigned int N, double * x, double * re, double * im);

void Ogib_N2(int N, double * x,double * out, double * tmp);




#endif