#ifndef ACORE_LIB_CPP
#define ACORE_LIB_CPP


#include <math.h>
//#include <fstream.h>
#include <stdio.h>
#include "ACore_lib.h"
#include "WAV_lib.h"
#include "FFT_lib2.cpp"
#include "GDI_lib.h"


double Relevant(double fr, double * ampl, int N, double step)
{
	// среднее значение амплитуды
	int i;

	if(fr/step>=N)
		return 0;
	double res=0;
	for(i=1;i<=5;i++)// (i=1; (i<=6) && int(i*fr/step)<N; i++)//13-14; 146-147; 
		res+=ampl[int(i*fr/step)];
//	res/=double(i-1); // т.к. i на последнем шаге вышло за пределы

	return res;
}

double MainACore(double * ampl, int N, double step)
{
	double core=0;
	double core_rel=0;
	double c,rel;

	int min_i=(int)(MIN_HZ/step);
	int max_i=(int)(MAX_HZ/step + 1);
	if(max_i>=N)max_i=N-1; // N/2? // и вообще, в Relevant() вроде проверяется
	int i,j;

	for(i=min_i;i<max_i;i++)
	{
		for(j=0;j<NUM_SECT;j++)
		{
			c = (i + section[j] + EPS_ADD)*step;
			rel=Relevant(c,ampl,N,step);
			if(rel>core_rel)
			{ // проверяем кратные пики
			//	int d = c/core +0.5;
			//	if( (d==1) || (    fabs((core*d-c)/step) > d     ) )
			//	{
					core=c;
					core_rel=rel;
			//	}
			}
		}
	}

	return core;
}

/////////// STATISTIC //////////////////////////////////////////////////////////////////////////////
/*
int Statistic(WavSound w, char * output, int gdi_tmp, double gdi_sc, bool gdi_draw=true)
{
//	int und=0; // undefined acoustic cores
	// пока что не учитывается возможность на каком-то участке сигнала не определить главное ядро

	double *w_st,*w_en;
	int i,n;
	int no_n;

	w_st=w.data;
	w_en=w.data+w.length;
	// fabs(*w_) <?.0 - дискуссионный вопрос
	while(fabs(*w_st)<5.0 &&
		w_st<w_en)
		w_st++;
	while(fabs(*w_en)<5.0 && 
		w_st<w_en)
		w_en--;
	n=(w_en-w_st)/FR_BUF_SIZE;

	if(n==0)
		return -1;

	double * Fr = (double*)malloc(n*sizeof(double)); // frequency
	double * A1 = (double*)malloc(n*sizeof(double)); // амплитуда первого пика
	double * E12 = (double*)malloc(n*sizeof(double)); // отношения пиков А1/А2
	double * E23 = (double*)malloc(n*sizeof(double));
	double * E34 = (double*)malloc(n*sizeof(double)); // и  вообще,
	double * E45 = (double*)malloc(n*sizeof(double)); // выделить сразу всем
	//EIJ - все в один 2д массив

	double FrA,A1A,E12A,E23A,E34A,E45A; // average
	double FrS,A1S,E12S,E23S,E34S,E45S; // sqrt(dispersion)

	
	double step=double(w.sampleRate)/double(FR_BUF_SIZE);

#define GDITMP gdi_tmp

	if(gdi_draw)GDI_Clear();

	//if(gdi_draw)GDI_Text(output,0,0);


	no_n=0;
	for(i=0;i<n;i++)
	{
		memcpy((void*)fr_buf, (void*)(w_st+i*FR_BUF_SIZE), FR_BUF_SIZE*sizeof(double));
		memset((void*)(fr_buf+FR_BUF_SIZE),0,FR_BUF_SIZE*sizeof(double));

		if(gdi_draw)GDI_DrawArray(FR_BUF_SIZE,fr_buf,RGB(0,150,0),1.0,START-GDITMP+150*i);
	
		fft_N2_real(FR_BUF_SIZE,fr_buf,fr_buf,fr_buf+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf);

		if(gdi_draw)GDI_DrawArray(FR_BUF_SIZE/2,fr_buf,RGB(150,0,0),gdi_sc,START-GDITMP+150*i);


		int numP=FindPeaks(fr_buf,FR_BUF_SIZE/2,peak_buf,PEAK_BUF_SIZE);
		double core = FindACore(peak_buf,numP,fr_buf,FR_BUF_SIZE/2); // fr_buf,FR_BUF_SIZE/2,step);

		
		if(core==0.0)
		{
			no_n++;
			continue;
		}

		Fr[i-no_n]=core*step;
		A1[i-no_n]=fr_buf[int(core)];
/*
		E12[i-no_n]=fr_buf[int(2*core)] / fr_buf[int(core)];
		E23[i-no_n]=fr_buf[int(3*core)] / fr_buf[int(core)];
		E34[i-no_n]=fr_buf[int(4*core)] / fr_buf[int(core)];
		E45[i-no_n]=fr_buf[int(5*core)] / fr_buf[int(core)];
/*//*
		E12[i-no_n]=fr_buf[int(core)] / fr_buf[int(2*core)];
		E23[i-no_n]=fr_buf[int(2*core)] / fr_buf[int(3*core)];
		E34[i-no_n]=fr_buf[int(3*core)] / fr_buf[int(4*core)];
		E45[i-no_n]=fr_buf[int(4*core)] / fr_buf[int(5*core)];
///*//*

		if(gdi_draw)GDI_DrawACore(core*step,step,FR_BUF_SIZE/2,fr_buf,RGB(100,100,0),gdi_sc,START-GDITMP+150*i);
	}
	n-=no_n;

	if(gdi_draw)GDI_DrawArray(n,Fr,RGB(150,150,150),1.0,START+250);

	// собственно статистика
	FrA=0;A1A=0;E12A=0;E23A=0;E34A=0;E45A=0;
	for(i=0;i<n;i++)
	{
		FrA+=Fr[i];
		A1A+=A1[i];
		E12A+=E12[i];
		E23A+=E23[i];
		E34A+=E34[i];
		E45A+=E45[i];
	}
	FrA/=double(n);
	A1A/=double(n);
	E12A/=double(n);
	E23A/=double(n);
	E34A/=double(n);
	E45A/=double(n);

	FrS=0;A1S=0;E12S=0;E23S=0;E34S=0;E45S=0;
	for(i=0;i<n;i++)
	{
		FrS+=(Fr[i]-FrA)*(Fr[i]-FrA);
		A1S+=(A1[i]-A1A)*(A1[i]-A1A);
		E12S+=(E12[i]-E12A)*(E12[i]-E12A);
		E23S+=(E23[i]-E23A)*(E23[i]-E23A);
		E34S+=(E34[i]-E34A)*(E34[i]-E34A);
		E45S+=(E45[i]-E45A)*(E45[i]-E45A);
	}
	FrS/=double(n);
	FrS=sqrt(FrS);
	A1S/=double(n);
	A1S=sqrt(A1S);
	E12S/=double(n);
	E12S=sqrt(E12S);
	E23S/=double(n);
	E23S=sqrt(E23S);
	E34S/=double(n);
	E34S=sqrt(E34S);
	E45S/=double(n);
	E45S=sqrt(E45S);


	// output:
	/*
	ofstream os(output);
	os << "Number of segments: " << n << "\n";
	os << "Frequency : " << FrA << "		" << FrS << "\n";
	os << "Energy 1  : " << A1A << "		" << A1S << "\n";
	os << "Energy 1/2: " << E12A << "		" << E12S << "\n";
	os << "Energy 2/3: " << E23A << "		" << E23S << "\n";
	os << "Energy 3/4: " << E34A << "		" << E34S << "\n";
	os << "Energy 4/5: " << E45A << "		" << E45S << "\n";

	for(i=0;i<n;i++) //
		os << Fr[i] << "\n";
	os.close();
	//*//*


	free(Fr);
	free(A1);
	free(E12);
	free(E23);
	free(E34);
	free(E45);

	return 0;
}//*/



///////// FIND /////////////////////////////////////////////////////////////////////////////////
#define COEF_PEAK (2.0)
int FindPeaks(double *frArr, int N, int *peakArr, int maxP)
{
	int i;
	int numP;

	int state=0; // 0-поиск пика, 1-поиск минимума
	double min,max;
	double fr;
	
	min=max=frArr[1];
	numP=0;
	i=2;

	while(i<N && numP<maxP)
	{
		
		while(i<N) // ищем пик
		{
			fr=frArr[i];
			if(fr<min)
			min=fr;
			if(numP>0)
			{
				if(fr>COEF_PEAK*min && fr>min+15) // для 4096 было 100 вроде // нашли пик
				{
					max=fr;
					peakArr[numP]=i;
					i++;
					numP++;
					break;
				}
			}
			else
				if(fr>COEF_PEAK*min && fr>300) // для 4096 было 2000 и 2*коэф // дубль, позже убрать
				{
					max=fr;
					peakArr[numP]=i;
					i++;
					numP++;
					break;
				}

			i++;
		}

		while(i<N) // ищем максимум пика + спад
		{
			fr=frArr[i];
			if(fr>max)
			{
				max=fr;
				peakArr[numP-1]=i;
			}
			if(fr<max/COEF_PEAK)
			{
				min=fr;
				i++;
				break;
			}
			i++;
		}

	}

	return numP;
}

/*
// OLD
void FindFormants(int N, double * x, int maxP, int *peakArr)// more eps, less osred1	// eps=0.001 , osred1-300
{
	FILE * file_out=fopen("out.txt","w");
	double eps=100;
	double Eps=100;

	double max;
	int mi;
	int numP=0;
	for(int i=1;i<N-1 && numP<maxP;i++)
	{
		if(x[i]<x[i-1]+eps)
			continue;
		if(x[i]<x[i+1]-eps)
			continue;
		if(x[i]>x[i+1]+eps)
		{
			GDI_DrawFormant(i,x[i]);
			fprintf(file_out,"\nformant: %i",i);
			peakArr[numP]=i;
			numP++;
			continue;
		}
//*		
		max=x[i];
		mi=i;
		do
		{
			i++;
			if(x[i]>max){max=x[i];mi=i;}
		}
		while( (fabs(max-x[i+1])<Eps) && i<N-1 && numP<maxP );

		if(i==N-1)
		{
			GDI_DrawFormant(mi,max);
			fprintf(file_out,"\nformant: %i",mi);
			peakArr[numP]=mi;
			numP++;
			break;
		}
		if(max<x[i+1])
			continue;

		GDI_DrawFormant(mi,max);
		fprintf(file_out,"\nformant: %i",mi);//*
		peakArr[numP]=mi;
		numP++;
	}	

	fclose(file_out);
	return;
}
//*/


inline double Max(double a, double b)
{
	return  (a>b)?a:b;
}

inline double Min(double a, double b)
{
	return  (a>b)?b:a;
}


double FindACore(int *peakArr, int maxP, double *frArr, int frN)
{
	if(maxP==0)
		return 0;

	double	fr_st=peakArr[0]-0.5,
			fr_en=peakArr[0]+0.5;


	int i_st,i_en;
	int i,
		n1,n2;

	n1=1;
	for(i=2;i*peakArr[0]<frN && n1<maxP;i++)
	{
		i_st = (int)(i*fr_st+0.5);
		i_en = (int)(i*fr_en+0.5);

		while(n1<maxP && peakArr[n1]<i_st)n1++;
		n2=n1;
		while(n2<maxP && peakArr[n2]<=i_en)n2++;
		if(n1>=maxP)
			break;
		if(peakArr[n1]>i_en)
			continue;

		int max;
		for(max=n1,n1++;n1<n2;n1++)
		{
			if(frArr[peakArr[n1]]>=frArr[peakArr[max]])
				max=n1;
		}

		fr_st=Max( double(peakArr[max]-0.5)/double(i) , fr_st );
		fr_en=Min( double(peakArr[max]+0.5)/double(i) , fr_en );

	}

	return (fr_st+fr_en)/2.0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

double rand10()
{
	return double(rand())/double(RAND_MAX);
}

double ACoreRandomGenerate(double* gen, int N)
{
	int i,j;
	double fr = N*(0.005 + 0.005*rand10());
	for(i=0;i<N;i++)
		gen[i]=0;
	int nf=(int)(N/(2*fr)-1);
	if(nf>30)nf=30;
//	nf=2;
	for(j=1;j<nf;j++)
	{
		double re=N*(2.0+rand10())/log(double(j+1));//N*(2.5)/log(double(j+1));
		double im=N*(1.0+rand10())/log(double(j+1));
		for(i=0;i<N;i++)
			gen[i]+=re*cos(PI_2*i*j*fr/double(N))+im*sin(PI_2*i*j*fr/double(N));
	}
	for(i=0;i<N;i++)
		gen[i]/=double(N)/2.0;

	return fr;
}


inline void compl_div(double z1re, double z1im, double z2re, double z2im, double *re, double *im)
{
	double z2mod2=z2re*z2re+z2im*z2im;

	*re=(z1re*z2re+z1im*z2im)/z2mod2;
	*im=(z1im*z2re-z2im*z1re)/z2mod2;

}


void ACoreGenerate(double freq, double * gen, double *re, double *im, int N)
// freq - выражена в отсчетах, не в герцах
{
	if(freq<=1.0)
	{
		memset(gen,0,N*sizeof(double));
		return;
	}


	int i,j;
	int nf= (int)(double(N/2-0.5)/freq);
	double *coefs = (double*)malloc(2*nf*sizeof(double)); // coefs[2*i]=re[i] coefs[2*i+1]=im[i]
	double *buf   = (double*)malloc(2*N*sizeof(double));

	for(i=0;i<N;i++)
	{
		gen[i]=0;
		for(j=1;j<nf;j++)
			gen[i]+=   re[int(freq*j+0.5)]*cos(PI_2*i*j*freq/double(N))-im[int(freq*j+0.5)]*sin(PI_2*i*j*freq/double(N));
		gen[i] *= 2.0/double(N);
	}

	fft_N2_real(FR_BUF_SIZE,gen,buf,buf+N);

	for(i=1;i<=nf;i++)
	{
		int ind = int(i*freq+0.5);
		compl_div(re[ind],im[ind],buf[ind],buf[N+ind],coefs+2*i-2,coefs+2*i-1);
	}


	double c_r,c_i;
	for(i=0;i<N;i++)
		gen[i]=0;
	for(j=1;j<nf;j++)
	{
		c_r=coefs[2*j-2]*re[int(freq*j+0.5)] - coefs[2*j-1]*im[int(freq*j+0.5)];
		c_i=coefs[2*j-1]*re[int(freq*j+0.5)] + coefs[2*j-2]*im[int(freq*j+0.5)];
		for(i=0;i<N;i++)
		{
			gen[i]+=   c_r*cos(PI_2*i*j*freq/double(N))-c_i*sin(PI_2*i*j*freq/double(N));
		}
	}
	for(i=0;i<N;i++)
		gen[i] /= double(N)/2.0;

/*	//  выравнивание краев сигнала - пробуем
	double diff=gen[N-1]-gen[0];
	for(i=0;i<N;i++)
		gen[i]-=diff*i/double(N-1);//*/


	free(buf);
	free(coefs);
	return;
}

void ACoreGenerateEX(double freq, double * gen, double *ar, int N)
// freq - выражена в отсчетах, не в герцах
{
	if(freq<=1.0)
	{
		memset(gen,0,N*sizeof(double));
		return;
	}

	

	int i,j;
	int nf= (int)(double(N/2-0.5)/freq);

	int nn= (int)(double(N)/freq * int(freq));// предельный отсчет для разложения
	
	double *coefs = (double*)malloc(2*nf*sizeof(double)); // coefs[2*i]=re[i] coefs[2*i+1]=im[i]
	double *buf   = (double*)malloc(2*N*sizeof(double));

	for(j=0;j<nf;j++)
	{
		coefs[2*j]=coefs[2*j+1]=0;
		for(i=0;i<nn;i++) // i<nn,N
		{
			coefs[2*j]+=cos(PI_2*(j*freq)*i/double(N)) * ar[i];
			coefs[2*j+1]+=sin(PI_2*(j*freq)*i/double(N)) * ar[i];
		}
	}

	for(i=0;i<N;i++)
	{
		gen[i]=0;
		for(j=1;j<nf;j++)
			//gen[i]+= sqrt( coefs[2*j]*coefs[2*j]+coefs[2*j+1]*coefs[2*j+1] ) * sin(PI_2*i*j*freq/double(N));
			// - половина метода синус-модульной стыковки
			gen[i]+=   coefs[2*j]*cos(PI_2*i*j*freq/double(N))+coefs[2*j+1]*sin(PI_2*i*j*freq/double(N));
		gen[i] /= double(nn)/2.0; // nn/2,N/2
	}



	free(buf);
	free(coefs);
	return;
}



double fr_buf2[2*FR_BUF_SIZE];
double gener[FR_BUF_SIZE],gener2[FR_BUF_SIZE];
double gen_fr[2*FR_BUF_SIZE];
/*
int ACoreSynthes(WavSound w, char * output, int gdi_tmp=0, double gdi_sc=SCALE)
{
//	int und=0; // undefined acoustic cores
	// пока что не учитывается возможность на каком-то участке сигнала не определить главное ядро

	double *w_st,*w_en;
	int i,n;
	int no_n;

	w_st=w.data;
	w_en=w.data+w.length;
	// fabs(*w_) <?.0 - дискуссионный вопрос
	while(fabs(*w_st)<5.0 &&
		w_st<w_en)
		w_st++;
	while(fabs(*w_en)<5.0 && 
		w_st<w_en)
		w_en--;
	n=(w_en-w_st)/FR_BUF_SIZE;

	if(n==0)
		return -1;

	
	double step=double(w.sampleRate)/double(FR_BUF_SIZE);
	double *gener=(double*)malloc(FR_BUF_SIZE*sizeof(double));



	GDI_Clear();
	//GDI_Text(output,0,0);

	no_n=0;
	i=(n>25)?25:5;
	memcpy((void*)gener, (void*)(w_st+i*FR_BUF_SIZE), FR_BUF_SIZE*sizeof(double));
	//double precise_fr=ACoreRandomGenerate(gener,FR_BUF_SIZE);


	// PROVERKA PEREHODA PIKA --- FAIL ///////////////////////////
	
	//	for(i=0;i<FR_BUF_SIZE;i++)
	//		gener[i]=10*cos(11.871+PI_2*i*16.49245564997638/2048.0);	
	//	GDI_DrawArray(FR_BUF_SIZE,gener,RGB(0,150,0),1.0,START+100);

		memcpy((void*)fr_buf, (void*)gener, FR_BUF_SIZE*sizeof(double));
		memset((void*)(fr_buf+FR_BUF_SIZE),0,FR_BUF_SIZE*sizeof(double));

	/*
		fft_N2_real(FR_BUF_SIZE,fr_buf,fr_buf,fr_buf+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf2);
		GDI_DrawArrayEX(FR_BUF_SIZE,fr_buf2,RGB(150,0,0),0.01,START+100);
		
		double max_d=0;
		int max_i=0;
		for(i=0;i<FR_BUF_SIZE/2;i++)
		{
			if(fr_buf2[i]>max_d)
			{
				max_d=fr_buf2[i];
				max_i=i;
			}
		}

		ampl_count(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf2);
	//*////////////////////////////////////PROVERKA PEREHODA
		

/*
		///////////////
		// GO CEPSTRUM
		int j;
		for(j=0;j<FR_BUF_SIZE;j++)
			fr_buf2[j]=log(fr_buf2[j]*fr_buf2[j]);
		GDI_DrawArray(FR_BUF_SIZE/2,fr_buf2,RGB(50,50,0),1.0,START+200);
		fft_N2_real(FR_BUF_SIZE,fr_buf2,fr_buf2,fr_buf2+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,fr_buf2,fr_buf2+FR_BUF_SIZE,fr_buf2);
		GDI_DrawArray(FR_BUF_SIZE/2,fr_buf2,RGB(150,150,0),gdi_sc,START+200);
		// END CEPSTRUM
		///////////////
//*//*		


////////// PROVERKA FORMUL PERESCHETA ////////////
//*		
		double w0= 271.4;
		double ww = 44100.0/double(FR_BUF_SIZE);
		for(i=1;i<FR_BUF_SIZE;i++)
		{
			fr_buf[i]=0;
			fr_buf[i+FR_BUF_SIZE]=0;
			for(int j=0;j<50;j++)
			{
				fr_buf[i]				+=	300*w0*sin(2*PI*w0/ww) / (PI/ww*(w0*w0-double(i+j*FR_BUF_SIZE)*double(i+j*FR_BUF_SIZE)*ww*ww));
				fr_buf[i+FR_BUF_SIZE]	+=	300*double(i+j*FR_BUF_SIZE)*ww*(cos(2*PI*w0/ww)-1.0) / (PI/ww*(w0*w0-double(i+j*FR_BUF_SIZE)*double(i+j*FR_BUF_SIZE)*ww*ww));
			}
		}
		fr_buf[0]=fr_buf[FR_BUF_SIZE]=fr_buf[FR_BUF_SIZE-1]=fr_buf[2*FR_BUF_SIZE-1]=0;
		GDI_DrawArray(FR_BUF_SIZE,fr_buf,RGB(0,150,0),1.0,START+100);
		///////////////////////////
		// formula!
		for(i=0;i<FR_BUF_SIZE;i++)
		{
			fr_buf[i]=1.0*cos(2.0*PI*w0*i/44100.0);
		}
		fft_N2_real(FR_BUF_SIZE,fr_buf,fr_buf,fr_buf+FR_BUF_SIZE);
		int _max=0;double max=0;
		for(i=1;i<FR_BUF_SIZE && fabs(fr_buf[i])>max;i++)
		{
			max=fr_buf[i];
			_max=i;
		} // _max=13!!!
		double w_0=(fr_buf[_max]*_max*_max - fr_buf[_max-1]*(_max-1.0)*(_max-1.0))/(fr_buf[_max]-fr_buf[_max-1]);
		w_0=sqrt(w_0);
		w_0*=ww;
		// rabotaet! tolko problema s nasloeniem iz-za DFT -> schitat nuzhno na pike, togda minimum pogreshnosti
		///////////////////////////*//*

	//	fr_buf[0]=fr_buf[FR_BUF_SIZE-1]=fr_buf[FR_BUF_SIZE]=fr_buf[2*FR_BUF_SIZE-1]=0;
		fft_N2_rev(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE);
		GDI_DrawArray(FR_BUF_SIZE,fr_buf,RGB(150,0,0),100.0,START+100);
		ampl_count(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf2);

		return 0;











		int numP=FindPeaks(fr_buf2,FR_BUF_SIZE/2,peak_buf,PEAK_BUF_SIZE);

		double core = FindACore(peak_buf,numP,fr_buf2,FR_BUF_SIZE/2); // fr_buf,FR_BUF_SIZE/2,step);
		GDI_DrawACore(core*step,step,FR_BUF_SIZE/2,fr_buf2,RGB(100,100,0),gdi_sc,START+100);

		
		if(core==0.0)
		{
			return 0;
		}

		//ACoreGenerate(core,gener2,fr_buf,fr_buf+FR_BUF_SIZE,FR_BUF_SIZE);
		ACoreGenerateEX(core,gener2,gener,FR_BUF_SIZE);
		//				core - precise_fr

		/* синус-модульное преобразование
		fft_N2_real(FR_BUF_SIZE,gener2,fr_buf,fr_buf+FR_BUF_SIZE);
		fr_buf[0]=fr_buf[FR_BUF_SIZE]=0.0;
		for(int j=1;j<=FR_BUF_SIZE/2;j++)
		{
			fr_buf[FR_BUF_SIZE+j] = sqrt(fr_buf[j]*fr_buf[j] + fr_buf[j+FR_BUF_SIZE]*fr_buf[j+FR_BUF_SIZE]) ;
			//		im[j]		=	sqrt	(re[j]*re[j]+im[j]*im[j])
			fr_buf[2*FR_BUF_SIZE-j] =  -fr_buf[FR_BUF_SIZE+j];

			fr_buf[j]=0.0;
			fr_buf[FR_BUF_SIZE-j]=0.0;
		}
		fft_N2_rev(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,gener2,fr_buf+FR_BUF_SIZE);
		//*//*



		GDI_DrawArray(FR_BUF_SIZE,gener2,RGB(0,150,0),1.0,START+300);

		fft_N2_real(FR_BUF_SIZE,gener2,gen_fr,gen_fr+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,gen_fr,gen_fr+FR_BUF_SIZE,fr_buf2);
		GDI_DrawArray(FR_BUF_SIZE,fr_buf2,RGB(150,0,0),gdi_sc,START+300); // fbs/2!
		GDI_DrawACore(core*step,step,FR_BUF_SIZE/2,fr_buf2,RGB(100,100,0),gdi_sc,START+300);
//		// работаем с остатком сигнала
		SaveWAV("gener.wav",220500,44100,FR_BUF_SIZE,gener2);
		int k;
		for(k=0;k<FR_BUF_SIZE;k++)
			gener[k]=gener[k]-gener2[k];
		GDI_DrawArray(FR_BUF_SIZE,gener,RGB(0,150,0),1.0,START+500);
		fft_N2_real(FR_BUF_SIZE,gener,gen_fr,gen_fr+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,gen_fr,gen_fr+FR_BUF_SIZE,fr_buf2);
		GDI_DrawArray(FR_BUF_SIZE/2,fr_buf2,RGB(150,0,0),gdi_sc,START+500);
		GDI_DrawACore(core*step,step,FR_BUF_SIZE/2,fr_buf2,RGB(100,100,0),gdi_sc,START+500);
		SaveWAV("gener_ost.wav",220500,44100,FR_BUF_SIZE,gener);

		

	free(gener);
	return 0;
}
//*/

void ACoreFileMake(WavSound w)
{
	double *w_st,*w_en;
	int i,n;
//	int no_n;

	w_st=w.data;
	w_en=w.data+w.length;
	n=(w_en-w_st)/FR_BUF_SIZE;

	if(n==0)
		return;

	
	double step=double(w.sampleRate)/double(FR_BUF_SIZE);
	double *gener=(double*)malloc(n*FR_BUF_SIZE*sizeof(double));
	double *gener_ost=(double*)malloc(n*FR_BUF_SIZE*sizeof(double));

	for(i=0;i<n;i++)
	{
		memcpy((void*)fr_buf, (void*)(w_st+i*FR_BUF_SIZE), FR_BUF_SIZE*sizeof(double));
		memset((void*)(fr_buf+FR_BUF_SIZE),0,FR_BUF_SIZE*sizeof(double));
	
		fft_N2_real(FR_BUF_SIZE,fr_buf,fr_buf,fr_buf+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf2);

		int numP=FindPeaks(fr_buf2,FR_BUF_SIZE/2,peak_buf,PEAK_BUF_SIZE);

		double core = FindACore(peak_buf,numP,fr_buf2,FR_BUF_SIZE/2);
		ACoreGenerateEX(core,gener+i*FR_BUF_SIZE,w_st+i*FR_BUF_SIZE,FR_BUF_SIZE);

		//  стыковка синус-модулем - получается плохо
		/*
		fft_N2_real(FR_BUF_SIZE,gener+i*FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE);
		fr_buf[0]=fr_buf[FR_BUF_SIZE]=0.0;
		for(int j=1;j<=FR_BUF_SIZE/2;j++)
		{
			fr_buf[j+FR_BUF_SIZE] = sqrt(fr_buf[j]*fr_buf[j] + fr_buf[j+FR_BUF_SIZE]*fr_buf[j+FR_BUF_SIZE]) ;
			//		im[j]		=	sqrt	(re[j]*re[j]+im[j]*im[j])
			fr_buf[2*FR_BUF_SIZE-j] = -fr_buf[j+FR_BUF_SIZE];

			fr_buf[j]=0.0;
			fr_buf[FR_BUF_SIZE-j]=0.0;
		}
		fft_N2_rev(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,gener+i*FR_BUF_SIZE,fr_buf+FR_BUF_SIZE);//*/




		int k;
		for(k=0;k<FR_BUF_SIZE;k++)
			gener_ost[i*FR_BUF_SIZE+k]=w_st[i*FR_BUF_SIZE+k]-gener[i*FR_BUF_SIZE+k];
	}

	// стыковка смещением - несколько лучше
	//*
	double smesh=0;
	double *sme=(double*)malloc(n*sizeof(double));

	for(i=0;i<n-1;i++)
		sme[i]= (gener[(i+1)*FR_BUF_SIZE] - gener[(i+1)*FR_BUF_SIZE-1]) - (w_st[(i+1)*FR_BUF_SIZE] - w_st[(i+1)*FR_BUF_SIZE-1]);
	for(i=0;i<n-1;i++)
	{
		smesh+=sme[i];
		for(int j=0;j<FR_BUF_SIZE;j++)
			gener[(i+1)*FR_BUF_SIZE+j]-=smesh;
	}
	for(i=0;i<n*FR_BUF_SIZE;i++)
	{
		gener[i]+=smesh/2.0;
	//	gener_ost[i]=w_st[i]-gener[i]; - пока будем считать gener_ost по несмещенному gener
	}

	for(i=0;i<n-1;i++)
		sme[i]= (gener_ost[(i+1)*FR_BUF_SIZE] - gener_ost[(i+1)*FR_BUF_SIZE-1]);
	for(i=0;i<n-1;i++)
	{
		smesh+=sme[i];
		for(int j=0;j<FR_BUF_SIZE;j++)
			gener_ost[(i+1)*FR_BUF_SIZE+j]-=smesh;
	}
	for(i=0;i<n*FR_BUF_SIZE;i++)
	{
		gener_ost[i]+=smesh/2.0;
	}

	free(sme);
	//*/

	SaveWAV("_origin.wav",w.length,w.sampleRate,w.length,w.data);
	SaveWAV("_gener.wav",n*FR_BUF_SIZE,w.sampleRate,n*FR_BUF_SIZE,gener);
	SaveWAV("_gener_ost.wav",n*FR_BUF_SIZE,w.sampleRate,n*FR_BUF_SIZE,gener_ost);


	free(gener_ost);
	free(gener);

}



/*// проверка диапазона отсчета для синусоиды вещественного периода
	double f;
	int frMax;
	FILE *fl=fopen("out_tmp.txt","w");
	for(f=170.0*step;f<175*step;f+=0.02*step)
	{
		for(int ii=0;ii<FR_BUF_SIZE;ii++)
			fr_buf[ii]=25*sin(ii*2*PI*f/double(w.sampleRate));
		memset((void*)(fr_buf+FR_BUF_SIZE),0,FR_BUF_SIZE*sizeof(double));

		GDI_DrawArray(FR_BUF_SIZE,fr_buf,RGB(0,150,0),1.0,START);
		fft_N2_real(FR_BUF_SIZE,fr_buf,fr_buf,fr_buf+FR_BUF_SIZE);
		ampl_count(FR_BUF_SIZE,fr_buf,fr_buf+FR_BUF_SIZE,fr_buf);

		GDI_DrawArray(FR_BUF_SIZE/2,fr_buf,RGB(150,0,0),0.01,START);

		frMax=0;
		for(ii=1;ii<FR_BUF_SIZE/2;ii++)
			if(fr_buf[frMax]<fr_buf[ii])
				frMax=ii;

		fprintf(fl,"%i %f\n",frMax,3*f/step);
	}
	fclose(fl);//*/


#endif