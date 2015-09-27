#ifndef ACORE_LIB2_CPP
#define ACORE_LIB2_CPP


#include <math.h>
#include <stdlib.h>

#include "ACore_lib2.h"
#include "ACore_lib3.h"
#include "FFT_lib2.h"




#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))


void ACSave(FILE * f, ACore * c)
{
	fprintf(f,"time %f\n",c->t0);
	fprintf(f,"frequency %f\n",c->baseFr);
	fprintf(f,"number of amplitudes %i\n",c->numAmpl);
	for(int i=0;i<c->numAmpl;i++)
	{
		fprintf(f,"ampl %i		%f		%f\n",i+1,c->A[i],c->B[i]);
	}
	fprintf(f,"\n");
}

void ACSave(const char * fname, ACore * c)
{
	FILE * f;
	fopen_s(&f,fname,"w");

	ACSave(f,c);
	
	fclose(f);
}

void ACSave(const char * fname, ACoresSignal * cs)
{
	FILE * f;
	fopen_s(&f,fname,"w");
	fprintf(f,"length %f\n",cs->length);
	fprintf(f,"number of cores %i\n\n",cs->numCores);

	for(int i=0;i<cs->numCores;i++)
	{
		ACNorm(cs->cores+i);
		ACSave(f,cs->cores+i);
	}

	fclose(f);
}


void ACCopy(ACore * c, ACore * c1)
{
	c->baseFr=c1->baseFr;
	c->numAmpl=c1->numAmpl;
	c->t0=c1->t0;
	c->t1=c1->t1;
	c->t2=c1->t2;

	for(int i=0;i<ACORE_MAX_ELEM;i++)
	{
		c->A[i]=c1->A[i];
		c->B[i]=c1->B[i];
	}
}


void ACMix(ACore * c, ACore * c1, ACore * c2, double d1, double d2)
{
	c->baseFr=d1*c1->baseFr+d2*c2->baseFr;
	c->numAmpl=MAX(c1->numAmpl,c2->numAmpl);
	c->t0=d1*c1->t0+d2*c2->t0; // оставляем, т.к. по сути так оно и есть
	//c->t1=d1*c1->t1+d2*c2->t1; // не нужно, тк теряется смысл понятия отрезка
	//c->t2=d1*c1->t2+d2*c2->t2;

	for(int i=0;i<ACORE_MAX_ELEM;i++)
	{
		c->A[i]=d1*c1->A[i]+d2*c2->A[i];
		c->B[i]=d1*c1->B[i]+d2*c2->B[i];
	}
}


void	ACNorm(ACore * c)
{
	for(int i=0;i<ACORE_MAX_ELEM;i++)
	{
		c->A[i]=sqrt(c->A[i]*c->A[i]+c->B[i]*c->B[i]);
		c->B[i]=0;
	}
}


double Angle(double dCos, double dSin)
{
	double dd= sqrt(dCos*dCos+dSin*dSin);
	if(dd==0)
		return 0;
	dCos/=dd;
	dSin/=dd;
	double an=acos(dCos);
	if(dSin<0)
		an=PI_2-an;

	return an;
}



// добавить вариант реализации для ядра с коэф-тами вида амплитуда/фаза
void GenerateSignalStatic(int N, double * signal, double dt, ACore * core)
{
	int i,k;
	long double ld;

	for(i=0;i<N;i++)
	{
		ld=0;
		//signal[i]=0;
		for(k=0;k<core->numAmpl;k++)
		{
			ld+=core->A[k]*cos(PI_2*double(k+1)*dt*i*core->baseFr)+core->B[k]*sin(PI_2*double(k+1)*dt*i*core->baseFr);
			//signal[i]+=core->A[k]*cos(PI_2*double(k+1)*dt*i*core->baseFr)+core->B[k]*sin(PI_2*double(k+1)*dt*i*core->baseFr);
		}
		signal[i]=ld;
	}
}



// Пока что линейно, надо добавить сплайны
// а интегрирование фазы - ? попробуем трапецией?
// Соответствие фаз для разных ядер не считаем - оставляем только амплитуды

// считаем, что фазы в ядрах считаются отн-но точки t1 те в первом ядре фазы для ноля
void GenerateSignalDynamic(int N, double * signal, double dt, ACoresSignal * cs)
{
	int i,j,k;
	int end;
	ACore tempcor,tc1,tc2;

	double phase=0;
	double phShift[ACORE_MAX_ELEM]; // пока не используется

	double l=N*dt;
	

	if(cs->numCores==0)
		return;

	// начало - core[0]
	ACCopy(&tempcor,cs->cores); // *(cs->cores[0])
	for(i=0;i<tempcor.numAmpl;i++)
		phShift[i]=Angle(tempcor.A[i],tempcor.B[i]);
	ACNorm(&tempcor);

	end= (int)(cs->cores[0].t0/dt);//+0.5?
	for(i=0;i<end && i<N;i++)
	{
		signal[i]=0;
		for(k=0;k<tempcor.numAmpl;k++)
			signal[i]+=tempcor.A[k]*cos(PI_2*double(k+1)*phase + phShift[k]);
		phase+=dt*tempcor.baseFr; // интегрирование фазы - прямоуг.
		// но пока частота - конст. так что норм.
	}
	

	// цикл с промежуточными ядрами
	// i уже установлено, - ...
	double oldBaseFr=tempcor.baseFr;
	for(j=1;j<cs->numCores;j++)
	{
		 ACCopy(&tc1,cs->cores+j-1);
		 ACCopy(&tc2,cs->cores+j);
		 ACNorm(&tc1);
		 ACNorm(&tc2);

		 end= (int)(tc2.t0/dt);
		 for(;i<end && i<N;i++)
		 {
			ACMix(&tempcor,&tc1,&tc2,  (tc2.t0-i*dt)/(tc2.t0-tc1.t0)  ,  (i*dt-tc1.t0)/(tc2.t0-tc1.t0) );
			signal[i]=0;
			for(k=0;k<tempcor.numAmpl;k++)
				signal[i]+=tempcor.A[k]*cos(PI_2*double(k+1)*phase + phShift[k]);
			phase+=dt*(tempcor.baseFr+oldBaseFr)*0.5;
			oldBaseFr=tempcor.baseFr;
		 }
	}

	// заканчиваем, последнее core
	ACCopy(&tempcor,cs->cores+(cs->numCores-1));
	for(;i<N;i++)
	{
		signal[i]=0;
		for(k=0;k<tempcor.numAmpl;k++)
			signal[i]+=tempcor.A[k]*cos(PI_2*double(k+1)*phase + phShift[k]);
		phase+=dt*tempcor.baseFr;
	}


	// while noise
	//for(i=0;i<N;i++)
	//	signal[i]+=2*(rand()/((double)RAND_MAX) - 0.5);
	// while noise
	//gauss noise
	//for(i=0;i<N;i++)
	//	signal[i]+=cos(PI_2*rand()/((double)RAND_MAX))*sqrt(-2.0*log((rand()/((double)RAND_MAX))));
	//gauss noise

}

void GenerateSignalDynamicNO(int N, double * signal, double dt, ACoresSignal * cs)
{
	int i,j,k;
	int end;
	ACore tc;

	double phase;
//	double phShift[ACORE_MAX_ELEM]; // пока не используется

	double l=N*dt;
	

	if(cs->numCores==0)
		return;


	i=0;
	for(j=0;j<cs->numCores;j++)
	{
		 ACCopy(&tc,cs->cores+j);
		 end= (int)(tc.t2/dt);
		 phase=0;
		 for(;i<end && i<N;i++)
		 {
			signal[i]=0;
			for(k=0;k<tc.numAmpl;k++)
				signal[i]+=tc.A[k]*cos(PI_2*double(k+1)*phase)+tc.B[k]*sin(PI_2*double(k+1)*phase);
			phase+=dt*tc.baseFr;
		 }
	}

	for(;i<N;i++)
	{
		signal[i]=0;
		for(k=0;k<tc.numAmpl;k++)
			signal[i]+=tc.A[k]*cos(PI_2*double(k+1)*phase)+tc.B[k]*sin(PI_2*double(k+1)*phase);
		phase+=dt*tc.baseFr;

	}
	 

}




void RandomSignalStatic(int N, double * signal, double dt, double rand_fr=150.0+100.0*rand10())
{
	int i;
	double mx;

	ACore c;
	c.baseFr=rand_fr;
	c.numAmpl=20+rand()%10;
	c.t1=0;
	c.t2=N*dt;
	c.t0=(c.t1+c.t2)*0.5;


	mx=1.0+2.0*rand10();
	for(i=0;i<c.numAmpl;i++)
	{
		c.A[i]=0.5*(2.0+3.0*rand10())/Max(0.5,fabs(mx-i))+0.5*(3.0+2.0*rand10())/Max(0.5,fabs(4.0*mx-i)); // (2.0+rand10())/log(double(i+2));
		c.B[i]=0.5*(2.0+3.0*rand10())/Max(0.5,fabs(mx-i))+0.5*(3.0+2.0*rand10())/Max(0.5,fabs(4.0*mx-i)); // (1.0+rand10())/log(double(i+2));
	}
	for(;i<ACORE_MAX_ELEM;i++)
	{
		c.A[i]=0;
		c.B[i]=0;
	}

	GenerateSignalStatic(N,signal,dt,&c);
}



void RandomSignalDynamic(int N, double * signal, double dt, double rand_fr=150.0+100.0*rand10())
{
	int i,j;
	double mx;

	ACoresSignal cs;
	cs.length=N*dt;
	cs.numCores=int(N/8192); // 2 // 4+rand()%3 // N*dt*3.0
	if(cs.numCores<2)cs.numCores=2;

	cs.cores=new ACore[cs.numCores]; // malloc(cs.numCores*sizeof(ACore))
	for(j=0;j<cs.numCores;j++)
	{
		// !!!! Никакого согласования по фазам. Хорошо, что и так сгенерируется
		cs.cores[j].baseFr=180.0+50.0*rand10();
		cs.cores[j].numAmpl=20+rand()%10;
		cs.cores[j].t1=j*8192.0*dt; // N*dt*j/double(cs.numCores);
		cs.cores[j].t2=(j+1)*8192.0*dt; // N*dt*(j+1)/double(cs.numCores);
		cs.cores[j].t0=(cs.cores[j].t1+cs.cores[j].t2)*0.5;
	
		mx=5.0+5.0*rand10();
		for(i=0;i<cs.cores[j].numAmpl;i++)
		{
			cs.cores[j].A[i]=0.5*(2.0+3.0*rand10())/Max(0.5,fabs(mx-i))+0.5*(3.0+2.0*rand10())/Max(0.5,fabs(4.0*mx-i));
			cs.cores[j].B[i]=0.5*(2.0+3.0*rand10())/Max(0.5,fabs(mx-i))+0.5*(3.0+2.0*rand10())/Max(0.5,fabs(4.0*mx-i));
		}
		for(;i<ACORE_MAX_ELEM;i++)
		{
			cs.cores[j].A[i]=0;
			cs.cores[j].B[i]=0;
		}
	}

	GenerateSignalDynamic(N,signal,dt,&cs);
	delete cs.cores; // free(cs.cores)
}



void ModelSignalDynamic(int N, double * signal, double dt)
{
	int i,j;

	ACoresSignal cs;
	cs.length=N*dt;
	cs.numCores=2;

	cs.cores=new ACore[cs.numCores]; // malloc(cs.numCores*sizeof(ACore))
	for(j=0;j<cs.numCores;j++)
	{
		// !!!! Никакого согласования по фазам. Хорошо, что и так сгенерируется
		cs.cores[j].baseFr=200.0+100.0*j;//+0.0*j; // 200.0
		cs.cores[j].numAmpl=15;
		cs.cores[j].t1=j*N*dt;
		cs.cores[j].t2=j*N*dt;
		cs.cores[j].t0=j*N*dt;

		for(i=0;i<cs.cores[j].numAmpl;i++)
		{
			cs.cores[j].A[i]=3.0+0.0*j;//+j*0.0; // 10.0+j*5.0
			cs.cores[j].B[i]=4.0+0.0*j;//+j*0.0; // 10.0+j*5.0
		}
		for(;i<ACORE_MAX_ELEM;i++)
		{
			cs.cores[j].A[i]=0;
			cs.cores[j].B[i]=0;
		}
	}

	GenerateSignalDynamic(N,signal,dt,&cs);
	delete cs.cores; // free(cs.cores)
}




// пока берем старое, надо сделать заново, с качественным выделением пиков и 
// грамотным выделением нескольких ак ядер из массива пиков
double GetBaseFreq(int N, double * signal, double * re, double * im, double dt)
{
	//////////////////////
	// !!! Вынести массив ampl в аргументы всех функций вместе с re и im ?????
	double * ampl  = (double*)malloc(N*sizeof(double));
	//////////////////////
	int * peaks = (int*)malloc((N/4)*sizeof(int));
	ampl_count(N,re,im,ampl);
	int numP=FindPeaks(ampl,N,peaks,N/4);

	double res = FindACore(peaks,numP,ampl,N);
	free(peaks);
	free(ampl);
	return res/(dt*N);
}



// нормальное интегрирование?
// будем делать трапецией
void GetAmplitudes(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB) // maxAB=30 для наших рандомных сигналов
																// можно maxAB = ACORE_MAX_ELEM
// maxAB - сравнивать полученные амплитуды с "нулевым шумом" сигнала и определять количество значимых
// maxAB - по указателю? Определять внути функции?
// Пока что оставим, позднее надо реализовать
{

	//////////////////
	// вставка из ACoreGenerateEX + поправки
	double fr=baseFr*N*dt;
	int i,j;
	int nf=maxAB; // потом будет атоопределение maxAB, а nf - предельное количество возможных амплитуд 
							// + нужен размер массивов A и B
						// + автозануление лишних амплитуд в A и B
	//int nn=double(N-1)/fr * int(fr); // предельный отсчет для разложения
	double nn=int(double(N-1)*dt * baseFr); // кол-во полных периодов, укладывающихся в сигнал
	nn=nn/(baseFr*dt); // а здесь уже номер отсчета окончания
	
	for(j=0;j<nf;j++)
	{
		A[j]=0;
		B[j]=0;
		for(i=0;i<nn;i++) // простейшее интегрирование, поменять на лучшие методы!
		{
			A[j]+=cos(PI_2*((j+1))*i/double(N)*fr) * signal[i];
			B[j]+=sin(PI_2*((j+1))*i/double(N)*fr) * signal[i];
		}
		// i==nn // 
		A[j]+=cos(PI_2*(j+1)*int(nn)/double(N)*fr) * signal[int(nn)] * (nn-int(nn));
		B[j]+=sin(PI_2*(j+1)*int(nn)/double(N)*fr) * signal[int(nn)] * (nn-int(nn));

		A[j]*=2.0/double(nn);
		B[j]*=2.0/double(nn);
	}
	//////////////////
}

// трапеция + синусы
void GetAmplitudesEX(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB)
{
	double fr=baseFr*N*dt;
	int i,j;
	int nf=maxAB;
	double nn=int(double(N-1)*dt * baseFr); 
	nn=nn/(baseFr*dt); 

	for(j=0;j<nf;j++)
	{
		A[j]=0;
		B[j]=0;
		for(i=0;i<nn;i++) // простейшее интегрирование, поменять на лучшие методы!
		{
			double a,b,c,d;
			double a1,a2,a3;
			double t0,t1;			
			t0=i*dt;
			t1=(i+1)*dt;

			a=(signal[i+1]-signal[i])/dt;
			b=signal[i]*(i+1)-signal[i+1]*i;
			c=( cos(PI_2*(j+1)*baseFr*(i+1)*dt) - cos(PI_2*(j+1)*baseFr*(i)*dt) )/dt;
			d= cos(PI_2*(j+1)*baseFr*(i)*dt)*(i+1) - cos(PI_2*(j+1)*baseFr*(i+1)*dt)*i;
			a1=a*c/3.0;
			a2=(b*c+a*d)*0.5;
			a3=b*d;
			A[j]+=(a1*t1*t1*t1+a2*t1*t1+a3*t1) - (a1*t0*t0*t0+a2*t0*t0+a3*t0);

			c=( sin(PI_2*(j+1)*baseFr*(i+1)*dt) - sin(PI_2*(j+1)*baseFr*(i)*dt) )/dt;
			d= sin(PI_2*(j+1)*baseFr*(i)*dt)*(i+1) - sin(PI_2*(j+1)*baseFr*(i+1)*dt)*i;
			a1=a*c/3.0;
			a2=(b*c+a*d)*0.5;
			a3=b*d;
			B[j]+=(a1*t1*t1*t1+a2*t1*t1+a3*t1) - (a1*t0*t0*t0+a2*t0*t0+a3*t0);
		}

		if(nn < N-1)
		{
			double a,b,c,d;
			double a1,a2,a3;
			double t0,t1;
			t0=int(nn)*dt;
			t1=nn*dt;

			a=(signal[i+1]-signal[i])/dt;
			b=signal[i]*(i+1)-signal[i+1]*i;
			c=( cos(PI_2*(j+1)*baseFr*(i+1)*dt) - cos(PI_2*(j+1)*baseFr*(i)*dt) )/dt;
			d= cos(PI_2*(j+1)*baseFr*(i)*dt)*(i+1) - cos(PI_2*(j+1)*baseFr*(i+1)*dt)*i;
			a1=a*c/3.0;
			a2=(b*c+a*d)*0.5;
			a3=b*d;
			A[j]+=(a1*t1*t1*t1+a2*t1*t1+a3*t1) - (a1*t0*t0*t0+a2*t0*t0+a3*t0);

			c=( sin(PI_2*(j+1)*baseFr*(i+1)*dt) - sin(PI_2*(j+1)*baseFr*(i)*dt) )/dt;
			d= sin(PI_2*(j+1)*baseFr*(i)*dt)*(i+1) - sin(PI_2*(j+1)*baseFr*(i+1)*dt)*i;
			a1=a*c/3.0;
			a2=(b*c+a*d)*0.5;
			a3=b*d;
			B[j]+=(a1*t1*t1*t1+a2*t1*t1+a3*t1) - (a1*t0*t0*t0+a2*t0*t0+a3*t0);
		}

		A[j]*=2.0/(nn*dt);
		B[j]*=2.0/(nn*dt);
	}
}



void GetAmplitudesEX2(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB)
{
	double fr=baseFr*N*dt;
	int i,j;
	int nf=maxAB;
	double nn=int(double(N-1)*dt * baseFr);
	nn=nn/(baseFr*dt);

	for(j=0;j<nf;j++)
	{
		A[j]=0;
		B[j]=0;
		for(i=0;i<nn;i++) // простейшее интегрирование, поменять на лучшие методы!
		{
			double a,b,c;
			double t0,t1;			
			t0=i*dt;
			t1=(i+1)*dt;
			a=(signal[i+1]-signal[i])/dt;
			b=signal[i]*(i+1)-signal[i+1]*i;
			c=PI_2*(j+1)*baseFr;
			
			A[j]+=(a*cos(c*t1)/(c*c)+(a*t1+b)*sin(c*t1)/c) - (a*cos(c*t0)/(c*c)+(a*t0+b)*sin(c*t0)/c);

			B[j]+=(a*sin(c*t1)/(c*c)-(a*t1+b)*cos(c*t1)/c) - (a*sin(c*t0)/(c*c)-(a*t0+b)*cos(c*t0)/c);
		}

		if(nn < N-1)
		{
			double a,b,c;
			double t0,t1;			
			t0=((int)nn)*dt;
			t1=nn*dt;
			a=(signal[i+1]-signal[i])/dt;
			b=signal[i]*(i+1)-signal[i+1]*i;
			c=PI_2*(j+1)*baseFr;
			
			A[j]+=(a*cos(c*t1)/(c*c)+(a*t1+b)*sin(c*t1)/c) - (a*cos(c*t0)/(c*c)+(a*t0+b)*sin(c*t0)/c);

			B[j]+=(a*sin(c*t1)/(c*c)-(a*t1+b)*cos(c*t1)/c) - (a*sin(c*t0)/(c*c)-(a*t0+b)*cos(c*t0)/c);
		}

		A[j]*=2.0/(nn*dt);
		B[j]*=2.0/(nn*dt);
	}
}




void GetAmplitudesSPEC(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB)
{
	int i,j;

	int m=int(double(N-1)*dt * baseFr); // кол-во полных периодов, укладывающихся в сигнал
	double _l=double(m)/baseFr; // _l=m*l0=m/w0
	double w=(1.0/dt)/double(N);// шаг частоты в сигнале, w=44100/N

	double c1,c2;


	for(j=0;j<maxAB;j++)
	{
		A[j]=0;
		B[j]=0;
		for(i=1;i<N/2;i++) // i=0 или i=1 одинаково правильно
		{
			c1=i*w;
			c2=baseFr*double(j+1);  // double(m)*

			A[j]+=re[i]*(baseFr/(PI*m))*(c1/(c1*c1-c2*c2))*(sin(PI_2*c1*_l)) 
				+ im[i]*(baseFr/(PI*m))*(c1/(c1*c1-c2*c2))*(cos(PI_2*c1*_l)-1);
			B[j]+=re[i]*(baseFr/(PI*m))*(c2/(c1*c1-c2*c2))*(1-cos(PI_2*c1*_l))
				+ im[i]*(baseFr/(PI*m))*(c2/(c1*c1-c2*c2))*(sin(PI_2*c1*_l));



//			A[j]+=re[i]*sin(PI_2*baseFr*double(j+1)*N*dt)*double(j+1)*baseFr/(baseFr*double(j+1)*double(j+1)*baseFr-double(i)*double(i)/(dt*N*dt*N))/(PI*N*dt)
//				+im[i]*(1.0-cos(PI_2*baseFr*double(j+1)*N*dt))*baseFr*double(j+1)/(baseFr*double(j+1)*double(j+1)*baseFr-double(i)*double(i)/(dt*N*dt*N))/(PI*N*dt);
//			B[j]+=re[i]*(cos(PI_2*baseFr*double(j+1)*N*dt)-1)*double(i)/(dt*N)/(baseFr*double(j+1)*double(j+1)*baseFr-double(i)*double(i)/(dt*N*dt*N))/(PI*N*dt)
//				+im[i]*sin(PI_2*baseFr*double(j+1)*N*dt)*double(i)/(dt*N)/(baseFr*double(j+1)*double(j+1)*baseFr-double(i)*double(i)/(dt*N*dt*N))/(PI*N*dt);
		}
		A[j]*=(2.0/double(N));//(2.0/1024.0);//
		B[j]*=-(2.0/double(N));//-(2.0/1024.0);//
	//	double tmp=A[j];
	//	A[j]=B[j];
	//	B[j]=tmp;
		//A[j]*=2.0/_l;
		//B[j]*=2.0/_l;
	}

}



void GetAmplitudesSPECEX(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB)
{
	int i,j,k;

	int m=int(double(N-1)*dt * baseFr); // кол-во полных периодов, укладывающихся в сигнал
	double _l=double(m)/baseFr; // _l=m*l0=m/w0
	double w=(1.0/dt)/double(N);// шаг частоты в сигнале, w=44100/N

	double c1,c2;
	double d1,d2;


	for(j=0;j<maxAB;j++)
	{
		A[j]=0;
		B[j]=0;
		for(i=1;i<N/2;i++) // i=0 или i=1 одинаково правильно
		{
			d1=d2=0;
			c2=baseFr*double(j+1);  // double(m)*
			for(k=0;k<60;k++)
			{
				c1=(i+k*N)*w;
				A[j]+=re[i]*(baseFr/(PI*m))*(c1/(c1*c1-c2*c2))*(sin(PI_2*c1*_l)) 
					+ im[i]*(baseFr/(PI*m))*(c1/(c1*c1-c2*c2))*(cos(PI_2*c1*_l)-1);
				B[j]+=re[i]*(baseFr/(PI*m))*(c2/(c1*c1-c2*c2))*(1-cos(PI_2*c1*_l))
					+ im[i]*(baseFr/(PI*m))*(c2/(c1*c1-c2*c2))*(sin(PI_2*c1*_l));
				c1=(N-i+k*N)*w;
				A[j]+=re[i]*(baseFr/(PI*m))*(c1/(c1*c1-c2*c2))*(sin(PI_2*c1*_l)) 
					+ im[i]*(baseFr/(PI*m))*(c1/(c1*c1-c2*c2))*(cos(PI_2*c1*_l)-1);
				B[j]+=re[i]*(baseFr/(PI*m))*(c2/(c1*c1-c2*c2))*(1-cos(PI_2*c1*_l))
					+ im[i]*(baseFr/(PI*m))*(c2/(c1*c1-c2*c2))*(sin(PI_2*c1*_l));
			}


		}
		A[j]*=(2.0/double(N));//(2.0/1024.0);//
		B[j]*=-(2.0/double(N));//-(2.0/1024.0);//
	//	double tmp=A[j];
	//	A[j]=B[j];
	//	B[j]=tmp;
		//A[j]*=2.0/_l;
		//B[j]*=2.0/_l;
	}

}

void GetAmplitudesOLD(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB)
/// !!! Идет перезапись re и im!
{
	ACore cr;
	double * gen = new double[N];
	int i,j;
	double fr_d=baseFr*(dt*N);
	for(i=0;i<maxAB;i++)
	{
		A[i]=re[int((i+1)*fr_d+0.5)];
		B[i]=im[int((i+1)*fr_d+0.5)];
		cr.A[i]=A[i]*2/double(N);
		cr.B[i]=B[i]*2/double(N);
	}
	cr.baseFr=baseFr;
	cr.numAmpl=maxAB;

	for(j=0;j<1;j++) // ??????
	{
		GenerateSignalStatic(N,gen,dt,&cr);
	
		fft(N,gen,re,im);
		for(i=0;i<maxAB;i++)
		{
			double c_re,c_im,z2_re,z2_im,mod2;
			double tmp_re,tmp_im;
			z2_re=re[int((i+1)*fr_d+0.5)];
			z2_im=im[int((i+1)*fr_d+0.5)];
			mod2=z2_re*z2_re+z2_im*z2_im;
			c_re=A[i]*z2_re+B[i]*z2_im;	
			c_im=B[i]*z2_re-A[i]*z2_im;
			c_re/=mod2;
			c_im/=mod2;
			tmp_re=cr.A[i]*c_re-cr.B[i]*c_im;
			tmp_im=cr.A[i]*c_im+cr.B[i]*c_re;
			cr.A[i]=tmp_re;
			cr.B[i]=tmp_im;
		}
	}

	for(i=0;i<maxAB;i++)
	{
		A[i]=cr.A[i];
		B[i]=cr.B[i];
	}

	delete gen;
}

void GetACore(int N, double * signal, double * re, double * im, double dt,
			  ACore * core)
{
	core->baseFr = GetBaseFreq(N,signal,re,im,dt);
	if(core->baseFr==0)//core->baseFr=200; // ТОЛЬКО НА ВРЕМЯ ОТЛАДКИ
	{
		core->numAmpl=0;
		return;
	}
	core->baseFr = DichotomyBaseFreq(N,signal,re,im,dt,
		core->baseFr-0.5/double(dt*N),core->baseFr+0.5/double(dt*N));
	GetAmplitudesOPT(N,signal,/*re,im,*/dt,core->baseFr,core->A,core->B,AMPL_MAX_AB); // 30!!!!
	 // SPECEX?
	core->numAmpl=AMPL_MAX_AB; // Переделать вместе с GetAmplitudes!!!!
	for(int i=AMPL_MAX_AB;i<ACORE_MAX_ELEM;i++)
		core->A[i]=core->B[i]=0;
}




void SignalModel(int N, double * signal, double * re, double * im, double dt,
				double baseFr, double * result)
{
	ACore cr;
	cr.baseFr = baseFr;
	GetAmplitudesSPEC(N,signal,re,im,dt,cr.baseFr,cr.A,cr.B,AMPL_MAX_AB); // 30!!!!
	cr.numAmpl=AMPL_MAX_AB; // Переделать вместе с GetAmplitudes!!!!
	for(int i=AMPL_MAX_AB;i<ACORE_MAX_ELEM;i++)
		cr.A[i]=cr.B[i]=0;

	GenerateSignalStatic(N,result,dt,&cr);
}



void SignalRest(int N, double * signal, double * re, double * im, double dt,
				double baseFr, double * rest, bool spec)
{
	int i;
	ACore cr;
	cr.baseFr = baseFr;
	if(spec)
		//GetAmplitudesSPEC(N,signal,re,im,dt,cr.baseFr,cr.A,cr.B,AMPL_MAX_AB); // 30!!!!
		GetAmplitudesOPT(N,signal,dt,baseFr,cr.A,cr.B,50);
	else
		GetAmplitudesEX2(N,signal,re,im,dt,cr.baseFr,cr.A,cr.B,AMPL_MAX_AB); // 30!!!!
		//GetAmplitudesOPT(N,signal,dt,baseFr,cr.A,cr.B,50);
	cr.numAmpl=AMPL_MAX_AB; // Переделать вместе с GetAmplitudes!!!!
	for(i=AMPL_MAX_AB;i<ACORE_MAX_ELEM;i++)
		cr.A[i]=cr.B[i]=0;

	GenerateSignalStatic(N,rest,dt,&cr);

	for(i=0;i<N;i++)
		rest[i]=signal[i]-rest[i];
}



double SignalRestMod(int N, double * signal, double * re, double * im, double dt,
				double baseFr, bool spec)
{
	double res=0;
	double * rest = (double*)malloc(N*sizeof(double));

	SignalRest(N,signal,re,im,dt,baseFr,rest,spec);
	
	for(int i=0;i<N;i++) // и снова линейное интегр-е, заменить на трапеции
		res+=dt*rest[i]*rest[i]; // квадратичная норма
	res=sqrt(res);

	free(rest);
	return res;
}



void CreateACores(int N, int T, double * signal, double dt, ACoresSignal * cs)
{
	cs->length=N*dt;
	cs->numCores=N/T;
	cs->coreMethod=ACORE_METHOD_DEFAULT;
/*	if(cs->cores!=NULL)
	{
		free(cs->cores);
		cs->cores=NULL;
	}//*/
	cs->cores=new ACore[cs->numCores]; // malloc

	double * re = (double*)malloc(T*sizeof(double));
	double * im = (double*)malloc(T*sizeof(double));
	for(int i=0;i<cs->numCores;i++)
	{
		fft(T,signal+i*T,re,im);
		GetACore(T,signal+i*T,re,im,dt,cs->cores+i);
		cs->cores[i].t1=i*double(T)*dt;
		cs->cores[i].t2=(i+1)*double(T)*dt;
		cs->cores[i].t0=(i+0.5)*double(T)*dt;
	}

	free(re);
	free(im);

}


void CreateACoresEX(int N, double * signal, double dt, ACoresSignal * cs)
{
	int i,j,k;
	int T=1024;
	double fr;

	cs->length=N*dt;
	cs->numCores=0;
	cs->coreMethod=ACORE_METHOD_DEFAULT;
	cs->cores=new ACore[4*(N/T)]; // malloc

	double * re = (double*)malloc(T*sizeof(double));
	double * im = (double*)malloc(T*sizeof(double));
	for(i=0;i<N/T;i++)
	{
		fft(T,signal+i*T,re,im);
		fr=GetBaseFreq(T,signal+i*T,re,im,dt);
		if(fr<216)//fr<216)
		{
			fft(512,signal+i*T,re,im);
			cs->cores[cs->numCores].baseFr=DichotomyBaseFreq(512,signal+i*T,re,im,dt,
				fr-1.0/double(dt*1024),fr+1.0/double(dt*1024));
			if(fabs(cs->cores[cs->numCores].baseFr)<0.1)
			{
				cs->cores[cs->numCores].numAmpl=0;
			}
			else
			{
				//GetAmplitudesOPT(512,signal+i*T,dt,cs->cores[cs->numCores].baseFr,cs->cores[cs->numCores].A,cs->cores[cs->numCores].B,AMPL_MAX_AB); // 30!!!!
				GetAmplitudesSPEC(512,signal+i*T,re,im,dt,cs->cores[cs->numCores].baseFr,cs->cores[cs->numCores].A,cs->cores[cs->numCores].B,AMPL_MAX_AB); // 30!!!!
				// SPECEX?
				cs->cores[cs->numCores].numAmpl=AMPL_MAX_AB; // Переделать вместе с GetAmplitudes!!!!
				for(k=AMPL_MAX_AB;k<ACORE_MAX_ELEM;k++)
					cs->cores[cs->numCores].A[k]=cs->cores[cs->numCores].B[k]=0;
			}

			cs->cores[cs->numCores].t1=i*double(T)*dt;
			cs->cores[cs->numCores].t2=(i*double(T)+512)*dt;
			cs->cores[cs->numCores].t0=(i*double(T)+256)*dt;

			cs->numCores++;
/////////////////////////////////////////////////////
			fft(512,signal+i*T+512,re,im);
			cs->cores[cs->numCores].baseFr=DichotomyBaseFreq(512,signal+i*T+512,re,im,dt,
				fr-1.0/double(dt*1024),fr+1.0/double(dt*1024));
			if(fabs(cs->cores[cs->numCores].baseFr)<0.1)
			{
				cs->cores[cs->numCores].numAmpl=0;
			}
			else
			{
				//GetAmplitudesOPT(512,signal+i*T+512,dt,cs->cores[cs->numCores].baseFr,cs->cores[cs->numCores].A,cs->cores[cs->numCores].B,AMPL_MAX_AB); // 30!!!!
				GetAmplitudesSPEC(512,signal+i*T+512,re,im,dt,cs->cores[cs->numCores].baseFr,cs->cores[cs->numCores].A,cs->cores[cs->numCores].B,AMPL_MAX_AB); // 30!!!!
				// SPECEX?
				cs->cores[cs->numCores].numAmpl=AMPL_MAX_AB; // Переделать вместе с GetAmplitudes!!!!
				for(k=AMPL_MAX_AB;k<ACORE_MAX_ELEM;k++)
					cs->cores[cs->numCores].A[k]=cs->cores[cs->numCores].B[k]=0;
			}
			
			cs->cores[cs->numCores].t1=(i*double(T)+512)*dt;
			cs->cores[cs->numCores].t2=(i*double(T)+1024)*dt;
			cs->cores[cs->numCores].t0=(i*double(T)+768)*dt;

			cs->numCores++;
		}
		else
		{
			for(j=0;j<4;j++)
			{
				fft(256,signal+i*T+j*256,re,im);
				cs->cores[cs->numCores].baseFr=DichotomyBaseFreq(256,signal+i*T+j*256,re,im,dt,
				fr-1.0/double(dt*1024),fr+1.0/double(dt*1024));
				if(fabs(cs->cores[cs->numCores].baseFr)<0.1)
				{
					cs->cores[cs->numCores].numAmpl=0;
				}
				else
				{
					//GetAmplitudesOPT(256,signal+i*T+j*256,dt,cs->cores[cs->numCores].baseFr,cs->cores[cs->numCores].A,cs->cores[cs->numCores].B,AMPL_MAX_AB); // 30!!!!
					GetAmplitudesSPEC(512,signal+i*T+j*256,re,im,dt,cs->cores[cs->numCores].baseFr,cs->cores[cs->numCores].A,cs->cores[cs->numCores].B,AMPL_MAX_AB); // 30!!!!
					// SPECEX?
					cs->cores[cs->numCores].numAmpl=AMPL_MAX_AB; // Переделать вместе с GetAmplitudes!!!!
					for(k=AMPL_MAX_AB;k<ACORE_MAX_ELEM;k++)
						cs->cores[cs->numCores].A[k]=cs->cores[cs->numCores].B[k]=0;
				}

				cs->cores[cs->numCores].t1=(i*double(T)+j*256.0)*dt;
				cs->cores[cs->numCores].t2=(i*double(T)+j*256.0+256.0)*dt;
				cs->cores[cs->numCores].t0=(i*double(T)+j*256.0+128.0)*dt;

				cs->numCores++;
			}
		}

	}

	free(re);
	free(im);
}




double AdjustBaseFreq(int N, double * signal, double dt,
					 double oldFr)
					 // double * re, double * im, 
{
	ACore c1,c2;
	double * re = (double*)malloc((N/2)*sizeof(double));
	double * im = (double*)malloc((N/2)*sizeof(double));

	fft(N/2,signal,re,im);
	//GetACore(N/2,signal,re,im,dt,&c1);
	GetAmplitudes(N/2,signal,re,im,dt,oldFr,c1.A,c1.B);
	fft(N/2,signal+(N/2),re,im);
	//GetACore(N/2,signal+(N/2),re,im,dt,&c2);
	GetAmplitudes(N/2,signal+(N/2),re,im,dt,oldFr,c2.A,c2.B);
	
	double time=(N/2)*dt;
	double ang1,ang2;
	double angD;
	double err;

	// пока для первой гармоники.
	ang1=Angle(c1.A[0],c1.B[0]);
	ang2=Angle(c2.A[0],c2.B[0]);
	angD=PI_2*time*oldFr + ang1;
	
	//while(fabs(ang2-angD)>PI+0.001)angD-=PI_2; // 0.001 - для безопасности
	// считаем ошибку маленькой, близко
	while(angD>ang2)angD-=PI_2;
	if(fabs(ang2-angD)>fabs(angD+PI_2-ang2)) // fabs вроде не нужны
		err=angD+PI_2-ang2;
	else
		err=ang2-angD;

	oldFr+=err/PI_2*time; // помним, что частота передается в ф-ю по значению.


	free(re);
	free(im);


	return oldFr;
}


double dich_eps=0.000000001;
double DichotomyBaseFreq(int N, double * signal, double * re, double * im, double dt,
						  double fr_start, double fr_end)
{
//	int i;

	/*// MODIF
	double d=(fr_start+fr_end)*0.00125; //(1/200-я от средн. частоты)
	int n=(fr_end-fr_start)/d + 1;
	d=(fr_end-fr_start)/double(n);
	double fr_min=fr_start,re_min=SignalRestMod(N,signal,re,im,dt,fr_min),retmp;
	for(i=1;i<n;i++)
		if( (retmp=SignalRestMod(N,signal,re,im,dt,fr_start+i*d))<re_min )
		{
			fr_min=fr_start+i*d;
			re_min=retmp;
		}
	//*///

	double frs,fre,
		fr1,fr2,fr3;
	frs=fr_start;//fr_min-d;//fr_start;
	fre=fr_end;//fr_min+d;//fr_end;
	fr2=(frs+fre)/2.0;
	fr1=(frs+fr2)/2.0;
	fr3=(fr2+fre)/2.0;

	double res,ree,re1,re2,re3;
	res=SignalRestMod(N,signal,re,im,dt,frs);
	ree=SignalRestMod(N,signal,re,im,dt,fre);
	re1=SignalRestMod(N,signal,re,im,dt,fr1);
	re2=SignalRestMod(N,signal,re,im,dt,fr2);//re_min;//SignalRestMod(N,signal,re,im,dt,fr2);
	re3=SignalRestMod(N,signal,re,im,dt,fr3);

	bool spec;
	while(fabs(frs-fre)>dich_eps)
	{
		spec=(fre-frs>0.01)?false:true; // 0.01

		if(re1<re2 && re1<re3)
		{
			fre=fr2;
			fr2=fr1;
			fr1=(frs+fr2)/2.0;
			fr3=(fr2+fre)/2.0;
			ree=re2;
			re2=re1;
			re1=SignalRestMod(N,signal,re,im,dt,fr1,spec);
			re3=SignalRestMod(N,signal,re,im,dt,fr3,spec);
		}
		else if(re2<re1 && re2< re3)
		{
			frs=fr1;
			fre=fr3;
			fr1=(frs+fr2)/2.0;
			fr3=(fr2+fre)/2.0;
			res=re1;
			ree=re3;
			re1=SignalRestMod(N,signal,re,im,dt,fr1,spec);
			re3=SignalRestMod(N,signal,re,im,dt,fr3,spec);
		}
		else
		{
			frs=fr2;
			fr2=fr3;
			fr1=(frs+fr2)/2.0;
			fr3=(fr2+fre)/2.0;
			res=re2;
			re2=re3;
			re1=SignalRestMod(N,signal,re,im,dt,fr1,spec);
			re3=SignalRestMod(N,signal,re,im,dt,fr3,spec);
		}
	}

	return fr2;
}










bool ModifFr_Find(ACore * c1, ACore * c2, double * coef)
{
	int i;
	int er1=0;
	double mamp=0;

	for(i=0;i<c1->numAmpl;i++)
		if(c1->A[i] > mamp)
			mamp=c1->A[i];

	for(i=0;i<c1->numAmpl;i++)
	{
		double d=fabs(c2->A[i] - c1->A[i]);
		if(d>0.05*mamp)
			return false;
		else if(d>0.03*mamp)
			er1++;
	}
	if(er1<=3)
	{
		*coef=c2->baseFr/c1->baseFr;
		return true;
	}

	return false;
}


bool ModifAmp_Find(ACore * c1, ACore * c2, double * coef)
{
	int i;
	int er1=0;
	double mamp=0;

	double cf=0;
	int n=0;

	if(fabs((c1->baseFr-c2->baseFr)/c1->baseFr) >0.03 )
		return false;

	for(i=0;i<c1->numAmpl;i++)
		if(c1->A[i] > mamp)
			mamp=c1->A[i];
	for(i=0;i<c1->numAmpl;i++)
		if(c1->A[i] > 0.2*mamp)
		{
			cf+=c2->A[i]/c1->A[i];
			n++;
		}
	cf/=n;


	for(i=0;i<c1->numAmpl;i++)
	{
		double d=fabs(c2->A[i] - cf*c1->A[i]);
		if(d>0.05*cf*mamp)
			return false;
		else if(d>0.03*cf*mamp)
			er1++;
	}
	if(er1<=3)
	{
		*coef=cf;
		return true;
	}

	return false;
}




ACoreArray::ACoreArray()
{
	maxblocks=32;
	cores=new ACore*[32];
	cores[0]=new ACore[ACARRAY_BLOCKSIZE];
	numblocks=1;
	numcores=0;
}
ACoreArray::~ACoreArray()
{
	for(int i=0; i<numblocks;i++)
		delete cores[i];
	delete cores;
}
int ACoreArray::AddCore(ACore *c)
{
	if(numcores==maxcores())
	{
		cores[numblocks]=new ACore[ACARRAY_BLOCKSIZE]; // if ( ... == NULL) return 0;
		numblocks++;

		if(numblocks==maxblocks)
		{
			ACore** tmp=new ACore*[maxblocks+32]; // if ( ... == NULL) return 0;
			memcpy(tmp,cores,maxblocks*sizeof(ACore*));
			maxblocks+=32;
		}
	}
	ACCopy(&(cores[numcores/ACARRAY_BLOCKSIZE][numcores%ACARRAY_BLOCKSIZE]),c);
	numcores++;

	return 1;
}
ACore * ACoreArray::GetCore(int n)
{
	if( n<0 || n>=numcores )
		return NULL;
	return &(cores[n/ACARRAY_BLOCKSIZE][n%ACARRAY_BLOCKSIZE]);
}
void CreateACoresPeriodic(int N, double * signal, double dt, ACoreArray * ca)
{
	double re[1024],im[1024];
	ACore cr;

	if(N<1024)return;

	int p_st=0;
	double fr;
	int sz;

	while(p_st<N)
	{
		if(p_st<=N-1024)
		{
			fft(1024,&signal[p_st],re,im);
			fr=GetBaseFreq(1024,&signal[p_st],re,im,dt);
		}
		else
		{
			fft(1024,&signal[N-1024],re,im);
			fr=GetBaseFreq(1024,&signal[N-1024],re,im,dt);
		}

		if(fr<216)
			sz=512;
		else
			sz=256;

		if(p_st<=N-sz)
		{
			fft(sz,&signal[p_st],re,im);
			cr.baseFr=DichotomyBaseFreq(sz,&signal[p_st],re,im,dt,fr-1.0/double(dt*1024),fr+1.0/double(dt*1024));
			GetAmplitudesSPEC(sz,&signal[p_st],re,im,dt,cr.baseFr,cr.A,cr.B); // SPECEX?
			cr.numAmpl=AMPL_MAX_AB;
			cr.t1=p_st*dt;
			p_st+=int(1.0/(dt*cr.baseFr))+1;
			cr.t2=p_st*dt;
			cr.t0=(cr.t1+cr.t2)*0.5;
		}
		else
		{
			fft(sz,&signal[N-sz],re,im);
			cr.baseFr=DichotomyBaseFreq(sz,&signal[N-sz],re,im,dt,fr-1.0/double(dt*1024),fr+1.0/double(dt*1024));
			GetAmplitudesSPEC(sz,&signal[N-sz],re,im,dt,cr.baseFr,cr.A,cr.B); // SPECEX?
			cr.numAmpl=AMPL_MAX_AB;
			cr.t1=(N-sz)*dt;
			cr.t2=N*dt;
			cr.t0=(cr.t1+cr.t2)*0.5;
			p_st=N;//p_st+=int(1.0/(dt*cr.baseFr))+1;
		}

		ca->AddCore(&cr);
	}
}


void GetACoreInPoint(int N, double * signal, int pt, double dt, ACore * core)
{
	double re[1024],im[1024];
	int p_st;
	double fr;
	int sz;

	if(N<1024) return;


	if(pt<512)
		p_st=0;
	else if(pt>N-512)
		p_st=N-1024;
	else
		p_st=pt-512;
	fft(1024,&signal[p_st],re,im);
	fr=GetBaseFreq(1024,&signal[p_st],re,im,dt);
	if(fr<216)
		sz=512;
	else
		sz=256;

	if(pt<sz/2)
		p_st=0;
	else if(pt>N-sz/2)
		p_st=N-sz;
	else
		p_st=pt-sz/2;
	fft(sz,&signal[p_st],re,im);
	core->baseFr=DichotomyBaseFreq(sz,&signal[p_st],re,im,dt,fr-1.0/double(dt*1024),fr+1.0/double(dt*1024));
	GetAmplitudesSPEC(sz,&signal[p_st],re,im,dt,core->baseFr,core->A,core->B); // SPECEX?
	core->numAmpl=AMPL_MAX_AB;
	core->t0=core->t1=core->t2=pt*dt;

	double ang,mod,angadd=core->baseFr*dt*(pt-p_st);
	for(int i=0;i<AMPL_MAX_AB;i++)
	{
		mod=sqrt(core->A[i]*core->A[i] + core->B[i]*core->B[i]);
		ang=Angle(core->A[i],core->B[i]);
		ang+=angadd*(i+1);
		core->A[i]=mod*cos(ang);
		core->B[i]=mod*sin(ang);
	}
}


#endif