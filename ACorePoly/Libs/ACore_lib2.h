#ifndef ACORE_LIB2_H
#define ACORE_LIB2_H





#define ACORE_MAX_ELEM 100
struct ACore
{
	double baseFr; // базовая частота
	int    numAmpl;
	double A[ACORE_MAX_ELEM]; // коэфициенты амплитуд косинусов
	double B[ACORE_MAX_ELEM]; // коэфициенты амплитуд синусов
	// A[0],B[0] - амплитуды ПЕРВОЙ гармоники
	// A[i],B[i] - (i+1)-ой гармоники

	double t1,t2; // границы отрезка разложения
	double t0; // точка привязки
};


#define ACORE_METHOD_DEFAULT 1
struct ACoresSignal
{
	ACoresSignal():length(0),numCores(0),cores(nullptr),coreMethod(ACORE_METHOD_DEFAULT) {}
	~ACoresSignal() {if(cores)delete cores; cores=nullptr;}

	double length; // длительность сигнала
	int numCores; // число ядер
	ACore * cores; // массив ядер
	int coreMethod; // метод разложения сигнала на ядра
};


#define ACARRAY_BLOCKSIZE 32
struct ACoreArray
{
	ACore ** cores;
	int maxblocks;
	int numblocks;
	int numcores;
	int maxcores(){return numblocks*ACARRAY_BLOCKSIZE;}

	ACoreArray();
	~ACoreArray();
	int AddCore(ACore * c);
	ACore * GetCore(int n);
};


void ACSave(FILE * f, ACore * c);
void ACSave(const char * fname, ACore * c);
void ACSave(const char * fname, ACoresSignal * cs);

void	ACCopy(ACore * c, ACore * c1);
void	ACMix(ACore * c, ACore * c1, ACore * c2, double d1, double d2);
void	ACNorm(ACore * c);
double Angle(double dCos, double dSin);

void ACConvert(ACoresSignal* cs, ACoreArray* acr)
{
	cs->length=acr->GetCore(acr->numcores-1)->t2;
	cs->coreMethod=ACORE_METHOD_DEFAULT;
	cs->numCores=acr->numcores;
	cs->cores=new ACore[cs->numCores];
	for(int i=0;i<cs->numCores;i++)
		ACCopy(&cs->cores[i],acr->GetCore(i));
}


void GenerateSignalStatic(int N, double * signal, double dt, ACore * core);
// dt - шаг времени, N*dt - длительность сигнала.
void GenerateSignalDynamic(int N, double * signal, double dt, ACoresSignal * cs);
void GenerateSignalDynamicNO(int N, double * signal, double dt, ACoresSignal * cs);


void RandomSignalStatic(int N, double * signal, double dt, double rand_fr);  // rand_fr - nado li?
void RandomSignalDynamic(int N, double * signal, double dt, double rand_fr); // 

void ModelSignalDynamic(int N, double * signal, double dt);



#define AMPL_MAX_AB 50

double GetBaseFreq(int N, double * signal, double * re, double * im, double dt);
void GetAmplitudes(int N, double * signal, double * re, double * im, double dt,
				   double baseFreq, double * A, double * B, int maxAB=AMPL_MAX_AB);
void GetACore(int N, double * signal, double * re, double * im, double dt,
			  ACore * core);
// core строится без временных меток, только частота/амплитуды

void GetAmplitudesEX(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB=AMPL_MAX_AB);
void GetAmplitudesEX2(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB=AMPL_MAX_AB);
void GetAmplitudesSPEC(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB=AMPL_MAX_AB);
void GetAmplitudesOLD(int N, double * signal, double * re, double * im, double dt,
				   double baseFr, double * A, double * B, int maxAB=AMPL_MAX_AB);


void SignalModel(int N, double * signal, double * re, double * im, double dt,
				double baseFr, double * result);
void SignalRest(int N, double * signal, double * re, double * im, double dt,
				double baseFr, double * rest, bool spec=false);
double SignalRestMod(int N, double * signal, double * re, double * im, double dt,
				double baseFr, bool spec=false);
// напр. baseFr = GetBaseFreq();
// а вообще заточено под дихотомию


// разные методы?
void CreateACores(int N, int T, double * signal, double dt, ACoresSignal * cs);
// обратная к GenerateSignalDynamic
// T - длина отрезка для одного ядра.

// как вариант берем 2048 отсчетов, считаем ~базовую частоту, затем по ядру на виток

void CreateACoresEX(int N, double * signal, double dt, ACoresSignal * cs);




double AdjustBaseFreq(int N, double * signal, double dt,
					 double oldFr); // пока что для стац. сигнала
// берем две половины сигнала и сравниваем фазы

double DichotomyBaseFreq(int N, double * signal, double * re, double * im, double dt,
						  double fr_start, double fr_end);
// double * re, double * im, 

void CreateACoresPeriodic(int N, double * signal, double dt, ACoreArray * ca);
void GetACoreInPoint(int N, double * signal, int pt, double dt, ACore * core);




bool ModifFr_Find(ACore * c1, ACore * c2, double * coef);
bool ModifAmp_Find(ACore * c1, ACore * c2, double * coef);




#endif