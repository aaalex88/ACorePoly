
#include <stdlib.h>
#include <assert.h>

#include "ACoresMultiLib.h"


class cStackAllocatorManager
{
public:
	static cStackAllocatorManager & Instance()
	{
		static cStackAllocatorManager instance(ALLOCATOR_SIZE);
		return instance;
	}

	void *	Malloc(int size)
	{
		if (size < 0 || (m_ptr+size) > (m_buf+m_bufSize))
		{
			assert((m_ptr+size) <= (m_buf+m_bufSize));
			return NULL;
		}

		void * ptr    = m_ptr; 
		       m_ptr += size; 
		return ptr;
	}
	void	Free(void * ptr)
	{
		if ( ptr >= m_buf && ptr < m_ptr )
		m_ptr = (char*)ptr;
	}

protected:
	cStackAllocatorManager(int bufSize)
	{
		if (!(m_buf = (char*)malloc(bufSize)))
		{
			exit(-1);
		}
		m_bufSize = bufSize;
		m_ptr = m_buf;
	}
	~cStackAllocatorManager()
	{
		free(m_buf);
	}


	int m_bufSize;
	char * m_buf;
	char * m_ptr;

};

class cStackAllocator
{
public:
	cStackAllocator(){m_ptr = (char*)cStackAllocatorManager::Instance().Malloc(0);}
	~cStackAllocator(){cStackAllocatorManager::Instance().Free(m_ptr);}
	void * Malloc(int size){return cStackAllocatorManager::Instance().Malloc(size);}

private:
	char * m_ptr;

};
/*
 *
 *
 *
 */
void ResolveSegment(ACore * oldCores, ACore * newCores, int N, double dt, double * signal)
{
	cStackAllocator al;

	SBasis basis;
	ACore one_core;
	double * rest = (double*)al.Malloc(N*sizeof(double));
	double * one_core = (double*)al.Malloc(N*sizeof(double));
	double old_freq[MAX_CORES];

	int num_cores;
	// num_cores = ...

//	FormBasisByCores(basis, oldCores);
//	ResolveBasis(basis, N, dt, signal);
//	FormCoresByBasis(basis, newCores);

	do
	{
//		GetRest(N, dt, signal, rest, newCores);
		for(int i = 0; i < MAX_CORES && newCores[i].numAmpl != 0; i++)
		{
			old_freq[i] = newCores[i].baseFr;
			// isolated counting
			GenerateSignalStatic(N,one_core,dt,newCores+i);
//			AddSignal(N,one_core,rest); // one_core += rest;
			newCores[i].baseFr = DichotomyBaseFreq(N,one_core,0,0,dt,newCores[i].baseFr-5.0,newCores[i].baseFr+5.0);
			GetAmplitudesOPT(N,one_core,dt,newCores[i].baseFr,newCores[i].A,newCores[i].B,MAX_CORE_AMPL);
		}
	}
	while (0);//(MaxDiff(newCores,old_freq) > 0.00001);

	if (num_cores < MAX_CORES)
	{
		ACore new_cor;
//		if(TryGetCore(rest,new_cor))
		{
//			InserCoreToArr(new_cor,newCores);
		}
	}

}



void FormBasisByCores (SBasis & basis, ACores * cores)
{
}

bool ResolveBasis (SBasis & basis, int N, double dt, double * signal
{
}

void FormCoresByBasis (SBasis & basis, ACores * cores)
{
}

void GetRest (int N, double dt, double * signal, double * rest, ACores * cores)
{
}

void AddSignal(int N, double * signal, double * add_signal)
{
}

double MaxDiff(ACore * cores, double * oldFreq)
{
}

bool TryGetCore(int N, double dt, double * signal, ACore & core)
{
}

void InserCoreToArr(ACore & core, ACore * coreArr)
{
}




void GenerateTestSignal2Cores()
{

	double * arr = new double[441000];
	double dt = 1.0 / 44100.0;

	double a1[] = { 5, 2, 4, 1 };
	double a2[] = { 7, 3, 0.5, 2 };

	for (int i = 0; i < 441000; i++)
	{
		arr[i] = 0;

		for (int j = 0; j < 4; j++)
		{
			arr[i] += a1[j] * (1.0 + 0.5*sin(i*dt))  *  sin( PI_2*(j+1)*(300.0*i*dt + 25.0*(i*dt)*(i*dt)) );
			arr[i] += a2[j] * ( exp( -(i*dt-6.0)*(i*dt-6.0)/2.0 ) ) * sin( PI_2*(j+1)*(650.0*i*dt + 300.0*sin(0.5*i*dt)) );
		}
	}

	SaveWAV("test2cores.wav",441000,44100,441000,arr);

	delete arr;
}