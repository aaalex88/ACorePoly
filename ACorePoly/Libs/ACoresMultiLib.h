#pragma once


#include "libs/ACore_lib.cpp"
#include "libs/ACore_lib2.cpp"
#include "libs/ACore_lib3.cpp"


#define MAX_CORES 5
#define MAX_CORE_AMPL 20


#define ALLOCATOR_SIZE (128*1024*1024)
class cStackAllocatorManager;
class cStackAllocator;


struct SBasis
{
	int N;
	int K;

	int cores_numAmpl[MAX_CORES];

	double * ampl;
	double * basis;
};

void GenerateTestSignal2Cores();
void ResolveSegment(ACore * oldCores, ACore * newCores, int N, double dt, double * signal);

//	void FormBasisByCores(basis, oldCores);
//	bool ResolveBasis(basis, N, dt, signal);
//	void FormCoresByBasis(basis, newCores);
//	void GetRest(N, dt, signal, rest, newCores);
//	void AddSignal(N,one_core,rest); // elementary
//	double MaxDiff(newCores,old_freq); // elementary

//	bool TryGetCore(rest,new_cor)
//	void InserCoreToArr(new_cor,newCores); // elementary


