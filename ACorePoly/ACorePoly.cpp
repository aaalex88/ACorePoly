// ACorePoly.cpp: определяет точку входа для консольного приложения.
//


#include "stdafx.h"



#include "libs/FFT_lib2.cpp"
#include "libs/WAV_lib.cpp"
#include "libs/ACore_lib.cpp"
#include "libs/ACore_lib2.cpp"
#include "libs/ACore_lib3.cpp"


#include "PolyLib/ACorePolyLib.h"
#include "Tests/Tests.h"

// #include "gsl\linalg\gsl_linalg.h"


int _tmain(int argc, char* argv[])
{
	/*
	Solution<SignalSolver>::AnalyseSignal<shared_ptr<Signal>,shared_ptr<SegmentDescription>>
										( [](){ return shared_ptr<Signal>(new Signal(SignalDescription())); }, 
											[](shared_ptr<SegmentDescription> segDesc){ return; });///*/

	// ACorePolyTests::PolynomTest();
	// ACorePolyTests::ACoresTest();
	// ACorePolyTests::decomposeTest();
	// ACorePolyTests::decomposeTest2();
	// ACorePolyTests::decomposeTest1x1();

	// ACorePolyTests::PolyCoreTestPower("testPolyamp_N1k_pow3.txt", 1024, 15.0, 5, 30, 10); // todo fun : mk filename from params
	ACorePolyTests::PolyCoreTest2Cores();

	// getchar();
	
	return 0;
}

