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

using namespace ACorePolyLib;
using namespace ACorePolyTests;

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

	// PolyCoreTestPower("testPolyamp_N4k_pow3.txt", 4*1024, 15.0, 5, 30, 10); // todo fun : mk filename from params
	PrintToFile("TestPolyAlg2Cores_N4k.txt", PolyCoreTest2Cores);
	PolyCoreTestPower("testPolyamp_N4k_pow3.txt", 4*1024, randomRange(10.0, 20.0), 3, 30, 10);
	PolyCoreTestPower("testPolyamp_N4k_pow5.txt", 4*1024, randomRange(10.0, 20.0), 5, 30, 10);

	PolyCoreTest("testPolyamp_eqpow_N4k.txt", 4*1024, randomRange(10.0, 20.0), 30, 10);
	
	return 0;
}

