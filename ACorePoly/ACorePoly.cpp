// ACorePoly.cpp: ���������� ����� ����� ��� ����������� ����������.
//


#include "stdafx.h"



#include "libs/FFT_lib2.cpp"
#include "libs/WAV_lib.cpp"
#include "libs/ACore_lib.cpp"
#include "libs/ACore_lib2.cpp"
#include "libs/ACore_lib3.cpp"


#include "PolyLib/ACorePolyLib.h"
#include "Tests/Tests.h"




int _tmain(int argc, char* argv[])
{
	/*
	Solution<SignalSolver>::AnalyseSignal<shared_ptr<Signal>,shared_ptr<SegmentDescription>>
										( [](){ return shared_ptr<Signal>(new Signal(SignalDescription())); }, 
											[](shared_ptr<SegmentDescription> segDesc){ return; });///*/
	ACorePolyTests::ACoresTest();

	getchar();
	
	return 0;
}

