#pragma once

#include "PolyLib/ACore.h"

namespace ACorePolyTests
{

	void decomposeTest();
	void decomposeTest2();
	void decomposeTest1x1();
	void ACoresTest();
	void PolyCoreTestPower(const char * reportFileName, int N, double freq, int origAmpPow, int maxNumAmpl, int maxAmpPow);
	void PolyCoreTest2Cores();
	double CompareACoreAmpls(const ACorePolyLib::ACore c1, const ACorePolyLib::ACore c2);
	void PolynomTest();
	void PolynomIntegrTest();

}