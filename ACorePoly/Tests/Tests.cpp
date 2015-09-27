#include "stdafx.h"

#include "Tests.h"
#include "PolyLib/Optimisation.h"


namespace ACorePolyTests
{

	void decomposeTest()
	{
		int N = 5;
		int dim = 3;
		double basis[15];
		double res[3];
		double signal[5] = {1., 2., 3., 4., 5.};
		memset((void*)basis, 0, 15 * sizeof(double));
		for (int i = 0; i < 3; ++i) {
			basis[6 * i] = 1;
		}
		ACorePolyLib::decompose(N, dim, basis, signal, res);

		for (int i = 0; i < 3; ++i) {
			printf("%f ", res[i]);
		}

		// OutputDebugString();
	}

}