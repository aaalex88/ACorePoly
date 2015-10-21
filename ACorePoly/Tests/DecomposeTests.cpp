#include "stdafx.h"

#include <algorithm>

#include "Tests.h"
#include "PolyLib/Optimisation.h"
#include "PolyLib/ACore.h"
#include "Polylib/PolyLibHelpers.h"
#include "Polylib/SegmentDescription.h"
#include "Polylib/DescriptionBuilder.h"


namespace ACorePolyTests
{
	void decomposeTest()
	{
		int N = 5;
		int dim = 3;
		double basis[15];
		double res[3];
		double signal[5] = { 1., 2., 3., 4., 5. };
		memset((void*)basis, 0, 15 * sizeof(double));
		for (int i = 0; i < 3; ++i) {
			basis[6 * i] = 1;
		}
		ACorePolyLib::decompose(N, dim, basis, signal, res);

		for (int i = 0; i < 3; ++i) {
			printf("%f ", res[i]);
		}
		printf("\n");

		// OutputDebugString();
	}

	void decomposeTest2()
	{
		const int N = 5;
		const int dim = 2;
		double basis[dim * N];
		double res[dim];
		double signal[N];
		memset((void*)basis, 0, N * sizeof(double));
		for (int i = 0; i < N; ++i) {
			signal[i] = 10.0;
			basis[i] = 1.0;
			basis[N + i] = 0.0;
		}
		basis[N] = 1.0;

		ACorePolyLib::decompose(N, dim, basis, signal, res);

		for (int i = 0; i < dim; ++i) {
			printf("%f ", res[i]);
		}
		printf("\n");
	}

	void decomposeTest1x1()
	{
		const int N = 4;
		const int dim = 1;
		double basis[dim * N];
		double res[dim];
		double signal[N];
		memset((void*)basis, 0, N * sizeof(double));
		for (int i = 0; i < N; ++i) {
			signal[i] = 10.0;
			basis[i] = 1.0;
		}

		ACorePolyLib::decomposeAlt(N, dim, basis, signal, res);

		for (int i = 0; i < dim; ++i) {
			printf("%f ", res[i]);
		}
		printf("\n");
	}
}