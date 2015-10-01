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

	void ACoresTest()
	{
		ACorePolyLib::ACore core;
		ACorePolyLib::RandomCore(core, 10, 10, 3);
		ACorePolyLib::SegmentDescription desc;
		desc.cores.push_back(core);
		desc.desc.N = 1024;
		desc.desc.dt = 1.0 / 44100.0; // TODO : const!
		desc.desc.startTime = 0.0;


		ACorePolyLib::Signal s(desc.desc);
		core.AddToArray(s.GetData(), s.GetDesc().N);

		ACorePolyLib::SegmentDescription res;
		ACorePolyLib::DescriptionBuilder builder;
		builder.LinearOptimize(res, desc.GetOptParams(), s);

		// CompareSegmentDescriptions(res, desc);
		CompareACoreAmpls(core, res.cores[0]);
	}

	void CompareACoreAmpls(ACorePolyLib::ACore c1, ACorePolyLib::ACore c2)
	{
		double diff = 0;

		int minCores = (std::min)(c1.ampl.size(), c2.ampl.size());
		int maxCores = (std::max)(c1.ampl.size(), c2.ampl.size());
		ACorePolyLib::ACore & bigCore = (c1.ampl.size() == maxCores) ? c1 : c2;

		for (int i = 0; i < minCores; ++i) {
			diff += (c1.ampl[i] - c2.ampl[i]).Integrate()(1.0); // L1-norm? no: is correct only under assumption that polynom is positive!
			// sqrt( ().Sqr().Integrate() ) => L2-norm
		}

		for (int i = minCores; i < maxCores; ++i) {
			diff += bigCore.ampl[i].Integrate()(1.0);
		}
		printf("\nDifference between cores amplitudes is : %f \n", diff);

	}

}