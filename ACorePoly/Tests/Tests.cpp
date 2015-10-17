#include "stdafx.h"

#include <algorithm>

#include "Tests.h"
#include "PolyLib/Optimisation.h"
#include "PolyLib/ACore.h"
#include "Polylib/PolyLibHelpers.h"
#include "Polylib/SegmentDescription.h"
#include "Polylib/DescriptionBuilder.h"

using namespace ACorePolyLib;

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

	void ACoresTest()
	{
		ACorePolyLib::ACore core;
		ACorePolyLib::RandomCore(core, 18.24, 25, 7);
		ACorePolyLib::SegmentDescription desc;
		desc.cores.push_back(core);
		desc.desc.N = 1 * 1024;
		desc.desc.dt = stdDeltaTime;
		desc.desc.startTime = 0.0;


		ACorePolyLib::Signal s(desc.desc);
		core.AddToArray(s.GetData(), s.GetDesc().N);

		ACorePolyLib::SegmentDescription res;
		ACorePolyLib::DescriptionBuilder builder;
		builder.LinearOptimize(res, desc.GetOptParams(), s);

		
		// CompareSegmentDescriptions(res, desc);
		printf("\nDifference between cores amplitudes is : %f \n", CompareACoreAmpls(core, res.cores[0]));

		ACorePolyLib::Signal sRes(desc.desc);
		res.cores[0].AddToArray(sRes.GetData(), sRes.GetDesc().N);

		printf("\nOriginal signal norm is : %f\n", s.GetNorm());
		printf("\nResult signal norm is : %f\n", sRes.GetNorm());
		sRes.Substract(s.GetData());
		printf("\nDifference between signals is : %f \n", sRes.GetNorm());

	}

	void PolyCoreTestPower(const char * reportFileName, int N, double freq, 
		int origAmpPow, int maxNumAmpl, int maxAmpPow)
	{
		FILE * report = nullptr;
		errno_t err = fopen_s(&report, reportFileName, "w");
		if (!report)
			throw; // TODO : handle
		fprintf_s(report, 
			"PolyCore Linear Optimization full test\nSignal length is %d; power of amplitudes is %d; approximate frequency is %f\n\n", 
			N, origAmpPow, freq);


		for (int numAmpl = 1; numAmpl <= maxNumAmpl; ++numAmpl) {

			ACorePolyLib::ACore core; // = RandomCore(...)
			ACorePolyLib::RandomCore(core, freq, numAmpl, origAmpPow);
			fprintf_s(report, "Original Signal: number of amplitudes is %d, Ampl power is %d\n", numAmpl, origAmpPow);

			ACorePolyLib::SegmentDescription desc; // TODO : SegDesc(core, stdSignalDescription) / or vector<ACore>
			desc.cores.push_back(core);
			desc.desc.N = N;		
			desc.desc.dt = stdDeltaTime;
			desc.desc.startTime = 0.0;

			ACorePolyLib::Signal s(desc.desc); // TODO : Signal(ACore) / Signal(SegDesc)
			core.AddToArray(s.GetData(), s.GetDesc().N);

			ACorePolyLib::SegmentDescription res;
			ACorePolyLib::DescriptionBuilder builder;
			ACorePolyLib::Signal sRes(desc.desc);

			SegmentOptParams optParams = desc.GetOptParams(); // ? auto copy constructor

			for (int ampPow = 0; ampPow <= maxAmpPow; ++ampPow) {
				optParams.param[0].ampPower = ampPow;
				builder.LinearOptimize(res, optParams, s);

				sRes.Reset();
				res.cores[0].AddToArray(sRes.GetData(), sRes.GetDesc().N);

				fprintf_s(report, "Ampl power of result is %d", ampPow);
				fprintf_s(report, "\nRatio of difference between cores amplitudes is : %f", CompareACoreAmpls(core, res.cores[0]) / core.L2Norm());
				sRes.Substract(s.GetData());
				fprintf_s(report, "\nRatio of norm of difference to orig signal is : %f", sRes.GetNorm() / s.GetNorm()); // TODO : better english!

				fprintf_s(report, "\n");
			}
			
			fprintf(report, "\n\n");
		}

		
		fclose(report);
	}

	void PolyCoreTest2Cores()
	{
		int ampPow = 5;
		int N = 1024;

		ACorePolyLib::ACore core1;
		ACorePolyLib::RandomCore(core1, 21, 4, ampPow);
		ACorePolyLib::ACore core2;
		ACorePolyLib::RandomCore(core2, 15.2, 6, ampPow);

		ACorePolyLib::SegmentDescription desc; // TODO
		desc.cores.push_back(core1);
		desc.cores.push_back(core2);
		desc.desc.N = N;
		desc.desc.dt = stdDeltaTime;
		desc.desc.startTime = 0.0;


		ACorePolyLib::Signal s(desc.desc);
		desc.BuildSignal(s);

		ACorePolyLib::SegmentDescription res;
		ACorePolyLib::DescriptionBuilder builder;
		builder.LinearOptimize(res, desc.GetOptParams(), s);


		ACorePolyLib::Signal sDiff(desc.desc);
		res.BuildSignal(sDiff);
		sDiff.Substract(s.GetData());

		double diff1 = CompareACoreAmpls(core1, res.cores[0]) / core1.L2Norm();
		double diff2 = CompareACoreAmpls(core2, res.cores[1]) / core2.L2Norm();

		printf("\nRatios of difference between cores amplitudes are : %f, %f", diff1, diff2);
		printf("\nRatio of norm of difference to orig signal is : %f", sDiff.GetNorm() / s.GetNorm());

		getchar();
	}

	double CompareACoreAmpls(const ACorePolyLib::ACore c1, const ACorePolyLib::ACore c2)
	{
		double diff = 0;

		int minCores = (std::min)(c1.ampl.size(), c2.ampl.size());
		int maxCores = (std::max)(c1.ampl.size(), c2.ampl.size());
		const ACorePolyLib::ACore & bigCore = (c1.ampl.size() == maxCores) ? c1 : c2;

		for (int i = 0; i < minCores; ++i) {
			diff += sqr( (c1.ampl[i] - c2.ampl[i]).L2Norm() );
			// L2-norm, TODO : make Polynom Function
		}

		for (int i = minCores; i < maxCores; ++i) {
			diff += sqr( bigCore.ampl[i].L2Norm() );
		}
		
		return sqrt(diff);
	}


	void PolynomIntegrTest()
	{
		ACorePolyLib::Signal sig(ACorePolyLib::SignalDescription(1024, fdiv(1, 1024), 0.0));
		Polynom p(1.0);
		Polynom q(p.Integrate());
		p.FillSignal(sig);
		printf("\nSignal norm is : %f\n", sig.GetNorm());
		p.Integrate().FillSignal(sig);
		printf("\nSignal norm is : %f\n", sig.GetNorm());
	}

	void PolynomTest()
	{
		double x[3] = { 1, 0, 1 };
		double y[5] = { 1, 0, 2, 0, 1 };
		double z[6] = { 0, 1, 0, 2.0/3.0, 0, 0.2 };
		Polynom p(3, x);
		Polynom psqr(5, y);
		Polynom psqrint(6, z);
		printf("\npolys #1 : %s\n", p.Sqr() == psqr ? "true" : "false");
		printf("\npolys #2 : %s\n", psqr.Integrate() == psqrint ? "true" : "false");
	}

}