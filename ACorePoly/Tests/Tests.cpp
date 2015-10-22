#include "stdafx.h"

#include <algorithm>
#include <tuple>

#include "Tests.h"
#include "PolyLib/Optimisation.h"
#include "PolyLib/ACore.h"
#include "Polylib/PolyLibHelpers.h"
#include "Polylib/SegmentDescription.h"
#include "Polylib/DescriptionBuilder.h"

using namespace ACorePolyLib;

namespace ACorePolyTests
{

	
	// Not need because of next fun
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
		printf("\nDifference between cores amplitudes is : %f \n", get<0>(CompareACoreAmpls(core, res.cores[0])));

		ACorePolyLib::Signal sRes(desc.desc);
		res.cores[0].AddToArray(sRes.GetData(), sRes.GetDesc().N);

		printf("\nOriginal signal norm is : %f\n", s.GetNorm());
		printf("\nResult signal norm is : %f\n", sRes.GetNorm());
		sRes.Substract(s.GetData());
		printf("\nDifference between signals is : %f \n", sRes.GetNorm());

	}


	// TODO : split calculation and printing to diff funs?
	void PolyCoreTestPower(const char * reportFileName, int N, double freq, 
		int origAmpPow, int maxNumAmpl, int maxAmpPow)
	{
		FILE * report = nullptr; // TODO : use PrintToFile & FILE *
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


				double coresDiff;
				double maxAmpDiff;
				tie(coresDiff, maxAmpDiff) = CompareACoreAmpls(core, res.cores[0]);

				fprintf_s(report, "Ampl power of result is %d", ampPow);
				fprintf_s(report, "\nRatio of difference between cores is : %f", coresDiff / core.L2Norm());
				fprintf_s(report, "\nRatio of worst difference between cores amplitudes is : %f", maxAmpDiff / core.L2Norm());
				sRes.Substract(s.GetData());
				fprintf_s(report, "\nRatio of norm of difference to orig signal is : %f", sRes.GetNorm() / s.GetNorm()); // TODO : better english!

				fprintf_s(report, "\n\n");
			}
			
			fprintf(report, "\n\n");
		}

		
		fclose(report);
	}


	void PolyCoreTest(const char * reportFileName, int N, double freq,
		int maxNumAmpl, int maxAmpPow) // TODO : big refactor, code duplicate
	{
		FILE * report = nullptr; // TODO : use PrintToFile & FILE *
		errno_t err = fopen_s(&report, reportFileName, "w");
		if (!report)
			throw new runtime_error("couldn't open file in PolyCoreTest"); // TODO : handle
		fprintf_s(report,
			"PolyCore Linear Optimization full test\nSignal length is %d; approximate frequency is %f\n\n",
			N, freq);


		for (int numAmpl = 1; numAmpl <= maxNumAmpl; ++numAmpl) {

			fprintf_s(report, "Original Signal: number of amplitudes is %d\n", numAmpl);

			for (int ampPow = 0; ampPow <= maxAmpPow; ++ampPow) {

				ACorePolyLib::ACore core; // = RandomCore(...)
				ACorePolyLib::RandomCore(core, freq, numAmpl, ampPow);
				

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



				optParams.param[0].ampPower = ampPow;
				builder.LinearOptimize(res, optParams, s);

				sRes.Reset();
				res.cores[0].AddToArray(sRes.GetData(), sRes.GetDesc().N);


				double coresDiff;
				double maxAmpDiff;
				tie(coresDiff, maxAmpDiff) = CompareACoreAmpls(core, res.cores[0]);

				fprintf_s(report, "Ampl power of signal is %d", ampPow);
				fprintf_s(report, "\nRatio of difference between cores is : %f", coresDiff / core.L2Norm());
				fprintf_s(report, "\nRatio of worst difference between cores amplitudes is : %f", maxAmpDiff / core.L2Norm());
				sRes.Substract(s.GetData());
				fprintf_s(report, "\nRatio of norm of difference to orig signal is : %f", sRes.GetNorm() / s.GetNorm()); // TODO : better english!

				fprintf_s(report, "\n\n");
			}

			fprintf(report, "\n\n");
		}


		fclose(report);
	}


	void PolyCoreTest2Cores(FILE * file)
	{
		int ampPow = 5;
		int N = 4*1024; // TODO : parametr

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


		double diff1 = get<0>(CompareACoreAmpls(core1, res.cores[0])) / core1.L2Norm();
		double diff2 = get<0>(CompareACoreAmpls(core2, res.cores[1])) / core2.L2Norm();
		fprintf(file, "\nRatio of difference between cores amplitudes are : %f, %f", diff1, diff2);
		fprintf(file, "\nRatio of norm of difference to orig signal is : %f\n", sDiff.GetNorm() / s.GetNorm());

	}

	


}