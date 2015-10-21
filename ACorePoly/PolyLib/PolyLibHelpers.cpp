#include <stdafx.h>
#include <math.h>
#include <algorithm>
#include <random>

#include "PolyLibHelpers.h"
#include "Polynom.h"
#include "ACore.h"


namespace ACorePolyLib
{
	
	bool AlmostZero(double d)
	{
		return fabs(d) < eps;
	}

	/*
	int max(int x, int y)
	{
		return x > y ? x : y;
	}

	int min(int x, int y)
	{
		return x < y ? x : y;
	} // */

	double ComputePhase(double cosAmp, double sinAmp)
	{
		double d = sqrt(cosAmp*cosAmp + sinAmp*sinAmp);
		if (AlmostZero(d))
		{
			assert(!"ComputePhase : Zero complex amplitude, cant compute correct phase!");
			return 0;
		}
		cosAmp /= d;
		sinAmp /= d;
		double ph = acos(cosAmp);
		if (sinAmp < 0)
			ph = pi_2 - ph;
		return ph;
	}

	double Factorial(int n)
	{
		double res = 1;
		for(int i = 2; i <= n; ++i)
		{
			res *= i;
		}
		return res;
	}

	double random(double max)
	{
		// TODO : adequate realization of rand
		double d = double(rand()) / double(RAND_MAX);
		return d * max;

		/*
		std::default_random_engine eng{};
		static std::uniform_real_distribution<>	d{};
		using parm_t = decltype(d)::param_type;
		return d(eng(), parm_t{0.0, max});*/
	}

	double randomRange(double start, double end)
	{
		assert(start <= end && "Incorrect interval!");
		return random(end - start) + start;
	}

	void RandomPolynomLowVary(Polynom & res, int pow, double startVal, double delta)
	{
		assert(pow >= 0 && "function RandomPolynom : Polynom cant have power below zero!");

		res.Reset();
		res.SetCoef(0, startVal);
		for (int i = 1; i <= pow; ++i) {
			res.SetCoef(i, randomRange(-delta/double(i), delta/double(i)));
		}

	}

	void RandomPolynom(Polynom & res, int pow, double startVal, double delta)
	{
		assert(pow >= 0 && "function RandomPolynom : Polynom cant have power below zero!");

		res.Reset();
		res.SetCoef(0, startVal);
		for (int i = 1; i <= pow; ++i) {
			res.SetCoef(i, randomRange(-delta, delta));
		}

	}

	//void RandomPolynomByTimeValues(Polynom & res, int pow, vector<double> vals)
	//	// vals.size() == pow + 1 !
	//{
	//	//
	//}

	//void RandomPolynomByTimeValues(Polynom & res, int pow, double * vals)
	//	// size of vals == pow + 1 !
	//{

	//}

	void RandomCoreLowVary(ACore & res, double fr, int numAmpl, int ampPower)
	{
		assert(fr > 0.0 && numAmpl >= 0 && ampPower >= 0);

		RandomPolynomLowVary(res.freq, 2, fr, 0.05 * fr);
		res.maxAmpl = numAmpl; 
		res.ampPower = ampPower;
		res.ampl.clear();
		res.phases.clear();
		Polynom ampl;
		for (int i = 0; i < numAmpl; ++i) {
			res.phases.push_back(random(pi_2));
			RandomPolynomLowVary(ampl, ampPower, 10.0/double(i+1), 3.0);
			res.ampl.push_back(ampl);
		}
	}

	void RandomCore(ACore & res, double fr, int numAmpl, int ampPower)
	{
		assert(fr > 0.0 && numAmpl >= 0 && ampPower >= 0);

		RandomPolynom(res.freq, 2, fr, 0.05 * fr);
		res.maxAmpl = numAmpl;
		res.ampPower = ampPower;
		res.ampl.clear();
		res.phases.clear();
		Polynom ampl;
		for (int i = 0; i < numAmpl; ++i) {
			res.phases.push_back(random(pi_2));
			RandomPolynom(ampl, ampPower, 10.0, 3.0);
			res.ampl.push_back(ampl);
		}
	}

	double CompareACoreAmpls(const ACorePolyLib::ACore c1, const ACorePolyLib::ACore c2)
	{
		double diff = 0;

		int minCores = (std::min)(c1.ampl.size(), c2.ampl.size());
		int maxCores = (std::max)(c1.ampl.size(), c2.ampl.size());
		const ACorePolyLib::ACore & bigCore = (c1.ampl.size() == maxCores) ? c1 : c2;

		for (int i = 0; i < minCores; ++i) {
			diff += sqr((c1.ampl[i] - c2.ampl[i]).L2Norm());
			// L2-norm, TODO : make Polynom Function
		}

		for (int i = minCores; i < maxCores; ++i) {
			diff += sqr(bigCore.ampl[i].L2Norm());
		}

		return sqrt(diff);
	}

};