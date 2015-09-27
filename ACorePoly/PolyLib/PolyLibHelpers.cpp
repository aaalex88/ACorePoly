#include <stdafx.h>

#include "PolyLibHelpers.h"
#include <math.h>

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

};