

#include "PolyLibHelpers.h"
#include <math.h>

namespace ACorePolyLib
{
	
	bool AlmostZero(double d)
	{
		return abs(d) < eps;
	}

	double ComputePhase(double cosAmp, double sinAmp)
	{
		double d = sqrt(cosAmp*cosAmp + sinAmp*sinAmp);
		if (AlmostZero(d))
			return 0;
		cosAmp /= d;
		sinAmp /= d;
		double ph = acos(cosAmp);
		if (sinAmp < 0)
			ph = pi_2 - ph;
		return ph;
	}

};