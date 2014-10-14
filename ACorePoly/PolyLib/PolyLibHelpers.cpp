

#include "PolyLibHelpers.h"
#include <math.h>

namespace ACorePolyLib
{
	const double eps = 10e-15;
	
	bool AlmostZero(double d)
	{
		return abs(d) < eps;
	}

};