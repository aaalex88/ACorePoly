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
		double z[6] = { 0, 1, 0, 2.0 / 3.0, 0, 0.2 };
		Polynom p(3, x);
		Polynom psqr(5, y);
		Polynom psqrint(6, z);
		printf("\npolys #1 : %s\n", p.Sqr() == psqr ? "true" : "false");
		printf("\npolys #2 : %s\n", psqr.Integrate() == psqrint ? "true" : "false");
	}
}