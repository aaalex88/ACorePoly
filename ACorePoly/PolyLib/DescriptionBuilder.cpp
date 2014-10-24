#include <stdafx.h>

#include "DescriptionBuilder.h"
#include "SegmentSolver.h"
#include "Optimisation.h"


namespace ACorePolyLib
{

	DescriptionBuilder::DescriptionBuilder()
	{
	}

	DescriptionBuilder::~DescriptionBuilder()
	{
	}

	void DescriptionBuilder::LinearOptimize(SegmentDescription & result, const SegmentOptParams & opt, const Signal & signal)
	{
		result.Reset();
		SegmentSolver s(opt, signal.GetDesc());
		int dim = s.GetDim();

		if (dim == 0)
			return;

		int N = signal.GetDesc().N;
		double* basis =new double(dim * N);
		s.FillBasis(N, basis);
		double * res = optimise(N, dim, basis, signal.GetData());
		result = s.ReadResult(res);

		delete[] basis;
		delete[] res;

		return;
	}

}