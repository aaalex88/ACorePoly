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

	bool DescriptionBuilder::LinearOptimize(SegmentDescription & result, const SegmentOptParams & opt, const Signal & signal)
	{
		result.Reset();
		SegmentSolver s(opt, signal.GetDesc());
		int dim = s.GetDim();

		if (dim == 0)
			return false;

		int N = signal.GetDesc().N;
		double* basis =new double[dim * N];
		s.FillBasis(N, basis);

		double * res = new double[dim];

		bool retStatus = decompose(N, dim, basis, signal.GetData(), res);

		result = s.ReadResult(res);
		delete[] basis;
		delete[] res;

		return retStatus;
	}

}