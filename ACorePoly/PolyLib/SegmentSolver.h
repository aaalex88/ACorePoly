#pragma once

//#include "ILinearSolver.h"
#include "ACore.h"
#include "ACoreSolver.h"
#include "SegmentDescription.h"


namespace ACorePolyLib
{

	class SegmentSolver // : public ILinearSolver<SegmentDescription>
	{
	public:
		SegmentSolver(const SegmentOptParams & params, const SignalDescription & desc);
		~SegmentSolver();

		int GetDim() const;
		void FillBasis(int N, double * basis) const;
		SegmentDescription & ReadResult(double * res);

	private:
		vector<shared_ptr<ACoreSolver>> m_solvers;
		SegmentDescription  m_desc;
	};
};