#include <stdafx.h>


#include "SegmentSolver.h"

namespace ACorePolyLib
{

	SegmentSolver::SegmentSolver(const SegmentOptParams & params, const SignalDescription & desc)
	{
		m_desc.desc = desc;
		m_desc.cores.resize(params.param.size());
		for (int i = 0; i < params.param.size(); ++i)
		{
			m_solvers.push_back( shared_ptr<ACoreSolver>( new ACoreSolver(m_desc.cores[i], params.param[i], desc) ) );
		}
	}

	SegmentSolver::~SegmentSolver()
	{
	}

	int SegmentSolver::GetDim() const
	{
		int res = 0;
		for (int i = 0; i < m_solvers.size(); ++i)
		{
			res += m_solvers[i]->GetDim();
		}

		return res;
	}

	void SegmentSolver::FillBasis(int N, double * basis) const
	{
		for (int i = 0; i < m_solvers.size(); ++i)
		{
			m_solvers[i]->FillBasis(N, basis);
			basis += N * m_solvers[i]->GetDim();
		}
	}

	SegmentDescription & SegmentSolver::ReadResult(double * res)
	{
		for (int i = 0; i < m_solvers.size(); ++i)
		{
			m_solvers[i]->ReadResults(res);
			res += m_solvers[i]->GetDim();
		}

		return m_desc;
	}


}

