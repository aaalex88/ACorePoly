#include <stdafx.h>

#include "ACoreSolver.h"
#include <math.h>
#include <algorithm>

namespace ACorePolyLib
{

	ACoreSolver::ACoreSolver(ACore & core, const ACoreOptParams & optParams, const SignalDescription & desc)
		: m_core(core)
		, m_desc(desc)
	{
		int i;

		// TODO: refactoring
		m_core.SetOptParams(optParams);

		double maxFr = 1.0 / m_desc.dt * 2;
		int maxPossibleAmpl =  (int)(maxFr / optParams.freq.MaxValueOnSegment(0, 1, 0.01));
		int numAmpl = (std::min)(maxPossibleAmpl, optParams.maxAmpl);
		int numExistingAmpl = (std::min)((int)m_core.ampl.size(), numAmpl);

		m_core.ampl.resize(numAmpl);
		m_core.phases.resize(numAmpl, 0);


		for (i = 0 ; i < numExistingAmpl; ++i)
		{
			m_amplSolvers.push_back( shared_ptr<IAmplSolver>(new AmplSolver(optParams.ampPower, m_core.phases[i], m_core.ampl[i])) );
		}
		for ( ; i < numAmpl; ++i)
		{
			m_amplSolvers.push_back( shared_ptr<IAmplSolver>(new ConstAmplSolver(m_core.phases[i], m_core.ampl[i])) );
		}

	}

	ACoreSolver::~ACoreSolver()
	{
	}

	int ACoreSolver::GetDim() const
	{
		int res = 0;
		for (size_t i = 0; i < m_amplSolvers.size(); ++i)
		{
			res += m_amplSolvers[i]->GetDim();
		}

		return res;
	}

	void ACoreSolver::FillBasis(int N, double * basis) const
	{
		Signal phase(SignalDescription(N, 1.0 / double(N), 0));
		m_core.freq.Integrate().FillSignal(phase);

		for (size_t i = 0; i < m_amplSolvers.size(); ++i)
		{
			m_amplSolvers[i]->FillBasis(N, phase.GetData(), i+1, basis);
			basis += N * m_amplSolvers[i]->GetDim();
		}
	}

	void ACoreSolver::ReadResults(const double * res)
	{
		for (size_t i = 0; i < m_amplSolvers.size(); ++i)
		{
			m_amplSolvers[i]->ReadResults(res);
			res += m_amplSolvers[i]->GetDim();
		}
	}

};