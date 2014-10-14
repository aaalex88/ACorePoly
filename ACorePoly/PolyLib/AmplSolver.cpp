

#include "AmplSolver.h"
#include <math.h>


namespace ACorePolyLib
{

	AmplSolver::AmplSolver(int ampPow, double ph0, Polynom & resultAmpl)
		: m_ampPow(ampPow)
		, m_ph0(ph0)
		, m_result(resultAmpl)
	{
		assert(m_ampPow >= 0);
	}

	int AmplSolver::GetDim() const
	{
		return m_ampPow + 1;
	}

	void AmplSolver::FillBasis(int N, double dt, const double * phase, int k, double * basis) const
	{
		int i;
		for (i = 0; i < N; ++i)
		{
			basis[i] = cos(m_ph0 + k * phase[i]);
		}
		for (int m = 1; m <= m_ampPow; ++m)
		{
			for (i = 0; i < N; ++i)
			{
				basis[N*m + i] = pow(double(i)/double(N), m) * basis[i];
			}
		}
	}

	void AmplSolver::ReadResults(const double * res) const
	{
		m_result.Reset();
		for (int i = 0; i <= m_ampPow; ++i)
		{
			m_result.SetCoef(i, res[i]);
		}
	}



	ConstAmplSolver::ConstAmplSolver(double & ph0, Polynom & resultAmpl)
		: m_resPh(ph0)
		, m_result(resultAmpl)
	{
	}

	int ConstAmplSolver::GetDim() const
	{
		return 2;
	}

	void ConstAmplSolver::FillBasis(int N, double dt, const double * phase, int k, double * basis) const
	{
		for (int i = 0; i < N; ++i)
		{
			basis[i] = cos(k * phase[i]);
			basis[N + i] = sin(k * phase[i]);
		}
	}

	void ConstAmplSolver::ReadResults(const double * res) const
	{
		double c = res[0];
		double s = res[1];
		m_resPh = ComputePhase(c, s);
		m_result.Reset();
		m_result.SetCoef(0, sqrt(c*c + s*s));
	}

};