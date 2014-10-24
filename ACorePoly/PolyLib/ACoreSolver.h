#pragma once

#include <memory>

//#include "ILinearSolver.h"
#include "AmplSolver.h"
#include "ACore.h"
#include "Descriptions.h"


namespace ACorePolyLib
{

	class ACoreSolver // : public ILinearSolver<ACore>
	{
	public:
		ACoreSolver(ACore & core, const ACoreOptParams & optParams, const SignalDescription & desc);
		~ACoreSolver();

		int GetDim() const;
		void FillBasis(int N, double * basis) const;
		void ReadResults(const double * res);

	private:
		vector<shared_ptr<IAmplSolver> > m_amplSolvers;
		ACore & m_core;
		SignalDescription m_desc;

	};
};