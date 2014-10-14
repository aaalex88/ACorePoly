#pragma once

#include <vector>
#include <memory>

#include "Signal.h"
#include "PolyLibHelpers.h"
#include "Descriptions.h"


using namespace std;
using namespace std::tr1;

namespace ACorePolyLib
{
	class ISignalSolver : public Interface
	{
	public:
		virtual shared_ptr<SegmentDescription> Solve(shared_ptr<Signal> s) = 0;
	};
	
	class SignalSolver : public ISignalSolver
	{
	public:
		SignalSolver();
		~SignalSolver();

		shared_ptr<SegmentDescription> Solve(shared_ptr<Signal> s);
		static void PostProcess(shared_ptr<SegmentDescription>) {  } // вынести в параметры

	private:
		shared_ptr<SegmentDescription> m_oldDescription;

	};

	template<typename SegmentGenerator, typename ResultHandler>
	void AnalyseSignal(ISignalSolver * solver, SegmentGenerator & segGen, ResultHandler & handle)
	{
		shared_ptr<SegmentDescription> segDesc;
		while ( (segDesc = segGen()) != nullptr )
		{
			handle(solver->Solve(segDesc));
		}
	}
};