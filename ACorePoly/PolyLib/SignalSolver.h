#pragma once

#include <vector>
#include <memory>

#include "Signal.h"
#include "PolyLibHelpers.h"
#include "SegmentDescription.h"


using namespace std;
using namespace std::tr1;

namespace ACorePolyLib
{
	/*
	class ISignalSolver : public Interface
	{
	public:
		virtual shared_ptr<SegmentDescription> Solve(shared_ptr<Signal> s) = 0;
	};//*/
	
	class SignalSolver // : public ISignalSolver
	{
	public:
		SignalSolver();
		~SignalSolver();

		shared_ptr<SegmentDescription> SolveSegment(shared_ptr<Signal> s);

	private:
		shared_ptr<SegmentDescription> Solve(const SegmentOptParams & startParams, const Signal & s) const;
		void ClearDesc(shared_ptr<SegmentDescription> desc) const;
		void ClearCore(ACore & core) const;
		shared_ptr<SegmentDescription> BlindFind(shared_ptr<Signal> s) const;

	private:
		shared_ptr<SegmentDescription> m_oldDescription;

	};



	template<typename SegmentGenerator, typename ResultHandler>
	void AnalyseSignal(SignalSolver * solver, SegmentGenerator & segGen, ResultHandler & handle)
	{
		shared_ptr<SegmentDescription> segDesc;
		while ( (segDesc = segGen()) != nullptr )
		{
			handle(solver->Solve(segDesc));
		}
	}

};