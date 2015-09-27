#include <stdafx.h>

#include "SignalSolver.h"
#include "Optimisation.h"
#include "DescriptionBuilder.h"
#include "ACore.h"


namespace ACorePolyLib
{

	SignalSolver::SignalSolver()
		: m_oldDescription(new SegmentDescription())
	{
	}

	SignalSolver::~SignalSolver()
	{
	}

	shared_ptr<SegmentDescription> SignalSolver::SolveSegment(shared_ptr<Signal> s)
	{
		shared_ptr<SegmentDescription> desc = Solve(m_oldDescription->GetOptParams().TimeShift( s->GetDesc().startTime - m_oldDescription->desc.startTime ), *s);
		shared_ptr<Signal> rest(new Signal(s->GetDesc()));
		desc->BuildSignal(*rest);
		rest->Substract(s->GetData());
		rest->Invert();
		shared_ptr<SegmentDescription> restDesc = BlindFind(rest);
		ClearDesc(restDesc);
		if (restDesc->NotEmpty())
		{
			desc = Solve(desc->GetOptParams() += restDesc->GetOptParams(), *s);
		}
		ClearDesc(desc);
		return desc;
	}

	shared_ptr<SegmentDescription> SignalSolver::Solve(const SegmentOptParams & startParams, const Signal & s) const
	{
		/*
		SegmentOptParamsDesc desc = startParams.GetDescription();
		vector<double> st_x = startParams.GetOptData();
		SegmentDescription seg;
		DescriptionBuilder builder;
		Signal tmpSig(s.GetDesc());
		
		vector<double> res_x = minimise(st_x, [&] (vector<double> x) { // or double * x as well
			SegmentOptParams opt(desc, x);
			builder.LinearOptimize(seg, opt, s);
			seg.BuildSignal(tmpSig);
			tmpSig.Substract(s.GetData());
			return tmpSig.GetNorm();

		});
		shared_ptr<SegmentDescription> res(new SegmentDescription());
		builder.LinearOptimize(*res, SegmentOptParams(desc, res_x), s);
		return res;
		// */

		throw;
		return shared_ptr<SegmentDescription>();
	}

	void SignalSolver::ClearDesc(shared_ptr<SegmentDescription> desc) const
	{
		for(size_t i = 0; i < desc->cores.size(); ++i)
		{
			ClearCore(desc->cores[i]);
		}

/*		for(int i = 0; i < desc->cores.size(); ++i)
		{
			if(desc->cores[i].GetNorm() < minCore)
			{
				desc->cores.erase(desc->cores.begin()+i, desc->cores.begin()+i);
			}
		} //&*/

	}

	void SignalSolver::ClearCore(ACore & core) const
	{
	}


	//*
	shared_ptr<SegmentDescription> SignalSolver::BlindFind(shared_ptr<Signal> s) const
	{
		shared_ptr<SegmentDescription> res(new SegmentDescription);
		res->desc = s->GetDesc();
		return res;
	}// */


}