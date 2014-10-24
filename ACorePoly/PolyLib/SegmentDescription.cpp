#include <stdafx.h>

#include "SegmentDescription.h"
#include "ACore.h"



namespace ACorePolyLib 
{

	SegmentOptParams SegmentDescription::GetOptParams() const
	{
		SegmentOptParams par;
		for (int i = 0; i < cores.size(); ++i)
			par.param.push_back(cores[i].GetOptParams());
		return par;
	}

	void SegmentDescription::BuildSignal(Signal & signal) const
	{
		signal.Reset();
		for (int i = 0; i < cores.size(); ++i)
		{
			cores[i].AddToArray(signal.GetData(), signal.GetDesc().N);
		}
	}

	void SegmentDescription::Reset()
	{
		cores.clear();
	}

	bool SegmentDescription::NotEmpty() const
	{
		return cores.size() == 0;
	}

	SegmentDescription SegmentDescription::operator+ (const SegmentDescription & other) const
	{
		SegmentDescription res;
		res.desc = desc;
		res.cores = cores;
		res.cores.insert(res.cores.end(), other.cores.begin(), other.cores.end());
		return res;
	}

}