#include <stdafx.h>

#include "Descriptions.h"

namespace ACorePolyLib
{

	SegmentOptParams::SegmentOptParams(const SegmentOptParamsDesc & desc, double * data)
	{
		for (int i = 0; i < desc.coreDesc.size(); ++i)
		{
			ACoreOptParams opt(desc.coreDesc[i], data);
			param.push_back(opt);

		}
	}

	SegmentOptParams::SegmentOptParams(const SegmentOptParamsDesc & desc, vector<double> & data)
	{
		auto it = data.begin();
		for (int i = 0; i < desc.coreDesc.size(); ++i)
		{
			ACoreOptParams opt(desc.coreDesc[i], it);
			param.push_back(opt);
			it += (desc.coreDesc[i].freqPow + desc.coreDesc[i].numPhases);

		}
	}
		
	SegmentOptParamsDesc SegmentOptParams::GetDescription() const
	{
		SegmentOptParamsDesc desc;
		for (int i = 0; i < param.size(); ++i)
		{
			desc.coreDesc.push_back(param[i].GetDescription());
		}
		return desc;
	}

	vector<double> SegmentOptParams::GetOptData() const
	{
		vector<double> res;
		for (int i = 0; i < param.size(); ++i)
		{
			vector<double> addVec = param[i].GetOptData();
			res.insert(res.end(), addVec.begin(), addVec.end());
		}
		return res;
	}

	SegmentOptParams SegmentOptParams::TimeShift(double shift) const
	{
		SegmentOptParams par;
		for (int i = 0; i < param.size(); ++i)
		{
			par.param.push_back(param[i].TimeShift(shift));
		}
		return par;
	}

	SegmentOptParams & SegmentOptParams::operator+= (const SegmentOptParams & other) const
	{
	}

	ACoreOptParams::ACoreOptParams(const Polynom & _freq, const vector<double> & _phases, int _maxAmpl, int _ampPower)
		: freq(_freq)
		, phases(_phases)
		, maxAmpl(_maxAmpl)
		, ampPower(_ampPower)
	{
	}

	ACoreOptParams::ACoreOptParams(const ACoreOptParamsDesc & desc, double * data)
	{
		freq = Polynom(desc.freqPow, data);
		for (int i = 0; i < desc.numPhases; ++i)
			phases.push_back(data[desc.freqPow + i]);
		maxAmpl = desc.maxAmpl;
		ampPower = desc.ampPower;
	}

	ACoreOptParams::ACoreOptParams(const ACoreOptParamsDesc & desc, vector<double> data)
	{
		freq = Polynom(desc.freqPow, data.begin());
		phases = vector<double>(data.begin() + desc.freqPow, data.end());
		maxAmpl = desc.maxAmpl;
		ampPower = desc.ampPower;
	}

	ACoreOptParams::ACoreOptParams(const ACoreOptParams & opt)
	{
		freq = opt.freq;
		phases = opt.phases;
		maxAmpl = opt.maxAmpl;
		ampPower = opt.ampPower;
	}

	ACoreOptParamsDesc ACoreOptParams::GetDescription() const
	{
		return ACoreOptParamsDesc(freq.Power(), phases.size(), ampPower, maxAmpl);
	}

	vector<double> ACoreOptParams::GetOptData() const
	{
		vector<double> res;
		for(int i = 0; i <= freq.Power(); ++i)
			res.push_back(freq[i]);
		for(int i = 0; i < phases.size(); ++i)
			res.push_back(phases[i]);
		return res;
	}

	ACoreOptParams ACoreOptParams::TimeShift(double shift) const
	{
		double phaseShift = (freq.Integrate())(shift);
		vector<double> ph;

		for(int i = 0; i < phases.size(); ++i)
		{
			double newPh = phases[i] + (i+1) * phaseShift;
			while (newPh > pi_2)
				newPh -= pi_2;
			ph.push_back( newPh );
		}

		return ACoreOptParams(freq.TimeShift(shift), ph, maxAmpl, ampPower);
	}

	ACoreOptParamsDesc::ACoreOptParamsDesc(int _freqPow, int _numPhases, int _ampPower, int _maxAmpl)
		: freqPow(_freqPow)
		, numPhases(_numPhases)
		, ampPower(_ampPower)
		, maxAmpl(_maxAmpl)
	{
	}


}