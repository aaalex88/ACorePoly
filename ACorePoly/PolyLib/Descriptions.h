#pragma once


#include <vector>

#include "Polynom.h"


using namespace std;




namespace ACorePolyLib
{
	struct ACore;

	struct ACoreOptParamsDesc;

	struct ACoreOptParams
	{
		Polynom freq;
		vector<double> phases;
		int ampPower; // not in opt params, should be set in acore solver!
		int maxAmpl;

		ACoreOptParams(const Polynom & _freq, const vector<double> & _phases, int _maxAmpl, int _ampPower);
		ACoreOptParams(const ACoreOptParamsDesc & desc, double * data);
		ACoreOptParams(const ACoreOptParamsDesc & desc, vector<double> data);
		template<typename Iter>
		ACoreOptParams(const ACoreOptParamsDesc & desc, Iter it);
		ACoreOptParams(const ACoreOptParams & opt);
		
		ACoreOptParamsDesc GetDescription() const;
		vector<double> GetOptData() const;

		ACoreOptParams TimeShift(double shift) const;
	};

	struct ACoreOptParamsDesc
	{
		int freqPow;
		int numPhases;
		int ampPower;
		int maxAmpl;

		ACoreOptParamsDesc(int _freqPow, int _numPhases, int _ampPower, int _maxAmpl);
	};



	struct SegmentOptParamsDesc;

	struct SegmentOptParams
	{
		vector<ACoreOptParams> param;
		
		SegmentOptParams() {}
		SegmentOptParams(const SegmentOptParamsDesc & desc, double * data);
		SegmentOptParams(const SegmentOptParamsDesc & desc, vector<double> & data);
				
		SegmentOptParamsDesc GetDescription() const;
		vector<double> GetOptData() const;
		SegmentOptParams TimeShift(double shift) const;

		SegmentOptParams & operator+= (const SegmentOptParams & other);
	};
	
	struct SegmentOptParamsDesc
	{
		vector<ACoreOptParamsDesc> coreDesc;
	};
	


	template<typename Iter>
	ACoreOptParams::ACoreOptParams(const ACoreOptParamsDesc & desc, Iter it)
	{
		freq = Polynom(desc.freqPow, it);
		it += desc.freqPow;
		for (int i = 0; i < desc.numPhases; ++i)
		{
			phases.push_back(*it);
			it++;
		}
		maxAmpl = desc.maxAmpl;
		ampPower = desc.ampPower;
	}
};
