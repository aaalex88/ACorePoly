#pragma once


#include <vector>

#include "Polynom.h"
#include "ACore.h"

using namespace std;




namespace ACorePolyLib
{

	struct ACoreParams
	{
		int freqPower;
		int ampPower;
		int maxAmpl;
	};

	struct ACoreOptParams
	{
		Polynom freq;
		vector<double> phases;
	};

	struct SegmentDescription
	{
		SignalDescription desc;
		vector<ACore> cores;
	};

	struct SegmentOptParams
	{
		vector<ACoreOptParams> param;
	};

	struct ACoreOptParamsDesc
	{
		int freqPow;
		int numPhases;
	};

	struct SegmentOptParamsDesc
	{
		int numCores;
		vector<ACoreOptParamsDesc> coreDesc;
	};

};
