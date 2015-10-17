#pragma once

#include "Descriptions.h"

namespace ACorePolyLib 
{

	struct SegmentDescription
	{
		SignalDescription desc; // TODO : there we dont need N/dt; only start/end time. So double st/en will be enough
		vector<ACore> cores;

		SegmentOptParams GetOptParams() const;
		void BuildSignal(Signal & signal) const;

		bool NotEmpty() const;
		SegmentDescription operator+ (const SegmentDescription & other) const;

		void Reset();

	};


}