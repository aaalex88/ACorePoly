#pragma once

#include "Descriptions.h"

namespace ACorePolyLib 
{

	struct SegmentDescription
	{
		SignalDescription desc;
		vector<ACore> cores;

		SegmentOptParams GetOptParams() const;
		void BuildSignal(Signal & signal) const;

		bool NotEmpty() const;
		SegmentDescription operator+ (const SegmentDescription & other) const;

		void Reset();

	};


}