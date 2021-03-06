#pragma once

#include <memory>
#include "SegmentDescription.h"
#include "Signal.h"


namespace ACorePolyLib
{

	class DescriptionBuilder
	{
	public:
		DescriptionBuilder();
		~DescriptionBuilder();

		bool LinearOptimize(SegmentDescription & result, const SegmentOptParams & opt, const Signal & signal);

	private:

	};

}

