#pragma once


#include <vector>

#include "Polynom.h"

using namespace std;




namespace ACorePolyLib
{

	struct ACore
	{
		Polynom freq;
		int maxAmpl;
		vector<double> phases;
		vector<Polynom> ampl;

		ACore():maxAmpl(0) {}
		ACore(const Polynom & _freq, int _maxAmpl): freq(_freq), maxAmpl(_maxAmpl) {}

		// GetOptParams
		// GetOptParamsDesc
	};

};
