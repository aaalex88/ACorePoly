#pragma once


#include <vector>

#include "Polynom.h"
#include "Descriptions.h"

using namespace std;



// 1
namespace ACorePolyLib
{

	struct ACore
	{
		Polynom freq;
		vector<double> phases;
		vector<Polynom> ampl;
		int maxAmpl;
		int ampPower;

		ACore():maxAmpl(0) {}
		ACore(const Polynom & _freq, int _maxAmpl): freq(_freq), maxAmpl(_maxAmpl) {}

		
		ACoreOptParams GetOptParams() const;
		void SetOptParams(const ACoreOptParams & optParams);

		void AddToArray(double * arr, int N) const;
		double L2Norm() const;

	};

};
