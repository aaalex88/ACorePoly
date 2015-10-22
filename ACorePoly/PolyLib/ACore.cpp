#include <stdafx.h>

#include "ACore.h"


namespace ACorePolyLib
{

	ACoreOptParams ACore::GetOptParams() const
	{
		return ACoreOptParams(freq, phases, maxAmpl, ampPower);
	}

	void ACore::SetOptParams(const ACoreOptParams & optParams)
	{
		freq = optParams.freq;
		phases = optParams.phases;
		maxAmpl = optParams.maxAmpl;
		ampPower = optParams.ampPower;
	}

	void ACore::AddToArray(double * arr, int N) const
		// assume arr have the same time interval, as ACore was built from
	{
		Polynom phase = freq.Integrate();
		double dt = 1.0 / double(N);
		for (int i = 0; i < N; ++i) // TODO : change cycles order?
		{
			double ph = phase(i * dt);
			for (size_t k = 0; k < ampl.size(); ++k)
			{
				arr[i] += ampl[k](i * dt) * cos(phases[k] + (k+1)*ph);
			}
		}
	}

	double ACore::L2Norm() const
	{
		double res = 0;
		for (size_t i = 0; i < ampl.size(); ++i) {
			res += sqr( ampl[i].L2Norm() );
		}

		return sqrt(res);
	}

}