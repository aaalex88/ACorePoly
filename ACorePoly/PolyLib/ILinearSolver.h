#pragma once

namespace ACorePolyLib
{
	
	template<typename T>
	class ILinearSolver : public Interface
	{
	public:
		virtual int GetDim() const = 0;
		virtual void FillBasis(int N, double dt, double * basis) const = 0;
		virtual T ReadResults(const double * res) const = 0;
	};
};