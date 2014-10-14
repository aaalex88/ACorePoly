#pragma once

#include "PolyLibHelpers.h"
#include "Polynom.h"


namespace ACorePolyLib
{
	
	class IAmplSolver : public Interface
	{
	public:
		virtual int GetDim() const = 0;
		virtual void FillBasis(int N, double dt, const double * phase, double * basis) const = 0;
		virtual void ReadResults(const double * res) const = 0;
	};


	class AmplSolver : public IAmplSolver
	{
	public:
		AmplSolver(int ampPow, double ph0, Polynom & resultAmpl);
		int GetDim() const;
		void FillBasis(int N, double dt, const double * phase, double * basis) const;
		void ReadResults(const double * res) const;

	private:

	};

	
	class ConstAmplSolver : public IAmplSolver
	{
	public:
		ConstAmplSolver(double & ph0, Polynom & resultAmpl);
		int GetDim() const;
		void FillBasis(int N, double dt, const double * phase, double * basis) const;
		void ReadResults(const double * res) const;

	private:

	};


};