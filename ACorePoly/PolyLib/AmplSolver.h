#pragma once

#include "PolyLibHelpers.h"
#include "Polynom.h"


namespace ACorePolyLib
{
	
	class IAmplSolver : public Interface
	{
	public:
		virtual int GetDim() const = 0;
		virtual void FillBasis(int N, const double * phase, int k, double * basis) const = 0;
		virtual void ReadResults(const double * res) const = 0;
	};


	class AmplSolver : public IAmplSolver
	{
	public:
		AmplSolver(int ampPow, double ph0, Polynom & resultAmpl);
		~AmplSolver() {}
		int GetDim() const;
		void FillBasis(int N, const double * phase, int k, double * basis) const;
		void ReadResults(const double * res) const;

	private:
		int m_ampPow;
		int m_ph0;
		Polynom & m_result;

	};

	
	class ConstAmplSolver : public IAmplSolver
	{
	public:
		ConstAmplSolver(double & ph0, Polynom & resultAmpl);
		~ConstAmplSolver() {}
		int GetDim() const;
		void FillBasis(int N, const double * phase, int k, double * basis) const;
		void ReadResults(const double * res) const;

	private:
		double & m_resPh;
		Polynom & m_result;

	};


};