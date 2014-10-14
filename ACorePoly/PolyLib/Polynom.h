#pragma once

#include <vector>
#include <memory>

#include "Signal.h"

using namespace std;



namespace ACorePolyLib
{

	class Polynom
	{
	private:

		static const double default_dt;

	public:
		Polynom();
		Polynom(const Polynom & poly);
		Polynom(const vector<double> & coef);
		~Polynom();

		void Reset();
		Polynom & operator=(const Polynom & poly);
		// ==,
		// +,+=,-,-=, *,*=,
		// + x, - x, ..
		// ...

		double operator[] (int ind) const;
		void SetCoef(int ind, double val);
		inline int Power() const;
		double operator() (double x) const;
		double MaxValueOnSegment(double start, double end, double dt = default_dt) const; // analitic computing should be realized
		bool IsNull() const;

		Polynom Integrate(); // defined integral with zero start value
		Polynom Differenciate();

		void FillSignal(Signal & signal) const;
		void AddToSignal(Signal & signal) const;
		// Transforms !!!
		// Serialise  !!!

	private:
		vector<double> m_coef;
		inline void DeleteTopZeros();
		inline void ExpandToIndex(int ind);
		
	};
		

};
