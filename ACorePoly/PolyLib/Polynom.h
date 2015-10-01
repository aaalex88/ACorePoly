#pragma once

#include <vector>
#include <memory>

#include "Signal.h"

using namespace std;



namespace ACorePolyLib
{

	class Polynom
	{
	public:
		Polynom();
		Polynom(const Polynom & poly);
		Polynom(double c);
		Polynom(const vector<double> & coef);
		Polynom(int n, const double * coef);
		template<typename Iter>
		Polynom(int n, Iter it);
		~Polynom();

		void Reset();
		Polynom & operator=(const Polynom & poly);
		Polynom & operator+=(const Polynom & poly);
		Polynom & operator-=(const Polynom & poly);
		Polynom operator+(const Polynom & poly);
		Polynom operator-(const Polynom & poly);
		// ==,
		// +,+=,-,-=, *,*=,
		// + x, - x, ..
		// ...

		double operator[] (int ind) const;
		void SetCoef(int ind, double val);
		inline int Power() const;
		double operator() (double x) const;
		double MaxValueOnSegment(double start, double end, double dt) const; // analitic computing should be realized
		bool IsNull() const;

		Polynom Integrate() const; // defined integral with zero start value
		Polynom Differenciate() const;
		Polynom TimeShift(double shift) const;

		void FillSignal(Signal & signal) const;
		void AddToSignal(Signal & signal) const;
		// Transforms !!!
		// Serialise  !!!

	private:
		vector<double> m_coef;
		inline void DeleteTopZeros();
		inline void ExpandToIndex(int ind);
		
	};


	template<typename Iter>
	Polynom::Polynom(int n, Iter it)
	{
		for (int i = 0; i < n; ++i)
		{
			m_coef.push_back(*it);
			it++;
		}
	}
		

};
