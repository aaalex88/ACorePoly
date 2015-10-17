#include <stdafx.h>

#include <math.h>
#include <algorithm>


#include "Polynom.h"
#include "PolyLibHelpers.h"



namespace ACorePolyLib
{

	Polynom::Polynom()
	{
	}

	Polynom::Polynom(const Polynom & poly)
	{
		m_coef = poly.m_coef;
	}

	Polynom::Polynom(double c)
	{
		m_coef.push_back(c);
	}

	Polynom::Polynom(const vector<double> & coef)
	{
		m_coef = coef;
		DeleteTopZeros();
	}

	Polynom::Polynom(int n, const double * coef)
	{
		for (int i = 0; i < n; ++i)
			m_coef.push_back(coef[i]);
	}

	Polynom::~Polynom()
	{
	}


	void Polynom::Reset()
	{
		m_coef.clear();
	}

	Polynom & Polynom::operator=(const Polynom & poly)
	{
		m_coef = poly.m_coef;

		return *this;
	}

	Polynom & Polynom::operator+=(const Polynom & poly)
	{
		for (int i = 0; i <= Power(); ++i) {
			m_coef[i] += poly[i];
		}
		for (int i = Power() + 1; i <= poly.Power(); ++i) {
			m_coef.push_back(poly[i]);
		}
		return *this;
	}

	Polynom & Polynom::operator-=(const Polynom & poly)
	{
		for (int i = 0; i <= Power(); ++i) {
			m_coef[i] -= poly[i];
		}
		for (int i = Power() + 1; i <= poly.Power(); ++i) {
			m_coef.push_back(-poly[i]);
		}
		return *this;
	}

	Polynom Polynom::operator+(const Polynom & poly) const
	{
		return Polynom(*this) += poly;
	}

	Polynom Polynom::operator-(const Polynom & poly) const
	{
		return Polynom(*this) -= poly;
	}

	Polynom Polynom::operator*(const Polynom & poly) const
	{
		Polynom res;
		for (int i = 0; i <= Power() + poly.Power(); ++i) {
			for (int j = 0; j <= i; ++j) {
				res.SetCoef(i, res[i] + (*this)[j] * poly[i - j]);
			}
		}
		return res;
	}

	bool Polynom::operator==(const Polynom & poly) const
	{
		int minPow = (std::min)(Power(), poly.Power());
		int maxPow = (std::max)(Power(), poly.Power());
		const Polynom & maxPoly = ((*this).Power() == maxPow) ? *this : poly;

		for (int i = 0; i <= minPow; ++i) {
			if (!AlmostZero((*this)[i] - poly[i]))
				return false;
		}
		for (int i = minPow + 1; i <= maxPow; ++i) {
			if (!AlmostZero(maxPoly[i]))
				return false;
		}

		return true;
	}

	double Polynom::operator[] (int ind) const
	{
		if (ind < 0)
		{
			assert(!"Polynom::[] : negative index");
			return 0;
		}

		if ((unsigned int)ind >= m_coef.size())
			return 0;

		return m_coef[ind];
	}

	void Polynom::SetCoefNZ(int ind, double val) // TODO : delete function?
	{
		if (ind < 0)
		{
			assert(!"Polynom::SetCoef : negative index");
			return;
		}

		if (ind > Power() && AlmostZero(val))
			return;

		ExpandToIndex(ind);
		m_coef[ind] = val;
		DeleteTopZeros();
	}

	void Polynom::SetCoef(int ind, double val)
	{
		if (ind < 0) { // TODO : dublicated code
			assert(!"Polynom::SetCoef : negative index");
			return;
		}
		if (ind > Power() && val == 0.0)
			return;

		ExpandToIndex(ind);
		m_coef[ind] = val;
	}

	int Polynom::Power() const
	{
		int pow = m_coef.size() - 1;
		return (std::max)(pow, 0);
	}

	double Polynom::operator() (double x) const
	{
		double res = 0;
		for (size_t i = 0; i < m_coef.size(); ++i)
		{
			res += m_coef[i] * pow(x, i);
		}
		return res;
	}

	double Polynom::MaxValueOnSegment(double start, double end, double dt) const
	{
		double res = (*this)(start);
		for (double x = start; x <= end; x += dt)
		{
			res = fmax(res, (*this)(x));
		}
		return res;
	}

	bool Polynom::IsNull() const
	{
		return m_coef.size() == 0;
	}


	Polynom Polynom::Sqr() const
	{
		return *this * *this;
	}

	Polynom Polynom::Integrate() const
	{
		vector<double> v;
		v.push_back(0);
		for (size_t i = 0; i < m_coef.size(); ++i)
		{
			v.push_back(m_coef[i] / double(i+1));
		}
		return Polynom(v);
	}

	Polynom Polynom::Differenciate() const
	{
		vector<double> v;
		for (size_t i = 1; i < m_coef.size(); ++i)
		{
			v.push_back(m_coef[i] * i);
		}
		return Polynom(v);
	}

	double Polynom::L2Norm() const
	{
		double sq = fabs( (*this).Sqr().Integrate()(1.0) );
		double res = sqrt(sq);
		return res;
	}

	Polynom Polynom::TimeShift(double shift) const
	{
		vector<double> v;
		Polynom p(*this);
		for (int i = 0; i <= Power(); ++i)
		{
			v.push_back(p(shift) / double(Factorial(i)));
			p = p.Differenciate();
		}
		return Polynom(v);
	}


	void Polynom::FillSignal(Signal & signal) const
	{
		int ind = 0;
		int N = signal.GetDesc().N;
		double st = signal.GetDesc().startTime;
		double dt = signal.GetDesc().dt;
		signal.FillData(0, N, [&](){ 
			double res = (*this)( st + ind * dt );
			ind++;
			return res;
		}); // TODO : Test! are we sure it works correct?
	}

	void Polynom::AddToSignal(Signal & signal) const
	{
		const SignalDescription & d = signal.GetDesc();
		double * data = signal.GetData();
		for (int i = 0; i < d.N; ++i)
		{
			data[i] += (*this)(d.startTime + i * d.dt);
		}
	}


	inline void Polynom::DeleteTopZeros()
	{
		while (m_coef.size() > 0 && AlmostZero(m_coef[m_coef.size() - 1]) )
			m_coef.pop_back();
	}
	inline void Polynom::ExpandToIndex(int ind)
	{
		if (ind < 0) {
			assert(!"Incorrect input in Polynom::ExpandToIndex function: ind < 0");
			return;
		}

		while (m_coef.size() <= (unsigned int)ind)
			m_coef.push_back(0);
	}

};
