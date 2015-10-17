#include <stdafx.h>

#include <math.h>
#include "Signal.h"


namespace ACorePolyLib
{

	Signal::Signal(const SignalDescription & desc) : m_desc(desc), m_autoDelete(false)
	{
		m_signal = new double[m_desc.N];

		if (m_signal != nullptr) {
			m_autoDelete = true;
			Reset();
		}
	}

	Signal::Signal(const SignalDescription & desc, double * signal) : m_desc(desc), m_signal(signal), m_autoDelete(false)
	{
	}

	Signal::~Signal()
	{
		if (m_autoDelete && m_signal)
			delete[] m_signal;
	}

	void Signal::FillData(int startIndex, int endIndex, const double * data)
	{
		if(startIndex < 0 || endIndex >= m_desc.N) {
			assert(!"Signal::FillData : Index Out Of Range, Operation Cancelled");
			return;
		}

		for (int i = 0; i < endIndex - startIndex; ++i)
			m_signal[startIndex + i] = data[i];
	}

	void Signal::Reset()
	{
		for (int i = 0; i < m_desc.N; ++i)
			m_signal[i] = 0;
	}

	void Signal::Add(const double * data)
	{
		for (int i = 0; i < m_desc.N; ++i)
			m_signal[i] += data[i];
	}

	void Signal::Substract(const double * data)
	{
		for (int i = 0; i < m_desc.N; ++i)
			m_signal[i] -= data[i];
	}

	void Signal::Invert()
	{
		for (int i = 0; i < m_desc.N; ++i)
			m_signal[i] = -m_signal[i];
	}

	inline double Signal::operator[] (int ind) const
	{
		if (ind < 0 || ind > m_desc.N)
		{
			assert(!"Signal::[] : Index Out Of Range, Operation Cancelled");
			return 0;
		}

		return m_signal[ind];
	}

	inline void Signal::Set(int ind, double val)
	{
		if (ind < 0 || ind >= m_desc.N)
		{
			assert(!"Signal::Set : Index Out Of Range, Operation Cancelled");
			return;
		}

		m_signal[ind] = val;
	}

	const SignalDescription & Signal::GetDesc() const
	{
		return m_desc;
	}

	double * Signal::GetData()
	{
		return m_signal;
	}

	const double * Signal::GetData() const
	{
		return m_signal;
	}

	double Signal::GetNorm() const
	{
		double res = 0;

		for (int i = 0; i < m_desc.N; ++i)
		{
			res += m_signal[i] * m_signal[i];
		}

		return sqrt( res * m_desc.dt );
	}


};