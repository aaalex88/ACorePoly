#pragma once

#include "PolyLibHelpers.h"



namespace ACorePolyLib
{

	struct SignalDescription
	{
		SignalDescription() : N(0), dt(0), startTime(0) {}
		SignalDescription(int _N, double _dt, double _st);

		int		N;
		double	dt;
		double	startTime;
	};

	class ISignal : public Interface
	{
	public:
		virtual double operator[] (int ind) const;
		virtual const SignalDescription & GetDescr() const;
	};

	class Signal
	{
	private:
		SignalDescription m_desc;
		double * m_signal;
		bool m_autoDelete; // delete m_signal in destructor;

	public:
		Signal(const SignalDescription & desc);						// create empty array
		Signal(const SignalDescription & desc, double * signal);	// use existing array as data
		~Signal();

		void FillData(int startIndex, int endIndex, const double * data); // copy data to signal: m_signal[st..en-1] = data[0..en-st-1]
		template<typename DataStream>
		void FillData(int startIndex, int endIndex, DataStream data); // m_signal[i] = data();
		void Reset();

		void Add(const double * data);
		void Substract(const double * data);
		void Invert();

		inline double operator[] (int ind) const;
		inline void Set(int ind, double val);
		const SignalDescription & GetDesc() const;
		double * GetData();
		const double * GetData() const;
		double GetNorm() const;
	};


	template<typename DataStream>
	void Signal::FillData(int startIndex, int endIndex, DataStream data)
	{
		if(startIndex < 0 || endIndex >= m_desc.N) {
			// assert!!!
			return;
		}

		for (int i = 0; i < endIndex - startIndex; ++i)
			m_signal[startIndex + i] = data();
	}

};
