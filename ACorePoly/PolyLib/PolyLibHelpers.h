#pragma once


namespace ACorePolyLib
{

	class Interface
	{
	public:
		virtual ~Interface() {}
	};


	template<typename RetType>
	class ISimpleFunctor : public Interface
	{
	public:
		RetType operator() () = 0;
	};


	template<typename RetType, typename Callable>
	class SimpleFunctor : public ISimpleFunctor<RetType>
	{
	private:
		Callable m_fun;

	public:
		SimpleFunctor(Callable fun) { m_fun = fun; }
		RetType operator() () { return m_fun(); }
	};


	const int n = 7; // testing no extern vars
	const double eps = 10e-15;
	const double pi = 3.1415926535897932384626433832795f;
	const double pi_2 = 6.283185307179596476925286766559f;

	bool AlmostZero(double d);
	double ComputePhase(double cosAmp, double sinAmp);

};
