#pragma once


namespace ACorePolyLib
{
	class Polynom;
	struct ACore;

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
	double Factorial(int n);

	inline double fdiv(int x, int y) { return double(x) / double(y); }
//	inline int max(int x, int y) { return x > y ? x : y; }
//	inline int min(int x, int y) { return x < y ? x : y; }


	double random(double max);
	double randomRange(double start, double end);
	
	void RandomPolynom(Polynom & res, int pow, double startVal, double delta);
	void RandomCore(ACore & res, double fr, int numAmpl, int ampPower);

};
