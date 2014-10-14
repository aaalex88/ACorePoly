#pragma once


namespace ACorePolyLib
{

	class Interface
	{
	public:
		virtual ~Interface();
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


	static const int n = 7; // testing no extern vars

	bool AlmostZero(double d);

};
