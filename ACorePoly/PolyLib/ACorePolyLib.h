#pragma once

#include <vector>
#include <memory>


#include "PolyLibHelpers.h"
#include "Polynom.h"
#include "Signal.h"
#include "ACore.h"
#include "Descriptions.h"
#include "SignalSolver.h"
#include "ILinearSolver.h"
#include "AmplSolver.h"
#include "ACoreSolver.h"
#include "SegmentSolver.h"
#include "SegmentDescription.h"
#include "DescriptionBuilder.h"
#include "Optimisation.h"

using namespace std;
using namespace std::tr1;




namespace ACorePolyLib
{

	

};










/*
 *
 *
 *
 * TEMPLATE SECTION
 *
 *
 *
 */

/*

template <typename Segment, typename SegDescr>
void AnalyseSignal(SegDescr (*solve)(Segment), Segment (*generateSegment)(void), void (handleResult)(SegDescr))
{
	while ( (auto seg = generateSegment()) != nullptr) // ?? shared_ptr<Segment> ??
	{
		auto res = solve(seg);
		handleResult(res);
	}
}
//*/


/*
template <typename Solver>
class Solution
{
public:
	template<typename Segment, typename SegDescr>
	static void AnalyseSignal(Segment (*generateSegment)(void), void (handleResult)(SegDescr))
	{
		Solver solver;
		::AnalyseSignal<Segment, SegDescr>([&](Segment seg) { return solver.Solve(seg);}, generateSegment, handleResult);
	}
};
//*/

/*
template<typename Segment, typename SegDescr, typename SolveFactory>
void AnalyseSignal(SolveFactory solveFactory, Segment (*generateSegment)(void), void (handleResult)(SegDescr))
{
	SegDescr (*solve)(Segment) = solveFactory();
	AnalyseSignal(solve, generateSegment, handleResult);
}
//*/


template <typename Segment, typename SegDescr>
struct SolveFactoryRetType 
{
	typedef SegDescr type (Segment);
};


template <typename Segment, typename SegDescr, typename Computations, typename Minimizator,typename PostProcesssor> // typename OnPointProcessor
// vector<double> -> Point
// BuildSignal, seg.params, Norm
typename SolveFactoryRetType<Segment,SegDescr>::type 
SolveFactory(Computations InitComputations, Minimizator minimize, PostProcesssor PostProcess) // OnPointProcessor OnPointSelected
{
	SegDescr oldDesc;
	return [=] (Segment seg) -> SegDescr {
		vector<double> st_x;
		LinearSolver solve;
		InitComputations(seg, oldDesc, st_x, solve);

		vector<double> res = minimize(st_x, [&](vector<double> x){
			SegDescr desc = solve(x);
	//		Segment resSeg = BuildSignal(desc/*, seg.Params*/) - seg;
	//		return Norm(resSeg);
		}); // ,OnPointSelected
		SegDescr desc = solve(res);
		PostProcess(desc);
		oldDesc = desc;
		return desc;
	};
}









