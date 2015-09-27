#include <stdafx.h>

#include "Optimisation.h"
#include "alglib/solvers.h"


using namespace alglib;

namespace ACorePolyLib
{

	double * optimise(int N, int dim, const double * basis, const double * signal)
	{
		return nullptr;
	}

	bool decompose(int N, int dim, const double * basis, const double * signal, double * res)
	{
		real_2d_array A;
		A.setcontent(dim, N, basis);

		real_1d_array B;
		B.setcontent(N, signal);

		real_1d_array X;
		X.setlength(dim);

		ae_int_t info;
		densesolverlsreport rep;

		rmatrixsolvels(A, dim, N, B, 0.0, info, rep, X);

		for (int i = 0; i < dim; ++i) {
			res[i] = X[i];
		}
		
		return (info == 1); // whether task is solved
	}
	
}