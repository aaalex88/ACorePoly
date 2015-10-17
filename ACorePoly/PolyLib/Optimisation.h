#pragma once


namespace ACorePolyLib
{
	double * optimise(int N, int dim, const double * basis, const double * signal);

	bool decompose(
		int N,					// length of signal
		int dim,				// number of basis vectors
		const double * basis,	// array [dim][N] of basis vectors (N*d + n)
		const double * signal,	// decomposing signal [N]
		double * res			// result vector [dim] of coeficients of basis vector that gives best approximation of signal
		);

	bool decomposeAlt(int N, int dim, const double * basis, const double * signal, double * res);

	bool decomposeGSL(int N, int dim, const double * basis, const double * signal, double * res);


	//template<typename Fun>
	//vector<double> minimise(vector<double> & st, Fun f); // double optimize(double *)
}