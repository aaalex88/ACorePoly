#pragma once


namespace ACorePolyLib
{
	double * optimise(int N, int dim, const double * basis, const double * signal);

	//template<typename Fun>
	//vector<double> minimise(vector<double> & st, Fun f); // double optimize(double *)
}