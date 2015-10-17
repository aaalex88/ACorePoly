#include <stdafx.h>

#include "Optimisation.h"
#include "alglib/solvers.h"


#include <gsl\gsl_linalg.h>

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

	bool decomposeAlt(int N, int dim, const double * basis, const double * signal, double * res)
	{
		real_2d_array AT;
		AT.setcontent(dim, N, basis);

		real_1d_array B;
		B.setcontent(N, signal);

		// X = (AT*A)^(-1)*AT*B

		real_1d_array X;
		X.setlength(dim);

		real_2d_array A;
		A.setlength(N, dim);
		rmatrixtranspose(dim, N, AT, 0, 0, A, 0, 0);

		real_2d_array ATA1;
		ATA1.setlength(dim, dim);
		rmatrixgemm(dim, dim, N, 1.0, AT, 0, 0, 0, A, 0, 0, 0, 0.0, ATA1, 0, 0);
		ae_int_t info1;
		matinvreport rep1;
		rmatrixinverse(ATA1, info1, rep1);
		if (info1 != 1) {
			throw;
			return false;
		}

		real_2d_array ATAAT;
		ATAAT.setlength(dim, N);
		rmatrixgemm(dim, N, dim, 1.0, ATA1, 0, 0, 0, AT, 0, 0, 0, 0.0, ATAAT, 0, 0);

		rmatrixmv(dim, N, ATAAT, 0, 0, 0, B, 0, X, 0);

		
		for (int i = 0; i < dim; ++i) {
			res[i] = X[i];
		}

		return true; // whether task is solved
	}
	
	bool decomposeGSL(int N, int dim, const double * basis, const double * signal, double * res)
	{
		assert(dim <= N);

		int i, j;
		gsl_matrix * A = gsl_matrix_alloc(N, dim);
		gsl_vector * B = gsl_vector_alloc(N);
		gsl_vector * X = gsl_vector_alloc(dim);
		gsl_vector * Residual = gsl_vector_alloc(N);
		gsl_vector * tau = gsl_vector_alloc(dim);

		for (i = 0; i < dim; ++i) {
			for (j = 0; j < N; ++j) {
				gsl_matrix_set(A, j, i, basis[i*N + j]);
			}
		}

		for (j = 0; j < N; ++j) {
			gsl_vector_set(B, j, signal[j]);
		}

		gsl_linalg_QR_decomp(A, tau); // TODO : error codes handling
		gsl_linalg_QR_lssolve(A, tau, B, X, Residual);


		for (i = 0; i < dim; ++i) {
			res[i] = gsl_vector_get(X, i);
		}
		gsl_matrix_free(A);
		gsl_vector_free(B);
		gsl_vector_free(X);
		gsl_vector_free(Residual);
		gsl_vector_free(tau);
		return true;
	}

}