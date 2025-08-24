#include "mex.h"
#include "matrix.h"
#include <math.h>

int binarySearch(double *K, int m, double u) {
	int a = 0, b = m;
	while (b > a) {
		int s = (a + b) / 2;
		if (u <= K[s]) {
			b = s;
		}
		else {
			a = s + 1;
		}
	}
	return (a + b) / 2;
}

// Helper function to calculate the periodic B-spline basis function
void basisfun(int s, double u, int d, double *K, double *N, double *N0) {
	int i, j;
	double left, right;

	N0[s - 1] = 0.0;
	N0[s] = 1.0;

	for (j = 1; j <= d; j++) {
		for (i = s; i <= s + j; i++) {
			left = (1.0 + u - K[i - j - 1]) / (K[i - 1] - K[i - j - 1]);
			right = (K[i] - u - 1.0) / (K[i] - K[i - j]);
			N[i] = left * N0[i - 1] + right * N0[i];
		}
		for (i = s; i <= s + j; i++) {
			N0[i] = N[i];
		}
	}
}

// Main MEX function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Input parameters
	int d = mxGetScalar(prhs[0]);
	int p = mxGetScalar(prhs[1]);
	double *K = mxGetPr(prhs[2]);
	double *U = mxGetPr(prhs[3]);

	// Dimensions and output
	mwSize nu = mxGetNumberOfElements(prhs[3]);
	mwSize m = mxGetNumberOfElements(prhs[2]);

	plhs[0] = mxCreateSparse(p, nu, p * nu, mxREAL);
	double *pr = mxGetPr(plhs[0]);
	mwIndex *ir = mxGetIr(plhs[0]);
	mwIndex *jc = mxGetJc(plhs[0]);

	// Allocate temporary basis function array
	double *N = (double *)mxCalloc(3*m, sizeof(double));
	double *N0 = (double *)mxCalloc(3*m, sizeof(double));


	// Expand K for periodic indexing
	double *K3 = (double *)mxCalloc(3 * m, sizeof(double));
	for (int i = 0; i < m; i++) {
		K3[i] = K[i];
		K3[i + m] = K[i] + 1.0;
		K3[i + 2 * m] = K[i] + 2.0;
	}

	// Fill sparse matrix
	int currentIdx = 0;
	for (mwSize col = 0; col < nu; col++) {
		int s = binarySearch(K3, m, U[col]);
		s += m;

		basisfun(s, U[col], d, K3, N, N0);
		jc[col] = currentIdx;
		for (int row = 2*m; row <= s + d; row++) {
			ir[currentIdx] = row - 2 * m;
			pr[currentIdx] = N[row];
			currentIdx++;
            N0[row] = 0.0;
		}
		for (int row = s; row < 2*m && row <= s + d; row++) {
			ir[currentIdx] = row - m;
			pr[currentIdx] = N[row];
			currentIdx++;
            N0[row] = 0.0;
		}
	}
	jc[nu] = currentIdx;

	// Free allocated memory
	mxFree(K3);
	mxFree(N);
	mxFree(N0);
}
