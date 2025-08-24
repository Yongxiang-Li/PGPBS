#include "mex.h"
#include "matrix.h"

int findspan(int n, int p, double u, double *U) {

	int low, high, mid;
	// special case 
	if (u == U[n + 1]) return(n);
	// do binary search 
	low = p;
	high = n + 1;
	mid = (low + high) / 2;
	while (u < U[mid] || u >= U[mid + 1])  {
		if (u < U[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return(mid);
}

void basisfun(int i, double u, int p, double *U, double *N) {
	int j, r;
	double saved, temp;
	// work space 
	double *left = (double*)mxMalloc((p + 1)*sizeof(double));
	double *right = (double*)mxMalloc((p + 1)*sizeof(double));

	N[0] = 1.0;
	for (j = 1; j <= p; j++) {
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		saved = 0.0;

		for (r = 0; r < j; r++) {
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;   // mexPrintf("%f\n", N[j]);
	}
	mxFree(left);
	mxFree(right);
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int d = mxGetScalar(prhs[0]);
	int p = mxGetScalar(prhs[1]);
	double *knots = mxGetPr(prhs[2]);
	double *U = mxGetPr(prhs[3]);
	mwSize n = mxGetM(prhs[3]);
	mwSize nzmax = n*(d + 1);
	plhs[0] = mxCreateSparse(p, n, nzmax, mxREAL);
	mwIndex *ir = mxGetIr(plhs[0]);
	mwIndex *jc = mxGetJc(plhs[0]);
	double *pr = mxGetPr(plhs[0]);
	int currentIdx = 0;
	for (int c = 0; c<n; c++) {
		jc[c] = currentIdx;
		int s = findspan(p - 1, d, U[c], knots);
		basisfun(s, U[c], d, knots, pr + currentIdx);
		for (int r = 0; r<d + 1; r++) {
			ir[currentIdx] = s - d + r;
			currentIdx++;
		}
	}
	jc[n] = currentIdx;
}