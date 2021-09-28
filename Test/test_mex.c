#include "matrix.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *A_mex = (mxArray *)prhs[0];
    // const size_t *nir = (size_t *) prhs[1];
    // const size_t *njc = (size_t *) prhs[2];
    // printf("*nir = %i,\t *njc = %i\n", *nir, *njc);
    
    mwIndex *jc = mxGetJc(A_mex);
    mwIndex *ir = mxGetIr(A_mex);
    mxDouble *x = mxGetPr(A_mex);
    printf("the list of jc: \n");
    for (size_t i = 0; i < 3 + 1; i++)
    {
        printf("jc[%i] = %i\n", i, jc[i]);
    }
    printf("the list of ir: \n");
    for (size_t i = 0; i < 6; i++)
    {
        printf("ir[%i] = %i\n", i, ir[i]);
    }
    printf("the list of x: \n");
    for (size_t i = 0; i < 10; i++)
    {
        printf("x[%i] = %f\n", i, x[i]);
    }
    
}