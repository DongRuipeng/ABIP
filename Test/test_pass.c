#include <stdio.h>
void add_mat(double *A, double *B, double *C, unsigned *n)
{
    for (unsigned i = 0; i < (*n) * (*n); i++)
    {
        C[i] = A[i] + B[i];
    }
}