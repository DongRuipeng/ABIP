#include <stdio.h>
#include <stdlib.h>

void print_mat(double *Sigma, unsigned *n)
{
    const unsigned len = *n;
    for (unsigned i = 0; i < len; i++)
    {
        printf("%f\n", Sigma[i]);
    }
    
}