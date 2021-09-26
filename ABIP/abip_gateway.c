#include <stdio.h>
#include <stdlib.h>

void sp_mat(double *Sigma, unsigned *p)
{
    double *A = (double *)malloc((2 * (*p) * (*p) + 3 * (*p)) * sizeof(double));
    unsigned *cInd = (unsigned *)malloc((2 * (*p) * (*p) + 3 * (*p)) * sizeof(unsigned));
    unsigned *rInd = (unsigned *)malloc((2 * (*p) * (*p) + 3 * (*p)) * sizeof(unsigned));
    unsigned k = 0;
    for (unsigned i = 0; i < (*p); i++)
    {
        for (unsigned j = 0; j < (*p); j++)
        {
            A[k] = Sigma[j * (*p) + i];
            A[k + (*p) * (*p)] = -A[k];
            cInd[k] = j;
            rInd[k] = i;
            cInd[k + (*p) * (*p)] = j + (*p);
            rInd[k + (*p) * (*p)] = i;
            k = k + 1;
        }
    }
    k = 2 * k;
    for (unsigned i = (2 * (*p) * (*p)); i < (2 * (*p) * (*p) + 3 * (*p)); i++)
    {
        A[i] = 1;
    }
    
    for (unsigned i = 0; i < (*p); i++)
    {
        rInd[k] = i;
        rInd[k + (*p)] = (*p) + i;
        rInd[k + 2 * (*p)] = (*p) + i;
        cInd[k] = 2 * (*p) + i;
        cInd[k + (*p)] = 2 * (*p) + i;
        cInd[k + 2 * (*p)] = 3 * (*p) + i;
        k = k + 1;
    }

    for (unsigned i = 0; i < (2 * (*p) * (*p) + 3 * (*p)); i++)
    {
        printf("location: (%i, %i);\t value: %f \n", rInd[i], cInd[i], A[i]);
    }

    free(A);
    free(cInd);
    free(rInd);
}