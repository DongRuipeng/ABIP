#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_mat(double *Sigma, unsigned *p)
{
    unsigned len_sig = (*p) * (*p);
    unsigned len_x = 2 * len_sig + 3 * (*p);

    // declare the memory
    double *x = malloc(len_x * sizeof(double));
    unsigned *jc = malloc(4 * (*p) * sizeof(unsigned) + 1);
    unsigned *ir = malloc(len_x * sizeof(unsigned));
    unsigned *tmp = malloc((*p) * sizeof(unsigned));

    // init temp vector
    for (unsigned i = 0; i < (*p); i++)
    {
        tmp[i] = i;
    }

    // copy the value
    memcpy(x, Sigma, len_sig * sizeof(double));
    for (unsigned i = len_sig; i < 2 * len_sig; i++)
    {
        x[i] = -Sigma[i - len_sig];
    }
    for (unsigned i = 2 * len_sig; i < len_x; i++)
    {
        x[i] = 1;
    }

    // copy the jc
    jc[0] = 0;
    for (unsigned i = 1; i < 2 * (*p) + 1; i++)
    {
        jc[i] = jc[i - 1] + (*p);
    }
    for (unsigned i = 2 * (*p) + 1; i < 3 * (*p) + 1; i++)
    {
        jc[i] = jc[i - 1] + 2;
    }
    for (unsigned i = 3 * (*p) + 1; i < 4 * (*p) + 1; i++)
    {
        jc[i] = jc[i - 1] + 1;
    }

    // copy the ir
    for (unsigned i = 0; i < 2 * (*p); i++)
    {
        memcpy(&ir[i * (*p)], tmp, (*p) * sizeof(unsigned));
    }
    for (unsigned i = 2 * len_sig; i < 2 * len_sig + 2 * (*p); i++)
    {
        ir[i] = i / 2 % (*p) + (*p) * (i % 2);
    }
    for (unsigned i = 2 * len_sig + 2 * (*p); i < len_x; i++)
    {
        ir[i] = i % (2 * len_sig + 2 * (*p)) + *p;
    }

    // print the result for test
    printf("the list of x: \n");
    for (unsigned i = 0; i < len_x; i++)
    {
        printf("%f\n", x[i]);
    }
    printf("the list of ir: \n");
    for (unsigned i = 0; i < len_x; i++)
    {
        printf("%i\n", ir[i]);
    }
    printf("the list of jc: \n");
    for (unsigned i = 0; i < 4 * (*p) + 1; i++)
    {
        printf("%i\n", jc[i]);
    }

    free(x);
    free(jc);
    free(ir);
    free(tmp);
}