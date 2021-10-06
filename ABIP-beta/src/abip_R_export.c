#include "glbopts.h"
#include "linalg.h"
#include "amatrix.h"
#include "abip.h"
#include "util.h"
#include <time.h>
// #include <R.h>
// #include <R_ext/Lapack.h>
// #include <R_ext/BLAS.h>

ABIPMatrix *trans_to_sparse_mat(abip_float *Sigma, abip_int *p)
{
    abip_int len_sig = (*p) * (*p);
    abip_int len_x = 2 * len_sig + 3 * (*p);

    // declare the memory
    abip_float *x = malloc(len_x * sizeof(abip_float));
    abip_int *jc = malloc(4 * (*p) * sizeof(abip_int) + 1);
    abip_int *ir = malloc(len_x * sizeof(abip_int));

    // copy the value
    memcpy(x, Sigma, len_sig * sizeof(abip_float));
    for (abip_int i = len_sig; i < 2 * len_sig; i++)
    {
        x[i] = -Sigma[i - len_sig];
    }
    for (abip_int i = 2 * len_sig; i < len_x; i++)
    {
        x[i] = 1;
    }

    // copy the jc
    jc[0] = 0;
    for (abip_int i = 1; i < 2 * (*p) + 1; i++)
    {
        jc[i] = jc[i - 1] + (*p);
    }
    for (abip_int i = 2 * (*p) + 1; i < 3 * (*p) + 1; i++)
    {
        jc[i] = jc[i - 1] + 2;
    }
    for (abip_int i = 3 * (*p) + 1; i < 4 * (*p) + 1; i++)
    {
        jc[i] = jc[i - 1] + 1;
    }

    // copy the ir
    abip_int *tmp = malloc((*p) * sizeof(abip_int));
    // init temp vector
    for (abip_int i = 0; i < (*p); i++)
    {
        tmp[i] = i;
    }

    for (abip_int i = 0; i < 2 * (*p); i++)
    {
        memcpy(&ir[i * (*p)], tmp, (*p) * sizeof(abip_int));
    }
    for (abip_int i = 2 * len_sig; i < 2 * len_sig + 2 * (*p); i++)
    {
        ir[i] = i / 2 % (*p) + (*p) * (i % 2);
    }
    for (abip_int i = 2 * len_sig + 2 * (*p); i < len_x; i++)
    {
        ir[i] = i % (2 * len_sig + 2 * (*p)) + *p;
    }

    ABIPMatrix *A = malloc(sizeof(ABIPMatrix));
    A->x = x;
    A->p = jc;
    A->i = ir;
    A->m = 2 * (*p);
    A->n = 4 * (*p);

    free(tmp);

    return A;
}

void free_abip_data(ABIPData *d)
{
    if (d)
    {
        if (d->A)
        {
            if (d->A->i)
            {
                free(d->A->i);
            }
            if (d->A->p)
            {
                free(d->A->p);
            }
            if (d->A->x)
            {
                free(d->A->x);
            }

            free(d->A);
        }
        if (d->b)
        {
            free(d->b);
        }
        if (d->c)
        {
            free(d->c);
        }
        if (d->stgs)
        {
            free(d->stgs);
        }

        free(d);
    }
}

void free_sol(ABIPSolution sol)
{
    if (sol.x)
    {
        free(sol.x);
    }
    if (sol.y)
    {
        free(sol.y);
    }
    if (sol.s)
    {
        free(sol.s);
    }
}

void abip_R(abip_float *Sigma, abip_int *p, abip_float *lambda, abip_int *nlambda, abip_float *Omega, abip_float *time)
{
    abip_int status;
    ABIPInfo info;

    // construct A
    ABIPMatrix *A = trans_to_sparse_mat(Sigma, p);

    // construct c
    abip_float *c = malloc(A->n * sizeof(abip_float));
    for (abip_int i = 0; i < 2 * (*p); i++)
    {
        c[i] = 1;
    }
    for (abip_int i = 2 * (*p); i < A->n; i++)
    {
        c[i] = 0;
    }

    // construct b_base
    abip_float b_base[A->m]; // = malloc(A->m * sizeof(abip_float));
    for (abip_int i = 0; i < (*p); i++)
    {
        b_base[i] = 1;
    }
    for (abip_int i = (*p); i < A->m; i++)
    {
        b_base[i] = 2;
    }

    ABIPSolution sol = {0};

    // declare the memory of b
    abip_float *b = malloc(A->m * sizeof(abip_float));

    // declare the memory of d and init it
    ABIPData *d = malloc(sizeof(ABIPData));
    d->A = A;
    d->b = b;
    d->c = c;
    d->m = A->m;
    d->n = A->n;
    d->sp = (abip_float)A->p[A->n] / (A->m * A->n);
    d->stgs = (ABIPSettings *)malloc(sizeof(ABIPSettings));
    ABIP(set_default_settings)
    (d);
    d->stgs->verbose = 0;
    d->stgs->warm_start = 0;

    clock_t start, end;
    // double time;
    start = clock();

    for (abip_int j = 0; j < (*p); j++)
    {
        for (abip_int l = 1; l < (*nlambda + 1); l++)
        {
            // update b
            memcpy(b, b_base, A->m * sizeof(abip_float));
            ABIP(scale_array)
            (b, lambda[l - 1], A->m);
            b[j] += 1;

            sol.x = malloc(d->n * sizeof(abip_float));
            sol.y = malloc(d->m * sizeof(abip_float));
            sol.s = malloc(d->n * sizeof(abip_float));

            status = ABIP(main)(d, &sol, &info);

            for (abip_int k = 0; k < (*p); k++)
            {
                Omega[k + j * (*p) + (l - 1) * (*p) * (*p)] = sol.x[k] - sol.x[k + (*p)];
            }
            free_sol(sol);
        }
    }
    end = clock();
    *time = ((double)(end - start)) / CLOCKS_PER_SEC;

    free_abip_data(d);
}
