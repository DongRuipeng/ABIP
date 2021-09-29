#include "glbopts.h"
#include "linalg.h"
#include "amatrix.h"
#include "abip.h"
#include "util.h"

abip_int parse_warm_start(const abip_float *p_mex, abip_float **p, abip_int len)
{
    *p = (abip_float *)abip_calloc(len, sizeof(abip_float));
    if (p_mex == ABIP_NULL)
    {
        return 0;
    }
    else
    {
        memcpy(*p, p_mex, len * sizeof(abip_float));
        return 1;
    }
}

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

ABIPData *construt_abip_data(ABIPMatrix *A, abip_float *b, abip_float *c)
{
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

    return d;
}

void free_r(ABIPData *d)
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

void abip_R(abip_float *Sigma, abip_int *p, abip_float *b, abip_float *c, ABIPSolution *init_sol)
{
    // abip_int i;
    abip_int status;

    ABIPMatrix *A = trans_to_sparse_mat(Sigma, p);
    ABIPData *d = construt_abip_data(A, b, c); //! A and c will not changed with different lambda, but the b does !

    ABIPSolution sol = {0}; //TODO: free the memory of x, y and s
    ABIPInfo info;

    d->stgs->warm_start = parse_warm_start(init_sol->x, &(sol.x), d->n); 
    //TODO: modify the function (a vector including pointer to sol.x,y,s)
    d->stgs->warm_start |= parse_warm_start(init_sol->y, &(sol.y), d->m);
    d->stgs->warm_start |= parse_warm_start(init_sol->s, &(sol.s), d->n);

    status = ABIP(main)(d, &sol, &info);

    free_r(d); //! Once free the memmory, all data will be destroyed in the next computation !
}
