#include "./include/glbopts.h"
#include "./include/linalg.h"
#include "./include/amatrix.h"
#include "./include/abip.h"
#include "./include/util.h"

typedef struct RData
{
    ABIPData *d;
    ABIPSolution *sol;
} RData;

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

void copy_sparse_A(abip_float *dist, const abip_float *sigma, abip_int *rInd, abip_int *cInd, abip_int dim)
{
    abip_int k = 0;
    abip_int len_Sig = dim * dim;
    abip_int len_dist = 2 * dim * dim + 3 * dim;

    for (abip_int i = 0; i < dim; i++)
    {
        for (abip_int j = 0; j < dim; j++)
        {
            dist[k] = sigma[j * dim + i];
            dist[k + len_Sig] = -dist[k];
            cInd[k] = j;
            rInd[k] = i;
            cInd[k + len_Sig] = j + dim;
            rInd[k + len_Sig] = i;
            k = k + 1;
        }
    }
    k = 2 * k;
    for (abip_int i = (2 * len_Sig); i < len_dist; i++)
    {
        dist[i] = 1;
    }
    for (abip_int i = 0; i < dim; i++)
    {
        rInd[k] = i;
        rInd[k + dim] = dim + i;
        rInd[k + 2 * dim] = dim + i;
        cInd[k] = 2 * dim + i;
        cInd[k + dim] = 2 * dim + i;
        cInd[k + 2 * dim] = 3 * dim + i;
        k = k + 1;
    }
}

RData *trans_R_data(abip_float *Sigma, abip_float *b, abip_float *c, abip_float *x, abip_float *y, abip_float *s, abip_int *p)
{
    RData *data;
    data = malloc(sizeof(RData));
    data->d = malloc(sizeof(ABIPData));
    data->sol = malloc(sizeof(ABIPSolution));

    abip_int m = 2 * (*p);
    abip_int n = 4 * (*p);

    abip_int len_Sig = (*p) * (*p);
    abip_int len_b = m;
    abip_int len_c = n;
    abip_int len_x = n;
    abip_int len_s = len_x;
    abip_int len_y = m;

    abip_int len_A = 2 * len_Sig + 3 * (*p);

    // copy sparse matrix A to data
    abip_float *A = (abip_float *)malloc(len_A * sizeof(abip_float));
    abip_int *rInd = (abip_int *)malloc(len_A * sizeof(abip_int));
    abip_int *cInd = (abip_int *)malloc(len_A * sizeof(abip_int));

    copy_sparse_A(data->d->A->x, Sigma, rInd, cInd, (*p));
    // data->d->A->x = A;
    data->d->A->m = m;
    data->d->A->n = n;
    data->d->A->i = rInd;
    data->d->A->p = cInd; //! the format of sparse matrix in ABIP is different

    //copy b and c to data
    data->d->b = (abip_float *)malloc(len_b * sizeof(abip_float));
    data->d->c = (abip_float *)malloc(len_c * sizeof(abip_float));
    memcpy(data->d->b, b, len_b * sizeof(abip_float));
    memcpy(data->d->c, c, len_c * sizeof(abip_float));

    data->sol->s = (abip_float *)s;
    data->sol->x = (abip_float *)x;
    data->sol->x = (abip_float *)y;

    return data;
}

void abip_R(RData *data)
{
    abip_int i;
    abip_int status;

    ABIPSolution sol = {0}; //TODO: free the memory of x, y and s
    ABIPInfo info;

    data->d->stgs->warm_start = parse_warm_start(data->sol->x, &(sol.x), data->d->n);
    data->d->stgs->warm_start |= parse_warm_start(data->sol->y, &(sol.y), data->d->m);
    data->d->stgs->warm_start |= parse_warm_start(data->sol->s, &(sol.s), data->d->n);

    status = ABIP(main)(data->d, &sol, &info);
}
