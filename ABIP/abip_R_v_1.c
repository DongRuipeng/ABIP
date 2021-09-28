#include "./include/glbopts.h"
#include "./include/linalg.h"
#include "./include/amatrix.h"
#include "./include/abip.h"
#include "./include/util.h"

void free_mex(ABIPData *d);

typedef struct rArray
{
    ABIPMatrix *A;
    double *b;
    double *c;
    const size_t length_b;
    const size_t length_c;
} rArray;

typedef struct rSetting
{
    unsigned normalize;
    double scale;
    double rho_y;
    double sparsity_ratio;

    unsigned max_ipm_iters;
    unsigned max_admm_iters;
    double eps;
    double alpha;
    double cg_rate; /* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */

    unsigned adaptive;
    double eps_cor;
    double eps_pen;

    unsigned verbose;    /* boolean, write out progress: 1 */
    unsigned warm_start; /* boolean, warm start (put initial guess in ABIPSolution struct): 0 */

    unsigned adaptive_lookback;
} rSetting;

abip_int parse_warm_start(const abip_float *p_mex, abip_float **p, abip_int len)
{
    *p = (abip_float *)abip_calloc(len, sizeof(abip_float));
    if (p_mex == ABIP_NULL)
    {
        return 0;
    }
    else if (mxIsSparse(p_mex) || (abip_int)*mxGetDimensions(p_mex) != len)
    {
        abip_printf("Error in warm-start (the input vectors should be dense and of correct size), running without full warm-start");
        return 0;
    }
    else
    {
        memcpy(*p, mxGetPr(p_mex), len * sizeof(abip_float));
        return 1;
    }
}

#if !(DLONG > 0)
abip_int *cast_to_abip_int_arr(mwIndex *arr, abip_int len)
{
    abip_int i;
    abip_int *arr_out = (abip_int *)abip_malloc(sizeof(abip_int) * len);
    for (i = 0; i < len; i++)
    {
        arr_out[i] = (abip_int)arr[i];
    }
    return arr_out;
}
#endif

#if SFLOAT > 0
abip_float *cast_to_abip_float_arr(double *arr, abip_int len)
{
    abip_int i;
    abip_float *arr_out = (abip_float *)abip_malloc(sizeof(abip_float) * len);
    for (i = 0; i < len; i++)
    {
        arr_out[i] = (abip_float)arr[i];
    }
    return arr_out;
}

double *cast_to_double_arr(abip_float *arr, abip_int len)
{
    abip_int i;
    double *arr_out = (double *)abip_malloc(sizeof(double) * len);
    for (i = 0; i < len; i++)
    {
        arr_out[i] = (double)arr[i];
    }
    return arr_out;
}
#endif

void set_output_field(mxArray **pout, abip_float *out, abip_int len)
{
    *pout = mxCreateDoubleMatrix(0, 0, mxREAL);
#if SFLOAT > 0
    mxSetPr(*pout, cast_to_double_arr(out, len));
    abip_free(out);
#else
    mxSetPr(*pout, out);
#endif
    mxSetM(*pout, len);
    mxSetN(*pout, 1);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    abip_int i;
    abip_int status;

    ABIPData *d;
    ABIPSolution sol = {0};
    ABIPInfo info;
    ABIPMatrix *A;

    const rArray *data;
    const ABIPMatrix *A_mex;
    const double *b_mex;
    const double *c_mex;

    const rSetting *settings;

    const mwSize one[1] = {1};
    const int num_info_fields = 13;
    const char *info_fields[] = {"status", "ipm_iter", "admm_iter", "mu", "pobj", "dobj", "resPri", "resDual", "relGap", "resInfeas", "resUnbdd", "setupTime", "solveTime"};

#if EXTRA_VERBOSE > 0
    abip_printf("size of mwSize = %i\n", (int)sizeof(mwSize));
    abip_printf("size of mwIndex = %i\n", (int)sizeof(mwIndex));
#endif

    d = (ABIPData *)mxMalloc(sizeof(ABIPData));
    d->stgs = (ABIPSettings *)mxMalloc(sizeof(ABIPSettings));
    data = prhs[0];

    A_mex = data->A;

    if (A_mex == ABIP_NULL)
    {
        abip_free(d);
        printf("ABIPData struct must contain a matrix 'A'. \n");
    }

    b_mex = data->b;

    if (b_mex == ABIP_NULL)
    {
        abip_free(d);
        printf("ABIPData struct must contain a vector 'b'. \n");
    }

    c_mex = data->c;

    if (c_mex == ABIP_NULL)
    {
        abip_free(d);
        printf("ABIPData struct must contain a vector 'c'.");
    }

    settings = prhs[1];

    d->n = (abip_int)data->length_c;
    d->m = (abip_int)data->length_b;

#if SFLOAT > 0
    d->b = castTo_abip_float_arr(b_mex, d->m);
    d->c = cast_to_abip_float_arr(c_mex, d->n);
#else
    d->b = (abip_float *)b_mex;
    d->c = (abip_float *)c_mex;
#endif

    ABIP(set_default_settings)
    (d);

    d->stgs->max_ipm_iters = (abip_int)settings->max_ipm_iters;
    d->stgs->max_admm_iters = (abip_int)settings->max_admm_iters;
    d->stgs->eps = (abip_float)settings->eps;
    d->stgs->cg_rate = (abip_float)settings->cg_rate;
    d->stgs->alpha = (abip_float)settings->alpha;
    d->stgs->rho_y = (abip_float)settings->rho_y;
    d->stgs->normalize = (abip_int)settings->normalize;
    d->stgs->scale = (abip_float)settings->scale;
    d->stgs->sparsity_ratio = (abip_float)settings->sparsity_ratio;
    d->stgs->adaptive = (abip_int)settings->adaptive;
    d->stgs->adaptive_lookback = (abip_int)settings->adaptive_lookback;
    d->stgs->verbose = (abip_int)settings->verbose;

    A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
    A->n = d->n;
    A->m = d->m;

#if DLONG > 0
    A->p = (abip_int *)A_mex->cInd;
    A->i = (abip_int *)A_mex->rInd;
#else
    A->p = cast_to_abip_int_arr(A_mex->p, A->n + 1);
    A->i = cast_to_abip_int_arr(A_mex->i, A->p[A->n]);
#endif

#if SFLOAT > 0
    A->x = cast_to_abip_float_arr(A_mex, A->p[A->n]);
#else
    A->x = (abip_float *)A_mex;
#endif

    d->A = A;
    d->sp = (abip_float)A->p[A->n] / (A->m * A->n);

    d->stgs->warm_start = parse_warm_start((mxArray *)mxGetField(data, 0, "x"), &(sol.x), d->n);
    d->stgs->warm_start |= parse_warm_start((mxArray *)mxGetField(data, 0, "y"), &(sol.y), d->m);
    d->stgs->warm_start |= parse_warm_start((mxArray *)mxGetField(data, 0, "s"), &(sol.s), d->n);

    status = ABIP(main)(d, &sol, &info);

    set_output_field(&plhs[0], sol.x, d->n);
    set_output_field(&plhs[1], sol.y, d->m);
    set_output_field(&plhs[2], sol.s, d->n);

    plhs[3] = mxCreateStructArray(1, one, num_info_fields, info_fields);

    mxSetField(plhs[3], 0, "status", mxCreateString(info.status));

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "ipm_iter", tmp);
    *mxGetPr(tmp) = (abip_float)info.ipm_iter;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "admm_iter", tmp);
    *mxGetPr(tmp) = (abip_float)info.admm_iter;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "pobj", tmp);
    *mxGetPr(tmp) = info.pobj;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "dobj", tmp);
    *mxGetPr(tmp) = info.dobj;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "resPri", tmp);
    *mxGetPr(tmp) = info.res_pri;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "resDual", tmp);
    *mxGetPr(tmp) = info.res_dual;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "relGap", tmp);
    *mxGetPr(tmp) = info.rel_gap;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "resInfeas", tmp);
    *mxGetPr(tmp) = info.res_infeas;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "resUnbdd", tmp);
    *mxGetPr(tmp) = info.res_unbdd;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "setupTime", tmp);
    *mxGetPr(tmp) = info.setup_time;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[3], 0, "solveTime", tmp);
    *mxGetPr(tmp) = info.solve_time;

    free_mex(d);
    return;
}

void free_mex(ABIPData *d)
{
    if (d)
    {
#if SFLOAT > 0
        if (d->b)
        {
            abip_free(d->b);
        }
        if (d->c)
        {
            abip_free(d->c);
        }
#endif

        if (d->A)
        {
#if !(DLONG > 0)
            if (d->A->p)
            {
                abip_free(d->A->p);
            }
            if (d->A->i)
            {
                abip_free(d->A->i);
            }
#endif

#if SFLOAT > 0
            if (d->A->x)
            {
                abip_free(d->A->x);
            }
#endif

            abip_free(d->A);
        }

        if (d->stgs)
        {
            abip_free(d->stgs);
        }

        abip_free(d);
    }
}
