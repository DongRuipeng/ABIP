void convolve(double *a, int *na, double *b, int *nb, double *ab)
{
    int nab = *na + *nb - 1;
    for (int i = 0; i < nab; i++)
        ab[i] = 0.0;
    for (int i = 0; i < *na; i++)
        for (int j = 0; j < *nb; j++)
            ab[i + j] += a[i] * b[j];
}