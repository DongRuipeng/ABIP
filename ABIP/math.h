// The header declare the matrix skeleton 
#ifndef MATH_H_GUARD
#define MATH_H_GUARD 

typedef struct MATRIX
{
    const int m; 
    const int n; 
    double *A; 
} Matrix;

typedef struct VECTOR
{
    const int m; 
    double *a;
} Vector;


// typedef struct PROBLEM_DATA
// {
//     Matrix *A; 
//     double *b; 
//     double *c; 

//     const int m; 
//     const int n; 
// };

// typedef struct WORKER_DATA
// {
//     /* data */
// };

#endif