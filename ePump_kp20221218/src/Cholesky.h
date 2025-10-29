#ifndef CHHDR
#define CHHDR

#include "VectorDef.h"

/*
 A set of functions used for solving the following:
 
 Given A, a positive definite matrix, calculate
 
 Q = V1^T * A^(-1) * V2,
 
 where V1 and V2 are two column vectors.
 
 The solution is obtained,
 using Cholesky decomposition with Backward Substitution.
 The Cholesky matrix L is the lower triangular matrix which solves A=L*L^T.
 Then Q is given by the dot product of L^(-1)*V1 and L^(-1)*V2, which
 can be obtained by backward substitution.
 
 Note that the routines do not check that the matrix A is of the right
 form.  Beware to anyone who would pass a matrix A that is not positive definite.
 
 Algorithms taken from Numerical Recipes in C.
*/

void Cholesky(const d_mat &A, d_mat &L, int N);
void BackSub(const d_mat &L, const d_vec &V, d_vec &LinvV, int N);
double dot(const d_vec &A, const d_vec &B, int N);

#endif
