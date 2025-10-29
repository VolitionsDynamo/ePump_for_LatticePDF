#ifndef EGHDR
#define EGHDR

#include "VectorDef.h"

/*
 The subroutine Eigen inds the eigenvectors and eigenvalues of the
 real, symmetric NxN matrix A[i][j].
 Indices run from 0 to N-1.
 The eigenvectors are output in the columns of z[i][j].
 The eigenvalues are output in d[i].

 Based on the algorithms tred2 and tqli from Numerical Recipes in C.
 
 The subroutine eigsrt sorts the eigenvectors and eigenvalues in order 
 of largest to smallest eigenvalue.  It is also based on an algorithm
 from Numerical Recipes in C.
*/

void Eigen(const d_mat &A, d_mat &z, d_vec &d, int N);
void eigsrt(d_mat &z, d_vec &d, int N);


#endif
