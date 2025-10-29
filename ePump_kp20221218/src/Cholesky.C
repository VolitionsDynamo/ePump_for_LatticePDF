#include "Cholesky.h"

void Cholesky(const d_mat &A, d_mat &L, int N)
{
    double sum;
    for (int i=0;i<N;i++) {
        for (int j=i;j<N;j++) {
            sum=A[i][j];
            for (int k=i-1;k>=0;k--) {
                sum-=L[i][k]*L[j][k];
            }
            if (i==j) {
                L[i][i]=sqrt(sum);
            } else {
                L[j][i]=sum/L[i][i];
            }
        }
    }
}

void BackSub(const d_mat &L, const d_vec &V, d_vec &LinvV, int N)
{
    double sum;

    for (int i=0;i<N;i++) {
        sum=V[i];
        for (int k=i-1;k>=0;k--) {
            sum-=L[i][k]*LinvV[k];
        }
        LinvV[i]=sum/L[i][i];
    }
    
}

double dot(const d_vec &A, const d_vec &B, int N)
{
    double sp = 0.0;
    for (int i=0;i<N;i++) {
        sp+=A[i]*B[i];
    }
    return sp;
}

