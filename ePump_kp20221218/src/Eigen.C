#include "Eigen.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <limits>

double pythag(const double a, const double b);

inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

void Eigen(const d_mat &A, d_mat &z, d_vec &d, int N)
{
// This subroutine finds the eigenvectors and eigenvalues of the
// real, symmetric NxN matrix A[i][j].
// Indices run from 0 to N-1.
// The eigenvectors are output in the columns of z[i][j].
// The eigenvalues are output in d[i].

    
// The first part of this subroutine is based on tred2 from numerical recipes.
// It uses Householder reduction of A[i][j] to put it in tridiagonal form.
// The matrix z[i][j] becomes the orthoganal matrix effecting the transformation,
// d[i] become the diagonal elements of the tridiagonal matrix and e[i] are the
// off-diagonal elements, with e[0]=0.
    
// The second part of the subroutine is based on tqli from numerical recipes.
// It uses the QL algorithm with implicit shifts to find the eigenvalues
// and eigenvectors of the tridiagonal matrix found in part one (and
// hence the eigenvectors of A[i][j]).
    
 
// Part one:
    
    int l,k,j,i;
    double scale,hh,h,g,f;
    d_vec e(N);
    
    for (i=0;i<N;i++) {
        for (j=0;j<i;j++) {
            if (A[i][j]!=A[j][i]) {
                cerr<<"Error: Matrix A is not symmetric.\n"
                <<"Check for mismatch between index range and number of input values."<<endl;
                exit(1);
            }
            z[i][j]=A[i][j];
            z[j][i]=A[j][i];
        }
        z[i][i]=A[i][i];
    }
    
	for (i=N-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += abs(z[i][k]);
			if (scale == 0.0)
				e[i]=z[i][l];
			else {
				for (k=0;k<i;k++) {
					z[i][k] /= scale;
					h += z[i][k]*z[i][k];
				}
				f=z[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z[i][l]=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
                    z[j][i]=z[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z[j][k]*z[i][k];
					for (k=j+1;k<i;k++)
						g += z[k][j]*z[i][k];
					e[j]=g/h;
					f += e[j]*z[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z[j][k] -= (f*e[k]+g*z[i][k]);
				}
			}
		} else
			e[i]=z[i][l];
		d[i]=h;
	}
	d[0]=0.0;
	e[0]=0.0;
	for (i=0;i<N;i++) {
        if (d[i] != 0.0) {
            for (j=0;j<i;j++) {
                g=0.0;
                for (k=0;k<i;k++)
                    g += z[i][k]*z[k][j];
                for (k=0;k<i;k++)
                    z[k][j] -= g*z[k][i];
            }
        }
        d[i]=z[i][i];
        z[i][i]=1.0;
        for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
	}
    
// Part two:
    
    int m,iter;
	double s,r,p,dd,c,b;
	const double EPS=numeric_limits<double>::epsilon();
	for (i=1;i<N;i++) e[i-1]=e[i];
	e[N-1]=0.0;
	for (l=0;l<N;l++) {
		iter=0;
		do {
			for (m=l;m<N-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
                if (abs(e[m]) <= EPS*dd) break;
			}
			if (m != l) {
                if (iter++ == 30) {
                    cerr<<"Too many iterations in tqli"<<endl;
                    exit(1);
                }
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
                    for (k=0;k<N;k++) {
                        f=z[k][i+1];
                        z[k][i+1]=s*z[k][i]+c*f;
                        z[k][i]=c*z[k][i]-s*f;
                    }
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

void eigsrt(d_mat &z, d_vec &d, int N)
{
// Sort the eigenvectors and eigenvalues in order of largest to smallest eigenvalue.
	int k,j,i;
    double p;
	for (i=0;i<N-1;i++) {
		p=d[k=i];
		for (j=i+1;j<N;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
            for (j=0;j<N;j++) {
                p=z[j][i];
				z[j][i]=z[j][k];
				z[j][k]=p;
            }
        }
	}
}

double pythag(const double a, const double b) {
	double absa=abs(a), absb=abs(b);
	return (absa > absb ? absa*sqrt(1.0+absb*absb/(absa*absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+absa*absa/(absb*absb))));
}

