#include "ePump.h"
#include "Cholesky.h"
#include "Eigen.h"
#include <stdio.h>
#include <sstream>
#include <memory>

double MaxComp(const double &a, const double &b);

void TaylorCoeffs(double &dX, double &ddX, const double &X0, const double &Xp, const double &Xm, const double &Tp, const double &Tm, const double &Tbar);


void ePump::SetQ(int k) {
    
    for (int a=0;a<Na[k];a++) {
        for (int b=0;b<=a;b++) {
            Q[k][a][b]=dot(dX[k][a],dX[k][b],Ni);
            if (a!=b) {Q[k][b][a]=Q[k][a][b];}
        }
    }
    
}

void ePump::SetXdyn(int k) {
    
    Xdyn[k].resize(Na[k]);
    XXdyn[k].resize(Na[k]);
    for (int a=0;a<Na[k];a++) {
        Xdyn[k][a].resize(Ni);
        XXdyn[k][a].resize(Ni);
    }
    
    double dx, ddx, tp, tm, tbar;
    for (int a=0;a<Na[k];a++) {
        for (int i=0;i<Ni;i++) {
            if (T_flag==1) {
                tp=TdynRMS[i];
                tm=TdynRMS[i];
                tbar=TdynRMS[i];
            } else {
                tp=TdynP[i];
                tm=TdynM[i];
                tbar=TdynBar[i];
            }
            TaylorCoeffs(dx,ddx,Xmat[k][a][0],Xmat[k][a][2*i+1],Xmat[k][a][2*i+2],tp,tm,tbar);
            Xdyn[k][a][i]=dx;
            XXdyn[k][a][i]=ddx;
        }
    }
    
}


void ePump::PrepareAMlinear(int k) {
    
    d_mat temp;
    d_vec Atemp;
    d_mat Mtemp;
    
    temp.resize(Ni);
    Atemp.resize(Ni);
    Mtemp.resize(Ni);
    for (int i=0;i<Ni;i++) {
        temp[i].resize(Na[k]);
        Mtemp[i].resize(Ni);
    }
    
    for (int i=0;i<Ni;i++) {
        Atemp[i]=0.0;
        for (int a=0;a<Na[k];a++) {
            temp[i][a]=0.0;
            for (int b=0;b<Na[k];b++) {
                temp[i][a]+=Cm[k][a][b]*Xdyn[k][b][i];
            }
            Atemp[i]+=(XE[k][a]-X0[k][a])*temp[i][a];
        }
        Atemp[i]*=(weight[k]/Tsq);
        A[i]+=Atemp[i];
    }
    
    
    for (int i=0;i<Ni;i++) {
        for (int j=i;j<Ni;j++) {
            Mtemp[i][j]=0.0;
            for (int a=0;a<Na[k];a++) {
                Mtemp[i][j]+=Xdyn[k][a][i]*temp[j][a];
            }
            Mtemp[i][j]*=(weight[k]/Tsq);
            M[i][j]+=Mtemp[i][j];
        }
    }
    
}


void ePump::PrepareAM(int k, double alpha) {
    
    d_mat temp;
    d_vec Atemp;
    d_mat Mtemp;
    d_vec Btemp;
    
    temp.resize(Ni);
    Atemp.resize(Ni);
    Mtemp.resize(Ni);
    Btemp.resize(Ni);
    for (int i=0;i<Ni;i++) {
        temp[i].resize(Na[k]);
        Mtemp[i].resize(Ni);
    }
    
    for (int i=0;i<Ni;i++) {
        Atemp[i]=0.0;
        Btemp[i]=0.0;
        for (int a=0;a<Na[k];a++) {
            temp[i][a]=0.0;
            for (int b=0;b<Na[k];b++) {
                temp[i][a]+=Cm[k][a][b]*(Xdyn[k][b][i]+XXdyn[k][b][i]*z0[i]);
                Btemp[i]+=(XE[k][a]-X0new[k][a])*Cm[k][a][b]*XXdyn[k][b][i];
            }
            Atemp[i]+=(XE[k][a]-X0new[k][a])*temp[i][a];
        }
        Btemp[i]*=alpha;
        Atemp[i]*=(weight[k]/Tsq);
        A[i]+=Atemp[i];
    }
    
    
    for (int i=0;i<Ni;i++) {
        for (int j=i;j<Ni;j++) {
            Mtemp[i][j]=0.0;
            for (int a=0;a<Na[k];a++) {
                Mtemp[i][j]+=(Xdyn[k][a][i]+XXdyn[k][a][i]*z0[i])*temp[j][a];
            }
            if (i==j) {
                Mtemp[i][j]-=Btemp[i];
            }
            Mtemp[i][j]*=(weight[k]/Tsq);
            M[i][j]+=Mtemp[i][j];
        }
    }
    
}

void ePump::PrepareM(int k) {
    
    d_mat temp;
    d_mat Mtemp;
    
    temp.resize(Ni);
    Mtemp.resize(Ni);
    for (int i=0;i<Ni;i++) {
        temp[i].resize(Na[k]);
        Mtemp[i].resize(Ni);
    }
    
    for (int i=0;i<Ni;i++) {
        for (int a=0;a<Na[k];a++) {
            temp[i][a]=0.0;
            for (int b=0;b<Na[k];b++) {
                temp[i][a]+=Cm[k][a][b]*dX[k][b][i];
            }
        }
    }
    
    for (int i=0;i<Ni;i++) {
        for (int j=i;j<Ni;j++) {
            Mtemp[i][j]=0.0;
            for (int a=0;a<Na[k];a++) {
                Mtemp[i][j]+=dX[k][a][i]*temp[j][a];
            }
            M[i][j]+=Mtemp[i][j];
        }
    }
}


void ePump::ConstructUpdateMatrix(void)
{
    z0.resize(Ni);
    z0sqr.resize(Ni);

    d_mat ident(Ni);
    d_mat LinvIdent(Ni);
    d_mat IplusM(Ni);
    
    for (int i=0;i<Ni;i++) {
        ident[i].resize(Ni);
        IplusM[i].resize(Ni);
        LinvIdent[i].resize(Ni);
    }
    for (int i=0;i<Ni;i++) {
        for (int j=0;j<Ni;j++) {
                if (i==j) {
                    ident[i][j]=1.0;
                } else {
                    ident[i][j]=0.0;
                }
        }
    }

    if (!DiagonalQuad) {
        
        cout<<" Only Linear terms included in updating of best-fit predictions."<<endl<<endl;
        ePout<<" Only Linear terms included in updating of best-fit predictions."<<endl<<endl;
    
    for (int i=0;i<Ni;i++) {
        A[i]=0.0;
        for (int j=i;j<Ni;j++) {
            M[i][j]=0.0;
        }
    }
    for (int k; k<Ndata;k++) {
        if (dataIncluded[k]) {
            PrepareAMlinear(k);
        }
    }
    for (int i=0;i<Ni;i++) {
        for (int j=0;j<i;j++) {
            M[i][j]=M[j][i];
            IplusM[i][j]=M[i][j];
            IplusM[j][i]=M[j][i];
        }
        IplusM[i][i]=1.0+M[i][i];
    }
    Cholesky(IplusM,L,Ni);
    BackSub(L,A,LinvA,Ni);
    double dz=0.0;
    for (int i=0;i<Ni;i++) {
        BackSub(L,ident[i],LinvIdent[i],Ni);
        z0[i]=dot(LinvIdent[i],LinvA,Ni);
        z0sqr[i]+=z0[i]*z0[i];
        dz+=z0sqr[i];
    }
        dz=sqrt(dz);
        cout<<" Length of best-fit z0 vector is "<<dz<<endl<<endl;
        ePout<<" Length of best-fit z0 vector is "<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

    for (int k=0;k<Ndata;k++) {
        for (int a=0;a<Na[k];a++) {
            X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
            if (DiagonalQuad) {
                X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
            }
        }
    }

    // New delta Chi^2 for each data set

    Chi2=0.0;
    Chi2new=0.0;
    for (int k=0;k<Ndata;k++) {
        if (dataIncluded[k]) {
            dChi2new[k]=DeltaChiSquare(X0new[k],k);
            Chi2+=weight[k]*dChi2[k];
            Chi2new+=weight[k]*dChi2new[k];
        }
    }
    
    //  New delta Chi^2 for original global analysis:
    
    Chi20new=0.0;
    for (int i=0;i<Ni;i++) {
        Chi20new+=z0sqr[i];
    }
    Chi20new*=Tsq;
    Chi2new+=Chi20new;
    
    cout<<" Total Chi2 (z=0)  = "<<fixed<<Chi2<<endl;
    cout<<" Total Chi2 (z=z0) = "<<fixed<<Chi2new<<endl<<endl;
    ePout<<" Total Chi2 (z=0)  = "<<fixed<<Chi2<<endl;
    ePout<<" Total Chi2 (z=z0) = "<<fixed<<Chi2new<<endl<<endl;
   
    
    } else {
        
        cout<<" Diagonal Quadratic terms included in all updated best-fit predictions."<<endl<<endl;
        ePout<<" Diagonal Quadratic terms included in all updated best-fit predictions."<<endl<<endl;
        
        d_vec eps(Ni);   // vector z0(n)-z0(n-1)
        double dz; // length of vector eps
        double alpha1=1.0;  // acceleration parameter for vector eps
        double alpha2=0.0;   // reduction parameter for non-positive-definite term in M
        
//    for (int n=1;n<4;n++) {
    
        for (int i=0;i<Ni;i++) {
            z0[i]=0.0;
            z0sqr[i]=0.0;
        }
        
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=0)       = "<<fixed<<Chi2<<endl<<endl;
        ePout<<" Total Chi2 (z=0)       = "<<fixed<<Chi2<<endl<<endl;
   

        int n=0;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<
        scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        n=1;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

        n=2;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=3;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=4;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=5;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=6;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=7;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=8;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

        alpha2=1.0;
    
        n=9;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
        n=10;
        
        for (int i=0;i<Ni;i++) {
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,alpha2);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
            A[i]-=z0[i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dz=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            eps[i]=alpha1*dot(LinvIdent[i],LinvA,Ni);
            dz+=eps[i]*eps[i];
            z0[i]+=eps[i];
            z0sqr[i]=z0[i]*z0[i];
        }
        dz=sqrt(dz);
            cout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
            ePout<<" Iteration "<<n+1<<":  Distance between z0(n="<<n<<") and z0(n="<<n+1<<") is "<<scientific<<dz<<endl<<endl;
    
    // Construct new best-fit predictions for each data point

        for (int k=0;k<Ndata;k++) {
            for (int a=0;a<Na[k];a++) {
                X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
                if (DiagonalQuad) {
                    X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
                }
            }
        }

    // New delta Chi^2 for each data set

        Chi2=0.0;
        Chi2new=0.0;
        for (int k=0;k<Ndata;k++) {
            if (dataIncluded[k]) {
                dChi2new[k]=DeltaChiSquare(X0new[k],k);
                Chi2+=weight[k]*dChi2[k];
                Chi2new+=weight[k]*dChi2new[k];
            }
        }
    
    //  New delta Chi^2 for original global analysis:
    
        Chi20new=0.0;
        for (int i=0;i<Ni;i++) {
            Chi20new+=z0sqr[i];
        }
        Chi20new*=Tsq;
        Chi2new+=Chi20new;
    
        cout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
        ePout<<" Total Chi2 (z=z0(n="<<n+1<<")) = "<<fixed<<Chi2new<<endl<<endl;
   
        

    
   }
    
//    }
        

    
    
    
    
    
    
        
 /*
        
        
        
        
        
        for (int i=0;i<Ni;i++) {
//            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
    }
        /*
   
     n=2;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
    
       n=3;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
  
       n=4;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
  
       n=5;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
  
       n=6;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
  
       n=7;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
  
       n=8;
        for (int i=0;i<Ni;i++) {
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
            A[i]=0.0;
            for (int j=i;j<Ni;j++) {
                M[i][j]=0.0;
            }
        }
        for (int k; k<Ndata;k++) {
            if (dataIncluded[k]) {
                PrepareAM(k,z0iter,z0sqriter);
            }
        }
        for (int i=0;i<Ni;i++) {
            for (int j=0;j<i;j++) {
                M[i][j]=M[j][i];
                IplusM[i][j]=M[i][j];
                IplusM[j][i]=M[j][i];
            }
            IplusM[i][i]=1.0+M[i][i];
        }
        Cholesky(IplusM,L,Ni);
        BackSub(L,A,LinvA,Ni);
        dziter=0.0;
        for (int i=0;i<Ni;i++) {
            BackSub(L,ident[i],LinvIdent[i],Ni);
            z0[i]=dot(LinvIdent[i],LinvA,Ni);
            z0sqr[i]=z0[i]*z0[i];
            dziter+=(z0[i]-z0iter[i])*(z0[i]-z0iter[i]);
            z0iter[i]=z0[i];
            z0sqriter[i]=z0sqr[i];
            ePout<<"   n  "<<n+1<<"   "<<z0iter[i]<<endl;
        }
        dziter=sqrt(dziter);
        cout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
        ePout<<" Iteration "<<n+1<<":  Distance between z0["<<n+1<<"] and z0["<<n<<"] is "<<dziter<<endl;
   
  
    
  */
    
    
    
    //
    
    
    Meig();
    
    
    
    
    cout<<"  Vector z0[i] from i=1,Ni  (original eigenvector basis)"<<endl;
    ePout<<"  Vector z0[i] from i=1,Ni  (original eigenvector basis)"<<endl;
    for (int i=0;i<Ni;i++) {
        cout<<"  "<<z0[i]<<endl;
        ePout<<"  "<<z0[i]<<endl;
    }
    cout<<endl;
    ePout<<endl;
    cout<<"  Vector z0[r] = sum(U[r][i]z0[i] over i=1,Ni) from r=1,Ni  (updated eigenvector basis)"<<endl;
    ePout<<"  Vector z0[r] = sum(U[r][i]z0[i] over i=1,Ni) from r=1,Ni  (updated eigenvector basis)"<<endl;
    for (int i=0;i<Ni;i++) {
        double z0r=0.0;
        for (int j=0;j<Ni;j++) {
            z0r+=z[j][i]*z0[j];
        }
        cout<<"  "<<z0r<<endl;
        ePout<<"  "<<z0r<<endl;
    }
    cout<<endl;
    ePout<<endl;
    
    
    
    // Construct new best-fit predictions for each data point

    for (int k=0;k<Ndata;k++) {
        for (int a=0;a<Na[k];a++) {
            X0new[k][a]=X0[k][a]+dot(Xdyn[k][a],z0,Ni);
            if (DiagonalQuad) {
                X0new[k][a]+=(dot(XXdyn[k][a],z0sqr,Ni)/2.0);
            }
        }
    }

    // New delta Chi^2 for each data set
    
    for (int k=0;k<Ndata;k++) {
        dChi2new[k]=DeltaChiSquare(X0new[k],k);
        cout<<"chi2new[k] "<<dChi2new[k]<<endl;
    }
    
    //  New delta Chi^2 for original global analysis:
    
    Chi20new=0.0;
    for (int i=0;i<Ni;i++) {
        Chi20new+=z0sqr[i];
    }
    Chi20new*=Tsq;
    
    // Drop in total Chi^2 after update:
    
//    Chi2diff=-Tsq*dot(LinvA,LinvA,Ni);
//        cout<<"Chi20new "<<Chi20new<<endl;
//        cout<<"Chi2diff "<<Chi2diff<<endl;

}


void ePump::ConstructOptimizeMatrix(void)
{
    for (int i=0;i<Ni;i++) {
        for (int j=i;j<Ni;j++) {
            M[i][j]=0.0;
        }
    }
    for (int k; k<Ndata;k++) {
        PrepareM(k);
    }
    for (int i=0;i<Ni;i++) {
        for (int j=0;j<i;j++) {
            M[i][j]=M[j][i];
        }
    }
    Meig();
}

double ePump::DeltaChiSquare(const d_vec &YY, int k)
                 // Returns dChi^2 for data set k
{
    if (!dataIncluded[k]) {
        return 0.0;
    } else {
        d_vec dxa(Na[k]);
        double dc2=0.0;
        for (int a=0;a<Na[k];a++) {
            dxa[a]=YY[a]-XE[k][a];
        }
        for (int a=0;a<Na[k];a++) {
            dc2+=Cm[k][a][a]*dxa[a]*dxa[a];
            for (int b=a+1;b<Na[k];b++) {
                dc2+=2.0*dxa[a]*Cm[k][a][b]*dxa[b];
            }
        }
        return dc2;
    }
}


void ePump::UpdateErrors0(const d_vec &ddY, double &DY0, double &DYnew)
                 // Give error vector dY as input to calculate DY0 and DYnew.
{
    d_vec ddYr(Ni);
    DY0=sqrt(dot(ddY,ddY,Ni));
    for (int i=0;i<Ni;i++) {
        ddYr[i]=0.0;
        for (int j=0;j<Ni;j++) {
            ddYr[i]+=ddY[j]*z[j][i];
        }
    }
    double DYnewSqr=0.0;
    for (int i=0;i<Ni;i++) {
            DYnewSqr+=ddYr[i]*ddYr[i]/(1.0+d[i]);
    }
    DYnew=sqrt(DYnewSqr);
}

void ePump::UpdateErrors(const d_vec &YY, double &DY0sym, double &DY0up, double &DY0down,
                   double &DYnewSym, double &DYnewUp, double &DYnewDown, double &Y0new)
                 // Give 2*Ni+1 dimensional vector of Y values as input, with
                 //       YY[0]=Y[0] (central value), YY[2*i+1]=Y[i+], YY[2*i+2]=Y[i-].
                 // Calculates Ynew[0], and both symmetric and asymmetric versions of the errors:
                 //      DY0sym, DY0up, DY0down,band DYnewSym, DYnewUp, DYnewDown.
{
//    double dum;
    d_vec dY(Ni);
    d_vec ddY(Ni);
    double dy, ddy, tp, tm, tbar;

    for (int k=0;k<Ni;k++) {
        if (T_flag==1) {
            tp=TdynRMS[k];
            tm=TdynRMS[k];
            tbar=TdynRMS[k];
        } else {
            tp=TdynP[k];
            tm=TdynM[k];
            tbar=TdynBar[k];
        }
        TaylorCoeffs(dy,ddy,YY[0],YY[2*k+1],YY[2*k+2],tp,tm,tbar);
        dY[k]=dy;
        ddY[k]=ddy;
    }
    Y0new=YY[0]+dot(dY,z0,Ni);
    if (DiagonalQuad) {
        Y0new+=(dot(ddY,z0sqr,Ni)/2.0);
    }
    
    for (int k=0;k<Ni;k++) {
        dY[k]=(YY[2*k+1]-YY[2*k+2])/2.0;
    }
    UpdateErrors0(dY,DY0sym,DYnewSym);
    for (int k=0;k<Ni;k++) {
        dY[k]=MaxComp(YY[2*k+1]-YY[0],YY[0]-YY[2*k+2]);
    }
    UpdateErrors0(dY,DY0up,DYnewUp);
    for (int k=0;k<Ni;k++) {
        dY[k]=MaxComp(YY[0]-YY[2*k+2],YY[2*k+1]-YY[0]);
    }
    UpdateErrors0(dY,DY0down,DYnewDown);
}

void ePump::UpdateCorrCosine0(const d_vec &dY1, const d_vec &dY2, double &DY1, double &DY2,  double &Cos12, double &DY1new, double &DY2new, double &Cos12new)
                 // Give error vectors dY1 and dY2 as input.
{
    d_vec dY1r(Ni);
    d_vec dY2r(Ni);
    for (int i=0;i<Ni;i++) {
        dY1r[i]=0.0;
        dY2r[i]=0.0;
        for (int j=0;j<Ni;j++) {
            dY1r[i]+=dY1[j]*z[j][i];
            dY2r[i]+=dY2[j]*z[j][i];
        }
    }
    double DY1newSqr=0.0;
    double DY2newSqr=0.0;
    double DY12new=0.0;
    for (int i=0;i<Ni;i++) {
            DY1newSqr+=dY1r[i]*dY1r[i]/(1.0+d[i]);
            DY2newSqr+=dY2r[i]*dY2r[i]/(1.0+d[i]);
            DY12new+=dY1r[i]*dY2r[i]/(1.0+d[i]);
    }
    DY1new=sqrt(DY1newSqr);
    DY2new=sqrt(DY2newSqr);
    Cos12new=DY12new/DY1new/DY2new;

    DY1=sqrt(dot(dY1,dY1,Ni));
    DY2=sqrt(dot(dY2,dY2,Ni));
    Cos12=dot(dY1,dY2,Ni)/DY1/DY2;
}

void ePump::UpdateCorrCosine(const d_vec &YY1, const d_vec &YY2, double &DY1, double &DY2,  double &Cos12, double &DY1new, double &DY2new, double &Cos12new)
                 // Give 2*Ni+1 dimensional vector of Y1 and Y2 values as input, with
                 //       YY1[0]=Y1[0] (central value), YY1[2*i+1]=Y[i+], YY1[2*i+2]=Y[i-].   (Similarly for Y2.)
                 // Uses symmetric error vectors.
{
    d_vec dY1(Ni);
    d_vec dY2(Ni);
    for (int k=0;k<Ni;k++) {
        dY1[k]=(YY1[2*k+1]-YY1[2*k+2])/2.0;
        dY2[k]=(YY2[2*k+1]-YY2[2*k+2])/2.0;
    }
    UpdateCorrCosine0(dY1,dY2,DY1,DY2,Cos12,DY1new,DY2new,Cos12new);
}

void ePump::ConstructCm1(int k)
{
    if (Nlam[k]==0){
        if (verb) {
            cout<<"Inverse Covariance Matrix:"<<endl;
            ePout<<"Inverse Covariance Matrix:"<<endl;
        }
        for (int a=0;a<Na[k];a++) {
            if (verb) {
                cout<<" a="<<a+1<<", b="<<a+1<<","<<Na[k]<<endl;
                ePout<<" a="<<a+1<<", b="<<a+1<<","<<Na[k]<<endl;
            }
            for (int b=a;b<Na[k];b++) {
                if (b==a) {
                    Cm[k][a][b]=1.0/s[k][a]/s[k][a];
                } else {
                    Cm[k][a][b]=0.0;
                }
                if (b!=a) {
                    Cm[k][b][a]=Cm[k][a][b];
                }
                if (verb) {
                    cout<<Cm[k][a][b]<<"   ";
                    ePout<<Cm[k][a][b]<<"   ";
                }
            }
            if (verb) {
                cout<<endl;
                ePout<<endl;
            }
        }
        if (verb) {
            cout<<endl;
            ePout<<endl;
        }
    } else {
        d_mat calA(Nlam[k]);
        d_mat calL(Nlam[k]);
        d_mat ident(Nlam[k]);
        d_mat calLinvBeta(Na[k]);
        d_mat calLinvIdent(Nlam[k]);
        for (int j;j<Nlam[k];j++) {
            calA[j].resize(Nlam[k]);
            calL[j].resize(Nlam[k]);
            ident[j].resize(Nlam[k]);
            calLinvIdent[j].resize(Nlam[k]);
        }
        for (int a=0;a<Na[k];a++) {
            calLinvBeta[a].resize(Nlam[k]);
        }

        for (int i=0;i<Nlam[k];i++) {
            for (int j=i;j<Nlam[k];j++) {
                if (i==j) {
                    ident[i][j]=1.0;
                    calA[i][j]=1.0;
                } else {
                    ident[i][j]=0.0;
                    ident[j][i]=0.0;
                    calA[i][j]=0.0;
                }
                for (int a=0;a<Na[k];a++) {
                    calA[i][j]+=beta[k][a][i]*beta[k][a][j]/s[k][a]/s[k][a];
                }
                if (i!=j) {
                    calA[j][i]=calA[i][j];
                    
                }
            }
        }
    
        // Calculate inverse covariance matrix:
    
        Cholesky(calA,calL,Nlam[k]);
    
        for (int a=0;a<Na[k];a++) {
            BackSub(calL,beta[k][a],calLinvBeta[a],Nlam[k]);
        }
    
        if (verb) {
            cout<<"Inverse Covariance Matrix:"<<endl;
            ePout<<"Inverse Covariance Matrix:"<<endl;
        }
        for (int a=0;a<Na[k];a++) {
            if (verb) {
                cout<<" a="<<a+1<<", b="<<a+1<<","<<Na[k]<<endl;
                ePout<<" a="<<a+1<<", b="<<a+1<<","<<Na[k]<<endl;
            }
            for (int b=a;b<Na[k];b++) {
                if (b==a) {
                    Cm[k][a][b]=1.0/s[k][a]/s[k][a];
                } else {
                    Cm[k][a][b]=0.0;
                }
                Cm[k][a][b]-=dot(calLinvBeta[a],calLinvBeta[b],Nlam[k])/s[k][a]/s[k][a]/s[k][b]/s[k][b];
                if (b!=a) {
                    Cm[k][b][a]=Cm[k][a][b];
                }
                if (verb) {
                    cout<<Cm[k][a][b]<<"   ";
                    ePout<<Cm[k][a][b]<<"   ";
                }
            }
            if (verb) {
                cout<<endl;
                ePout<<endl;
            }
        }
        if (verb) {
            cout<<endl;
            ePout<<endl;
        }
    
        // Calculate optimal systematic error shifts;
    
        for (int i=0;i<Nlam[k];i++) {
            BackSub(calL,ident[i],calLinvIdent[i],Nlam[k]);
            for (int a=0;a<Na[k];a++) {
                calAinvBeta[k][a][i]=dot(calLinvIdent[i],calLinvBeta[a],Nlam[k]);
            }
        }
    
        cout<<" Nuisance parameters, residuals, etc from original best-fit:"<<endl<<endl;
        ePout<<" Nuisance parameters, residuals, etc from original best-fit:"<<endl<<endl;

        double dchi2=NuisanceParameters(X0[k],k);
        dchi2+=Chi2Residuals(X0[k],k);
        cout<<"  Delta Chi^2 for this data set :   "<<dchi2<<endl;
        ePout<<"  Delta Chi^2 for this data set :   "<<dchi2<<endl;
    }
    
}

double ePump::Chi2Residuals(const d_vec &YY, int k)
{
    // Calculates Shifted Data, residuals, and reduced Chi-square for
    // each data point in data set k, for theory values given by YY[k].
    // (This is only correct if error_type=1.)
    // Returns total reduced Chi-square from the data set.
    //
    
    if (error_type[k]!=1) {
        cerr<<"Residual Calculations are only valid for error_type=1"<<endl;
        cerr<<"ePump::Chi2Residuals was called for data set ["<<k<<"] with error_type="<<error_type[k]<<"."<<endl;
        exit(1);
    }

    d_vec r(Na[k]);  //Residuals
    double dchi2=0.0;
    
    for (int a=0;a<Na[k];a++) {
        r[a]=0.0;
        for (int b=0;b<Na[k];b++) {
        r[a]+=Cm[k][a][b]*(XE[k][b]-YY[b]);
        }
        r[a]*=s[k][a];
    }
    
    cout<<setw(2)<<" a      Data        ShiftedData     residual     reducedChi2"<<endl;
    ePout<<setw(2)<<" a      Data        ShiftedData     residual     reducedChi2"<<endl;
    for (int a=0;a<Na[k];a++) {
        double rc2;
        if (r[a]<0){
            rc2=r[a]*r[a];
        } else {
            rc2=-r[a]*r[a];
        }
        dchi2+=r[a]*r[a];
        cout<<setw(2)<<a+1<<scientific<<setprecision(4)<<setw(13)<<XE[k][a]<<setw(15)<<YY[a]+r[a]*s[k][a]<<setw(15)<<r[a]<<setw(12)<<fixed<<setprecision(2)<<rc2<<endl;
        ePout<<setw(2)<<a+1<<scientific<<setprecision(4)<<setw(13)<<XE[k][a]<<setw(15)<<YY[a]+r[a]*s[k][a]<<setw(15)<<r[a]<<setw(12)<<fixed<<setprecision(2)<<rc2<<endl;
    }
    cout<<endl;
    ePout<<endl;
    return dchi2;
}

double ePump::NuisanceParameters(const d_vec &YY, int k)
{
    // Calculates optimal Nuisance Parameters for
    // data set k, using theory values given by YY[k].
    // (This is only correct if error_type=1.)
    // Returns contribution to Chi-square these nuisance parameters.
    //
    
    if (error_type[k]!=1) {
        cerr<<"Nuisance Parameters calculation are only set up for error_type=1"<<endl;
        cerr<<"ePump::NuisanceParameters was called for data set ["<<k<<"] with error_type="<<error_type[k]<<"."<<endl;
        exit(1);
    }

    d_vec lam(Nlam[k]);   // Optimal systematic error shifts
    double l2=0.0;
    for (int i=0;i<Nlam[k];i++) {
        lam[i]=0.0;
        for (int a=0;a<Na[k];a++) {
            lam[i]+=calAinvBeta[k][a][i]*(XE[k][a]-YY[a])/s[k][a]/s[k][a];
        }
        l2+=lam[i]*lam[i];
    }
    
    cout<<" Sum(l[i]^2)  : ";
    ePout<<" Sum(l[i]^2)  : ";
    for (int i=0;i<Nlam[k];i++) {
        cout<<"      l["<<setfill('0')<<setw(2)<<i+1<<"]   ";
        ePout<<"      l["<<setfill('0')<<setw(2)<<i+1<<"]   ";
    }
    cout<<setfill(' ')<<endl;
    ePout<<setfill(' ')<<endl;

    cout<<scientific<<setprecision(4)<<setw(12)<<l2<<"  :";
    ePout<<scientific<<setprecision(4)<<setw(12)<<l2<<"  :";
    for (int i=0;i<Nlam[k];i++) {
        cout<<scientific<<setprecision(4)<<setw(14)<<lam[i];
        ePout<<scientific<<setprecision(4)<<setw(14)<<lam[i];
    }
    cout<<endl;
    ePout<<endl;

    return l2;
}


void ePump::ConstructCm2(int k)
{
    d_mat cov(Na[k]);
    d_mat covL(Na[k]);
    d_mat ident(Na[k]);
    d_mat covLinvIdent(Na[k]);
    for (int a;a<Na[k];a++) {
        cov[a].resize(Na[k]);
        covL[a].resize(Na[k]);
        ident[a].resize(Na[k]);
        covLinvIdent[a].resize(Na[k]);
    }
    for (int a=0;a<Na[k];a++) {
        for (int b=a;b<Na[k];b++) {
            cov[a][b]=s[k][a]*rhocc[k][a][b]*s[k][b];
            for (int i=0;i<Nlam[k];i++) {
                cov[a][b]+=beta[k][a][i]*beta[k][b][i];
            }
            if (a==b) {
                ident[a][b]=1.0;
            } else {
                ident[a][b]=0.0;
                ident[b][a]=0.0;
                cov[b][a]=cov[a][b];
            }
        }
    }
    
    Cholesky(cov,covL,Na[k]);
    for (int a=0;a<Na[k];a++) {
        BackSub(covL,ident[a],covLinvIdent[a],Na[k]);
    }

    if (verb) {
        cout<<"Inverse Covariance Matrix:"<<endl;
        ePout<<"Inverse Covariance Matrix:"<<endl;
    }
    for (int a=0;a<Na[k];a++) {
        if (verb) {
            cout<<" a="<<a+1<<", b="<<a+1<<","<<Na[k]<<endl;
            ePout<<" a="<<a+1<<", b="<<a+1<<","<<Na[k]<<endl;
        }
        for (int b=a;b<Na[k];b++) {
            Cm[k][a][b]=dot(covLinvIdent[a],covLinvIdent[b],Na[k]);
            if (a!=b) {
                Cm[k][b][a]=Cm[k][a][b];
            }
            if (verb) {
                cout<<Cm[k][a][b]<<"   ";
                ePout<<Cm[k][a][b]<<"   ";
            }
        }
        if (verb) {
            cout<<endl;
            ePout<<endl;
        }
    }
    if (verb) {
        cout<<endl;
        ePout<<endl;
    }
    
}


void ePump::Meig(void)
{
    int i,j;
    d_mat A(Ni), X(Ni);
    d_vec n(Ni);
    for (i=0;i<Ni;i++) {
        z.resize(Ni);
        d.resize(Ni);
    }
    for (i=0;i<Ni;i++) {
        z[i].resize(Ni);
        X[i].resize(Ni);
        A[i].resize(Ni);
    }
    for (i=0;i<Ni;i++) {
        for (j=0;j<Ni;j++) {
            A[i][j]=M[i][j];
        }
    }
    
    if (verb) {
        cout<<" Input Matrix M:"<<endl<<endl;
        ePout<<" Input Matrix M:"<<endl<<endl;
        for (i=0;i<Ni;i++) {
            for (j=0;j<Ni;j++) {
                cout<<setw(12)<<fixed<<setprecision(9);
                ePout<<setw(12)<<fixed<<setprecision(9);
                cout<<A[i][j]<<"  ";
                ePout<<A[i][j]<<"  ";
            }
            cout<<endl;
            ePout<<endl;
        }
        cout<<endl<<endl;
        ePout<<endl<<endl;
    }
    Eigen(A,z,d,Ni);

    eigsrt(z,d,Ni);

//  Choose sign of eigenvector so it is predominantly positive:
    if (df_flag==0) {
       double zsum;
       for (i=0;i<Ni;i++) {
           zsum=0.0;
           for (j=0;j<Ni;j++) {
               zsum+=z[j][i];
           }
           if (zsum<0.0) {
               if (verb) {
                   cout<<" Changed sign of eigenvector "<<i<<endl;
                   ePout<<" Changed sign of eigenvector "<<i<<endl;
               }
               for (j=0;j<Ni;j++) {
                   z[j][i]=-z[j][i];
               }
           }
       }
    }

    cout<<endl;
    ePout<<endl;
    
    if (verb) {
        cout<<" Sorted (column) eigenvectors of M:"<<endl<<endl;
        ePout<<" Sorted (column) eigenvectors of M:"<<endl<<endl;
        for (i=0;i<Ni;i++) {
            for (j=0;j<Ni;j++) {
                cout<<setw(12)<<fixed<<setprecision(9);
                ePout<<setw(12)<<fixed<<setprecision(9);
                cout<<z[i][j]<<"  ";
                ePout<<z[i][j]<<"  ";
            }
            cout<<endl;
            ePout<<endl;
        }
        cout<<endl;
        ePout<<endl;
    }
    
    cout<<"-----------------------------------"<<endl;
    cout<<"     Sorted eigenvalues of M:"<<endl;
    cout<<"-----------------------------------"<<endl<<endl;
    ePout<<"-----------------------------------"<<endl;
    ePout<<"Sorted eigenvalues of M:"<<endl;
    ePout<<"-----------------------------------"<<endl<<endl;
    for (i=0;i<Ni;i++) {
        cout<<setw(12)<<fixed<<setprecision(8)<<right;
        ePout<<setw(12)<<fixed<<setprecision(8)<<right;
        cout<<d[i];
        ePout<<d[i];
    }
    cout<<endl<<endl;
    ePout<<endl<<endl;
    
    if (verb) {
        cout<<" Check normalization of eigenvectors:"<<endl<<endl;
        ePout<<" Check normalization of eigenvectors:"<<endl<<endl;
        for (j=0;j<Ni;j++) {
            n[j]=0.0;
            for (i=0;i<Ni;i++) {
                n[j]+=z[i][j]*z[i][j];
            }
            cout<<"e["<<j<<"].e["<<j<<"] = ";
            ePout<<"e["<<j<<"].e["<<j<<"] = ";
            cout<<setw(19)<<fixed<<setprecision(16)<<n[j]<<endl;
            ePout<<setw(19)<<fixed<<setprecision(16)<<n[j]<<endl;
        }
        cout<<endl<<endl;
        ePout<<endl<<endl;
        cout<<" Check that they are eigenvectors:"<<endl;
        cout<<" Calculate A.O-O.D, where O is column matrix of eigenvectors"<<endl;
        cout<<" and D is diagonal matrix of eigenvalues"<<endl<<endl;
        ePout<<" Check that they are eigenvectors:"<<endl;
        ePout<<" Calculate A.O-O.D, where O is column matrix of eigenvectors"<<endl;
        ePout<<" and D is diagonal matrix of eigenvalues"<<endl<<endl;
        for (i=0;i<Ni;i++) {
            for (j=0;j<Ni;j++) {
                X[i][j]=-z[i][j]*d[j];
                for (int k=0;k<Ni;k++) {
                    X[i][j]+=A[i][k]*z[k][j];
                }
            }
        }
        for (i=0;i<Ni;i++) {
            for (j=0;j<Ni;j++) {
                cout<<setw(10)<<scientific<<setprecision(2)<<X[i][j]<<"  ";
                ePout<<setw(10)<<scientific<<setprecision(2)<<X[i][j]<<"  ";
            }
            cout<<endl;
            ePout<<endl;
        }
        cout<<endl;
        ePout<<endl;
    }
    
//    for (i=0;i<Ni;i++) {
//        if (d[i]<-0.5) {
//            cout<<"Eigenvalue "<<i+1<<" is less than -0.5.  It has been replaced by -0.5."<<endl;
//            ePout<<"Eigenvalue "<<i+1<<" is less than -0.5.  It has been replaced by -0.5."<<endl;
//            d[i]=-0.5;
//        }
//        cout<<d[i];
//        ePout<<d[i];
//    }
//    cout<<endl<<endl;
//    ePout<<endl<<endl;
    
}

void ePump::OptimizedEVContributions(void)
{
    
    d_vec v2sum(Ni);
    
    for (int i=0;i<Ni;i++) {
        v2sum[i]=0.0;
    }
    
    cout<<"--------------------------------------------------------------------------------------"<<endl;
    cout<<" Fractional contribution of Optimized EV i (column) to variance of observable a (row) "<<endl;
    cout<<"--------------------------------------------------------------------------------------"<<endl<<endl;;
    ePout<<"--------------------------------------------------------------------------------------"<<endl;
    ePout<<" Fractional contribution of Optimized EV i (column) to variance of observable a (row) "<<endl;
    ePout<<"--------------------------------------------------------------------------------------"<<endl<<endl;;
    
    cout<<" DataSet Obs";
    ePout<<" DataSet Obs";
    for (int i=1;i<Ni+1;i++) {
        cout<<"        EV"<<setfill('0')<<setw(2)<<i;
        ePout<<"        EV"<<setfill('0')<<setw(2)<<i;
    }
    cout<<setfill(' ')<<endl;
    ePout<<setfill(' ')<<endl;
    
    for (int i=0;i<Ndata;i++) {
        for (int a=0;a<Na[i];a++) {
            cout<<setw(5)<<i+1<<setw(6)<<a+1<<"   ";
            ePout<<setw(5)<<i+1<<setw(6)<<a+1<<"   ";
            for (int r=0;r<Ni;r++) {
                double dxr=0.0;
                for (int j=0;j<Ni;j++) {
                    dxr+=dX[i][a][j]*z[j][r];
                }
                double v2=dxr*dxr/Q[i][a][a];
                cout<<setw(12)<<fixed<<setprecision(8)<<right;
                ePout<<setw(12)<<fixed<<setprecision(8)<<right;
                cout<<v2;
                ePout<<v2;
                v2sum[r]+=v2;
            }
            cout<<endl;
            ePout<<endl;
        }
    }
    cout<<endl;
    ePout<<endl;
    cout<<" Totals       ";
    ePout<<" Totals       ";
    for (int r=0;r<Ni;r++) {
            cout<<setw(12)<<fixed<<setprecision(8)<<right;
            ePout<<setw(12)<<fixed<<setprecision(8)<<right;
        cout<<v2sum[r];
        ePout<<v2sum[r];
    }
    cout<<endl;
    ePout<<endl;
    cout<<" Eigenvalues  ";
    ePout<<" Eigenvalues  ";
    for (int r=0;r<Ni;r++) {
            cout<<setw(12)<<fixed<<setprecision(8)<<right;
            ePout<<setw(12)<<fixed<<setprecision(8)<<right;
        cout<<d[r];
        ePout<<d[r];
    }
    cout<<endl;
    ePout<<endl;
    cout<<" (Totals and Eigenvalues should be equal.)    ";
    ePout<<" (Totals and Eigenvalues should be equal.)    ";
    cout<<endl<<endl;
    ePout<<endl<<endl;
}

void ePump::OptimizePDFs_CTEQ(void)
{
    stringstream line;
    vector<shared_ptr<ifstream>> infile;
    vector<shared_ptr<ofstream>> outfile;
    char infilename[80];
    char outfilename[80];
    char outfileFit[80];
    char num[3];
    double dumD;
    int dumI,nx,nt,nfmx,nfval,nx0,nt0,nfmx0,nfval0;
    int npts;
    string inputstring;
    d_vec UPD(2*Ni+1);

    cout<<"---------------------------------------"<<endl;
    cout<<"     Optimizing PDFs, CTEQ format "<<endl;
    cout<<"---------------------------------------"<<endl;
    ePout<<"---------------------------------------"<<endl;
    ePout<<"     Optimizing PDFs, CTEQ format "<<endl;
    ePout<<"---------------------------------------"<<endl;

    for (int k=0;k<Ni*2+1;k++) {
// Open the input files
        strcpy(infilename,PDFinfile);
        strcat(infilename,".");
        sprintf(num,"%02d",k);
        strcat(infilename,num);
        strcat(infilename,".pds");
        infile.push_back( make_shared<ifstream>(infilename));
        if ( !infile[k]->is_open() ) {
            cerr<<"Error opening input file: "<<infilename<<endl;
            exit(1);
        }
        cout<<infilename<<endl;
        ePout<<infilename<<endl;
// Open the output files
        strcpy(outfilename,PDFoutfile);
        strcat(outfilename,".");
        strcat(outfilename,num);
        strcpy(outfileFit,outfilename);
        strcat(outfilename,".pds");
        outfile.push_back( make_shared<ofstream>(outfilename));
        if ( !outfile[k]->is_open() ) {
            cerr<<"Error opening output file: "<<outfilename<<endl;
            exit(1);
        }
        cout<<outfilename<<endl;
        ePout<<outfilename<<endl;
        
// Read in headers from .pds file
        getline(*infile[k],inputstring);
        *outfile[k]<<"PDF Table for Fit #: "<<outfileFit<<endl;
        for (int i=0;i<4;i++) {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        }
        line.str(inputstring);

        line>>dumI>>dumD>>dumI>>dumI>>dumI>>nfmx>>nfval;
        line.clear();
        for (int i=0;i<2;i++) {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        }
        line.str(inputstring);
        line>>nx>>nt;
        line.clear();

        if (k==0) {
            nx0=nx;
            nt0=nt;
            nfmx0=nfmx;
            nfval0=nfval;
            npts=(nx+1)*(nt+1)*(nfval+nfmx+1);
// Assume no photon.  Must modify this if photon included.
        } else {
            if (nx0!=nx||nt0!=nt||nfmx0!=nfmx||nfval0!=nfval) {
                cerr<<"Error: nx,nt,nfmx,nval values do not match in file "<<infilename<<endl;
                exit(1);
            }
        }

        do {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        } while (inputstring!="Parton Distribution Table:");
    }


// Read out rotated table files to new .pds files:
    
    double df1,df2,f1,f2;
    for (int i=0;i<npts;i++) {
        for (int k=0;k<Ni*2+1;k++) {
            *infile[k]>>UPD[k];
            if ((i+1)%6==0) {
                getline(*infile[k],inputstring);
            }
        }
// Central .pds file
        *outfile[0]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<UPD[0];
        if ((i+1)%6==0) {
            *outfile[0]<<endl;
        }
// Rotated Eigenset .pds files
        for (int k=0;k<Ni;k++) {
            f1=f2=UPD[0];
            for (int j=0;j<Ni;j++) {
                if (df_flag==0||df_flag==2) {
                    df1=UPD[2*j+1]-UPD[0];
                    df2=UPD[0]-UPD[2*j+2];
                } else if (df_flag==1) {
                    df1=0.5*(UPD[2*j+1]-UPD[2*j+2]);
                    df2=df1;
                } else if (df_flag==3) {
                    if (z[j][k]>0.0) {
                        df1=UPD[2*j+1]-UPD[0];
                        df2=UPD[0]-UPD[2*j+2];
                    } else {
                        df2=UPD[2*j+1]-UPD[0];
                        df1=UPD[0]-UPD[2*j+2]; 
                    }
                } else if (df_flag==4) {
                    double df=0.5*(UPD[2*j+1]-UPD[2*j+2]);
                    double R=0.5*(UPD[2*j+1]+UPD[2*j+2]-2.0*UPD[0]);
                    df1=df+z[j][k]*R;
                    df2=df-z[j][k]*R;
                } else {
                    cerr<<"df_flag = "<<df_flag<<" is not implemented.";
                }
                f1+=z[j][k]*df1;
                f2-=z[j][k]*df2;
            }
            *outfile[2*k+1]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f1;
            *outfile[2*k+2]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f2;
            if ((i+1)%6==0) {
                *outfile[2*k+1]<<endl;
                *outfile[2*k+2]<<endl;
            }
        }
    }
    
    for (int k=0;k<Ni*2+1;k++) {
        infile[k]->close();
        outfile[k]->close();
    }
}

void ePump::OptimizePDFs_LHAPDF(void)
{
    ifstream infileInfo;
    ofstream outfileInfo;
    stringstream line;
    vector<shared_ptr<ifstream>> infile;
    vector<shared_ptr<ofstream>> outfile;
    char infilename[80];
    char outfilename[80];
    char outfileFit[80];
    char num[5];
    double dumD;
    int dumI,nx,nt,nfmx,nfval,nx0,nt0,nfmx0,nfval0;
    int npts;
    string inputstring;
    size_t b1, b2;
    string flavorstring;
    string::size_type nfstring=0;
    int nflavors=1;
    d_mat UPD(2*Ni+1);
    
    cout<<"---------------------------------------"<<endl;
    cout<<"     Optimizing PDFs, LHAPDF format "<<endl;
    cout<<"---------------------------------------"<<endl;
    ePout<<"---------------------------------------"<<endl;
    ePout<<"     Optimizing PDFs, LHAPDF format "<<endl;
    ePout<<"---------------------------------------"<<endl;
    
//  Copy input .info file to output .info file
    
    strcpy(infilename,PDFinfile);
    strcat(infilename,".info");
    infileInfo.open(infilename);
    if ( !infileInfo.is_open() ) {
        cerr<<"Error opening input file: "<<infilename<<endl;
        exit(1);
    }
    strcpy(outfilename,PDFoutfile);
    strcat(outfilename,".info");
    outfileInfo.open(outfilename);
    if ( !outfileInfo.is_open() ) {
        cerr<<"Error opening output file: "<<outfilename<<endl;
        exit(1);
    }
    
    getline(infileInfo,inputstring);
    outfileInfo<<inputstring.substr(0,inputstring.length()-1)<<" (Eigenvectors optimized by ePump)'"<<endl;
    getline(infileInfo,inputstring);
    while ( !infileInfo.eof() ) {
        if (inputstring.substr(0,11)=="NumMembers:") {
            outfileInfo<<"NumMembers: "<<2*Ni+1<<endl;
        } else {
            outfileInfo<<inputstring<<endl;
        }
        if (inputstring.substr(0,8)=="Flavors:") {
            b1=inputstring.find("[")+1;
            b2=inputstring.find("]");
            flavorstring=inputstring.substr(b1,b2-b1);
            while ((nfstring=flavorstring.find(',',nfstring)) !=string::npos )
            {
                nflavors++;
                flavorstring.erase(nfstring,1);
            }
            cout<<"nflavors = "<<nflavors<<endl;
            cout<<flavorstring<<endl;
            ePout<<"nflavors = "<<nflavors<<endl;
            ePout<<flavorstring<<endl;
        }
        getline(infileInfo,inputstring);
    }
    infileInfo.close();
    outfileInfo.close();
    
    for (int k=0;k<Ni*2+1;k++) {
        UPD[k].resize(nflavors);
// Open the input files
        strcpy(infilename,PDFinfile);
        strcat(infilename,"_");
        sprintf(num,"%04d",k);
        strcat(infilename,num);
        strcat(infilename,".dat");
        infile.push_back( make_shared<ifstream>(infilename));
        if ( !infile[k]->is_open() ) {
            cerr<<"Error opening input file: "<<infilename<<endl;
            exit(1);
        }
        cout<<infilename<<endl;
        ePout<<infilename<<endl;
// Open the output files
        strcpy(outfilename,PDFoutfile);
        strcat(outfilename,"_");
        strcat(outfilename,num);
        strcpy(outfileFit,outfilename);
        strcat(outfilename,".dat");
        outfile.push_back( make_shared<ofstream>(outfilename));
        if ( !outfile[k]->is_open() ) {
            cerr<<"Error opening output file: "<<outfilename<<endl;
            exit(1);
        }
        cout<<outfilename<<endl;
        ePout<<outfilename<<endl;
        
// Read in headers from .dat file

        do {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        } while (inputstring!=flavorstring);
    }
    
// Read out rotated table files to new .dat files:
    
    double df1,df2,f1,f2;
    getline(*infile[0],inputstring);
    while (inputstring.substr(0,3)!="---") {
        line.str(inputstring);
        for (int i=0;i<nflavors;i++) {
            line>>UPD[0][i];
            for (int k=1;k<Ni*2+1;k++) {
                *infile[k]>>UPD[k][i];
             }
        }
        line.clear();
        for (int k=1;k<Ni*2+1;k++) {
            getline(*infile[k],inputstring);
        }
// Central .dat file
        for (int i=0;i<nflavors;i++) {
            *outfile[0]<<" "<<scientific<<setprecision(8)<<UPD[0][i];
        }
        *outfile[0]<<endl;
// Rotated Eigenset .dat files
        for (int i=0;i<nflavors;i++) {
            for (int k=0;k<Ni;k++) {
                f1=f2=UPD[0][i];
                for (int j=0;j<Ni;j++) {
                  if (df_flag==0||df_flag==2) {
                      df1=UPD[2*j+1][i]-UPD[0][i];
                      df2=UPD[0][i]-UPD[2*j+2][i];
                  } else if (df_flag==1) {
                      df1=0.5*(UPD[2*j+1][i]-UPD[2*j+2][i]);
                      df2=df1;
                  } else if (df_flag==3) {
                      if (z[j][k]>0.0) {
                          df1=UPD[2*j+1][i]-UPD[0][i];
                          df2=UPD[0][i]-UPD[2*j+2][i];
                      } else {
                          df2=UPD[2*j+1][i]-UPD[0][i];
                          df1=UPD[0][i]-UPD[2*j+2][i];
                      }
                  } else if (df_flag==4) {
                      double df=0.5*(UPD[2*j+1][i]-UPD[2*j+2][i]);
                      double R=0.5*(UPD[2*j+1][i]+UPD[2*j+2][i]-2.0*UPD[0][i]);
                      df1=df+z[j][k]*R;
                      df2=df-z[j][k]*R;
                  } else {
                      cerr<<"df_flag = "<<df_flag<<" is not implemented.";
                  }
                  f1+=z[j][k]*df1;
                  f2-=z[j][k]*df2;
                }
                *outfile[2*k+1]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f1;
                *outfile[2*k+2]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f2;
            }
        }
        for (int k=0;k<Ni;k++) {
            *outfile[2*k+1]<<endl;
            *outfile[2*k+2]<<endl;
        }
        getline(*infile[0],inputstring);
    }
    for (int k=0;k<Ni*2+1;k++) {
        *outfile[k]<<"---"<<endl;
        infile[k]->close();
        outfile[k]->close();
    }
}

void ePump::UpdatePDFs_CTEQ(void)
{
    stringstream line;
    vector<shared_ptr<ifstream>> infile;
    vector<shared_ptr<ofstream>> outfile;
    char infilename[80];
    char outfilename[80];
    char outfileFit[80];
    char num[3];
    double dumD;
    int dumI,nx,nt,nfmx,nfval,nx0,nt0,nfmx0,nfval0;
    int npts;
    string inputstring;
    d_vec UPD(2*Ni+1);
    
    cout<<"---------------------------------------"<<endl;
    cout<<"     Updating PDFs, CTEQ format "<<endl;
    cout<<"---------------------------------------"<<endl;
    ePout<<"---------------------------------------"<<endl;
    ePout<<"     Updating PDFs, CTEQ format "<<endl;
    ePout<<"---------------------------------------"<<endl;

    
    for (int k=0;k<Ni*2+1;k++) {
// Open the input files
        strcpy(infilename,PDFinfile);
        strcat(infilename,".");
        sprintf(num,"%02d",k);
        strcat(infilename,num);
        strcat(infilename,".pds");
        infile.push_back( make_shared<ifstream>(infilename));
        if ( !infile[k]->is_open() ) {
            cerr<<"Error opening input file: "<<infilename<<endl;
            exit(1);
        }
        cout<<infilename<<endl;
        ePout<<infilename<<endl;
// Open the output files
        strcpy(outfilename,PDFoutfile);
        strcat(outfilename,".");
        strcat(outfilename,num);
        strcpy(outfileFit,outfilename);
        strcat(outfilename,".pds");
        outfile.push_back( make_shared<ofstream>(outfilename));
        if ( !outfile[k]->is_open() ) {
            cerr<<"Error opening output file: "<<outfilename<<endl;
            exit(1);
        }
        cout<<outfilename<<endl;
        ePout<<outfilename<<endl;
        
// Read in headers from .pds file
        getline(*infile[k],inputstring);
        *outfile[k]<<"PDF Table for Fit #: "<<outfileFit<<endl;
        for (int i=0;i<4;i++) {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        }
        line.str(inputstring);

        line>>dumI>>dumD>>dumI>>dumI>>dumI>>nfmx>>nfval;
        line.clear();
        for (int i=0;i<2;i++) {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        }
        line.str(inputstring);
        line>>nx>>nt;
        line.clear();

        if (k==0) {
            nx0=nx;
            nt0=nt;
            nfmx0=nfmx;
            nfval0=nfval;
            npts=(nx+1)*(nt+1)*(nfval+nfmx+1);
// Assume no photon.  Must modify this if photon included.
        } else {
            if (nx0!=nx||nt0!=nt||nfmx0!=nfmx||nfval0!=nfval) {
                cerr<<"Error: nx,nt,nfmx,nval values do not match in file "<<infilename<<endl;
                exit(1);
            }
        }

        do {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        } while (inputstring!="Parton Distribution Table:");
    }


// Read out rotated table files to new .pds files:
    
    double df1,df2,f1,f2,f0;
    double tp,tm,tbar;
    d_vec Ydyn(Ni);
    double dPDF, ddPDF;
    for (int i=0;i<npts;i++) {
        for (int k=0;k<Ni*2+1;k++) {
            *infile[k]>>UPD[k];
            if ((i+1)%6==0) {
                getline(*infile[k],inputstring);
            }
        }
// Central .pds file
        for (int j=0;j<Ni;j++) {
            
            if (T_flag==1) {
                tp=TdynRMS[j];
                tm=TdynRMS[j];
                tbar=TdynRMS[j];
            } else {
                tp=TdynP[j];
                tm=TdynM[j];
                tbar=TdynBar[j];
            }

            TaylorCoeffs(dPDF,ddPDF,UPD[0],UPD[2*j+1],UPD[2*j+2],tp,tm,tbar);
            Ydyn[j]=dPDF;
        }
        f0=UPD[0]+dot(Ydyn,z0,Ni);
        if (DiagonalQuad) {
            f0+=(dot(Ydyn,z0sqr,Ni)/2.0);
        }

        *outfile[0]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f0;
        if ((i+1)%6==0) {
            *outfile[0]<<endl;
        }
// Rotated Eigenset .pds files
        for (int k=0;k<Ni;k++) {
            f1=f2=0.0;
            for (int j=0;j<Ni;j++) {
                if (df_flag==0||df_flag==2) {
                    df1=UPD[2*j+1]-UPD[0];
                    df2=UPD[0]-UPD[2*j+2];
                } else if (df_flag==1) {
                    df1=0.5*(UPD[2*j+1]-UPD[2*j+2]);
                    df2=df1;
                } else if (df_flag==3) {
                    if (z[j][k]>0.0) {
                        df1=UPD[2*j+1]-UPD[0];
                        df2=UPD[0]-UPD[2*j+2];
                    } else {
                        df2=UPD[2*j+1]-UPD[0];
                        df1=UPD[0]-UPD[2*j+2]; 
                    }
                } else if (df_flag==4) {
                    double df=0.5*(UPD[2*j+1]-UPD[2*j+2]);
                    double R=0.5*(UPD[2*j+1]+UPD[2*j+2]-2.0*UPD[0]);
                    df1=df+z[j][k]*R;
                    df2=df-z[j][k]*R;
                } else {
                    cerr<<"df_flag = "<<df_flag<<" is not implemented.";
                }
                f1+=z[j][k]*df1;
                f2-=z[j][k]*df2;
            }
            f1/=sqrt(1.0+d[k]);
            f2/=sqrt(1.0+d[k]);
            f1+=f0;
            f2+=f0;
            *outfile[2*k+1]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f1;
            *outfile[2*k+2]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f2;
            if ((i+1)%6==0) {
                *outfile[2*k+1]<<endl;
                *outfile[2*k+2]<<endl;
            }
        }
    }
    
    for (int k=0;k<Ni*2+1;k++) {
        infile[k]->close();
        outfile[k]->close();
    }
}

void ePump::UpdatePDFs_LHAPDF(void)
{
    ifstream infileInfo;
    ofstream outfileInfo;
    stringstream line;
    vector<shared_ptr<ifstream>> infile;
    vector<shared_ptr<ofstream>> outfile;
    char infilename[80];
    char outfilename[80];
    char outfileFit[80];
    char num[5];
    double dumD;
    int dumI,nx,nt,nfmx,nfval,nx0,nt0,nfmx0,nfval0;
    int npts;
    string inputstring;
    size_t b1, b2;
    string flavorstring;
    string::size_type nfstring=0;
    int nflavors=1;
    d_mat UPD(2*Ni+1);
    d_vec f0;
    
    cout<<"---------------------------------------"<<endl;
    cout<<"     Updating PDFs, LHAPDF format "<<endl;
    cout<<"---------------------------------------"<<endl;
    ePout<<"---------------------------------------"<<endl;
    ePout<<"     Updating PDFs, LHAPDF format "<<endl;
    ePout<<"---------------------------------------"<<endl;

//  Copy input .info file to output .info file
    
    strcpy(infilename,PDFinfile);
    strcat(infilename,".info");
    infileInfo.open(infilename);
    if ( !infileInfo.is_open() ) {
        cerr<<"Error opening input file: "<<infilename<<endl;
        exit(1);
    }
    strcpy(outfilename,PDFoutfile);
    strcat(outfilename,".info");
    outfileInfo.open(outfilename);
    if ( !outfileInfo.is_open() ) {
        cerr<<"Error opening output file: "<<outfilename<<endl;
        exit(1);
    }
    
    getline(infileInfo,inputstring);
    outfileInfo<<inputstring.substr(0,inputstring.length()-1)<<" (Eigenvectors updated by ePump)'"<<endl;
    getline(infileInfo,inputstring);
    while ( !infileInfo.eof() ) {
        if (inputstring.substr(0,11)=="NumMembers:") {
            outfileInfo<<"NumMembers: "<<2*Ni+1<<endl;
        } else {
            outfileInfo<<inputstring<<endl;
        }
        if (inputstring.substr(0,8)=="Flavors:") {
            b1=inputstring.find("[")+1;
            b2=inputstring.find("]");
            flavorstring=inputstring.substr(b1,b2-b1);
            while ((nfstring=flavorstring.find(',',nfstring)) !=string::npos )
            {
                nflavors++;
                flavorstring.erase(nfstring,1);
            }
            cout<<"nflavors = "<<nflavors<<endl;
            cout<<flavorstring<<endl;
            ePout<<"nflavors = "<<nflavors<<endl;
            ePout<<flavorstring<<endl;
        }
        getline(infileInfo,inputstring);
    }
    infileInfo.close();
    outfileInfo.close();
    
    f0.resize(nflavors);
    
    for (int k=0;k<Ni*2+1;k++) {
        UPD[k].resize(nflavors);
// Open the input files
        strcpy(infilename,PDFinfile);
        strcat(infilename,"_");
        sprintf(num,"%04d",k);
        strcat(infilename,num);
        strcat(infilename,".dat");
        infile.push_back( make_shared<ifstream>(infilename));
        if ( !infile[k]->is_open() ) {
            cerr<<"Error opening input file: "<<infilename<<endl;
            exit(1);
        }
        cout<<infilename<<endl;
        ePout<<infilename<<endl;
// Open the output files
        strcpy(outfilename,PDFoutfile);
        strcat(outfilename,"_");
        strcat(outfilename,num);
        strcpy(outfileFit,outfilename);
        strcat(outfilename,".dat");
        outfile.push_back( make_shared<ofstream>(outfilename));
        if ( !outfile[k]->is_open() ) {
            cerr<<"Error opening output file: "<<outfilename<<endl;
            exit(1);
        }
        cout<<outfilename<<endl;
        ePout<<outfilename<<endl;
        
// Read in headers from .dat file

        do {
            getline(*infile[k],inputstring);
            *outfile[k]<<inputstring<<endl;
        } while (inputstring!=flavorstring);
    }
    
// Read out rotated table files to new .dat files:
    
    double df1,df2,f1,f2;
    double tp,tm,tbar;
    d_vec Ydyn(Ni);
    double dPDF, ddPDF;
    getline(*infile[0],inputstring);
    while (inputstring.substr(0,3)!="---") {
        line.str(inputstring);
        for (int i=0;i<nflavors;i++) {
            line>>UPD[0][i];
            for (int k=1;k<Ni*2+1;k++) {
                *infile[k]>>UPD[k][i];
             }
        }
        line.clear();
        for (int k=1;k<Ni*2+1;k++) {
            getline(*infile[k],inputstring);
        }
// Central .dat file
        for (int i=0;i<nflavors;i++) {
            for (int j=0;j<Ni;j++) {
                if (T_flag==1) {
                    tp=TdynRMS[j];
                    tm=TdynRMS[j];
                    tbar=TdynRMS[j];
                } else {
                    tp=TdynP[j];
                    tm=TdynM[j];
                    tbar=TdynBar[j];
                }

                TaylorCoeffs(dPDF,ddPDF,UPD[0][i],UPD[2*j+1][i],UPD[2*j+2][i],tp,tm,tbar);
                Ydyn[j]=dPDF;
            }
            f0[i]=UPD[0][i]+dot(Ydyn,z0,Ni);
            if (DiagonalQuad) {
                f0[i]+=(dot(Ydyn,z0sqr,Ni)/2.0);
            }
            *outfile[0]<<" "<<scientific<<setprecision(8)<<f0[i];
        }
        *outfile[0]<<endl;
// Rotated Eigenset .dat files
        for (int i=0;i<nflavors;i++) {
            for (int k=0;k<Ni;k++) {
                f1=f2=0.0;
                for (int j=0;j<Ni;j++) {
                  if (df_flag==0||df_flag==2) {
                      df1=UPD[2*j+1][i]-UPD[0][i];
                      df2=UPD[0][i]-UPD[2*j+2][i];
                  } else if (df_flag==1) {
                      df1=0.5*(UPD[2*j+1][i]-UPD[2*j+2][i]);
                      df2=df1;
                  } else if (df_flag==3) {
                      if (z[j][k]>0.0) {
                          df1=UPD[2*j+1][i]-UPD[0][i];
                          df2=UPD[0][i]-UPD[2*j+2][i];
                      } else {
                          df2=UPD[2*j+1][i]-UPD[0][i];
                          df1=UPD[0][i]-UPD[2*j+2][i];
                      }
                  } else if (df_flag==4) {
                      double df=0.5*(UPD[2*j+1][i]-UPD[2*j+2][i]);
                      double R=0.5*(UPD[2*j+1][i]+UPD[2*j+2][i]-2.0*UPD[0][i]);
                      df1=df+z[j][k]*R;
                      df2=df-z[j][k]*R;
                  } else {
                      cerr<<"df_flag = "<<df_flag<<" is not implemented.";
                  }
                  f1+=z[j][k]*df1;
                  f2-=z[j][k]*df2;
                }
                f1/=sqrt(1.0+d[k]);
                f2/=sqrt(1.0+d[k]);
                f1+=f0[i];
                f2+=f0[i];
                *outfile[2*k+1]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f1;
                *outfile[2*k+2]<<setw(13)<<scientific<<uppercase<<setprecision(5)<<f2;
            }
        }
        for (int k=0;k<Ni;k++) {
            *outfile[2*k+1]<<endl;
            *outfile[2*k+2]<<endl;
        }
        getline(*infile[0],inputstring);
    }
    for (int k=0;k<Ni*2+1;k++) {
        *outfile[k]<<"---"<<endl;
        infile[k]->close();
        outfile[k]->close();
    }
}

void ePump::UpdateObservables(void) {

    double Dsym,Dup,Ddown,DnewSym,DnewUp,DnewDown,D0new;
    double DY1,DY2,Cos12,DY1new,DY2new,Cos12new;


    cout<<endl;
    cout<<"------------------------------------------------------------"<<endl;
    cout<<"            Original and Updated Delta Chi-Squares"<<endl;
    cout<<"------------------------------------------------------------"<<endl;
    cout<<endl;
    ePout<<endl;
    ePout<<"------------------------------------------------------------"<<endl;
    ePout<<"            Original and Updated Delta Chi-Squares"<<endl;
    ePout<<"------------------------------------------------------------"<<endl;
    ePout<<endl;
    
    cout<<" Delta Chi2 for original global analysis"<<endl;
    cout<<" Old "<<fixed<<setprecision(4)<<0.0<<endl;
    cout<<" New "<<Chi20new<<endl;
    cout<<endl;
    ePout<<" Delta Chi2 for original global analysis"<<endl;
    ePout<<" Old "<<fixed<<setprecision(4)<<0.0<<endl;
    ePout<<" New "<<Chi20new<<endl;
    ePout<<endl;
    Chi2=0.0;
    Chi2new=0.0;
    int NaTot=0;
    for (int k=0;k<Ndata;k++) {
        if (dataIncluded[k]) {
            NaTot+=Na[k];
            cout<<" Delta Chi2 for Data Set: "<<dataset[k]<<" with "<<Na[k]<<" data points (excluding weight factor w[k] = "<<weight[k]<<")"<<endl;
            cout<<" Old "<<dChi2[k]<<endl;
            cout<<" New "<<dChi2new[k]<<endl;
            ePout<<" Delta Chi2 for Data Set: "<<dataset[k]<<" with "<<Na[k]<<" data points (excluding weight factor w[k] = "<<weight[k]<<")"<<endl;
            ePout<<" Old "<<dChi2[k]<<endl;
            ePout<<" New "<<dChi2new[k]<<endl;
            Chi2+=weight[k]*dChi2[k];
            Chi2new+=weight[k]*dChi2new[k];
        }
    }
    cout<<endl;
    ePout<<endl;
    cout<<" Total Chi2 for new data sets with "<<NaTot<<" data points (including weight factors)"<<endl;
    cout<<" Old "<<Chi2<<endl;
    cout<<" New "<<Chi2new<<endl;
    ePout<<" Total Chi2 for new data sets with "<<NaTot<<" data points (including weight factors)"<<endl;
    ePout<<" Old "<<Chi2<<endl;
    ePout<<" New "<<Chi2new<<endl;
    cout<<endl;
    ePout<<endl;
    
    Chi2new+=Chi20new;
//    double Chi2newAlt=Chi2+Chi2diff;
    cout<<" Total Chi2 for original global analysis plus new data sets (including weight factors)"<<endl;
    cout<<" Old "<<Chi2<<endl;
    cout<<" New "<<Chi2new<<endl;
//    cout<<" New "<<Chi2newAlt<<endl;
    ePout<<" Total Chi2 for original global analysis plus new data sets (including weight factors)"<<endl;
    ePout<<" Old "<<Chi2<<endl;
    ePout<<" New "<<Chi2new<<endl;
//    ePout<<" New "<<Chi2newAlt<<endl;
    cout<<endl;
    ePout<<endl;
    
    
    cout<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<" Original and Updated Shifted Data, Residuals, and Reduced Chi-Squares"<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<endl;
    ePout<<endl;
    ePout<<"-----------------------------------------------------------------------"<<endl;
    ePout<<" Original and Updated Shifted Data, Residuals, and Reduced Chi-Squares"<<endl;
    ePout<<"-----------------------------------------------------------------------"<<endl;
    ePout<<endl;
    
    
    for (int k=0;k<Ndata;k++) {
        if (dataIncluded[k]&&(error_type[k]==1)) {
            double dchi2;
            cout<<" For each data point in Data Set: "<<dataset[k]<<endl;
            cout<<" Old:"<<endl;
            ePout<<" For each data point in Data Set: "<<dataset[k]<<endl;
            ePout<<" Old:"<<endl;
            dchi2=Chi2Residuals(X0[k],k);
            cout<<" New:"<<endl;
            ePout<<" New:"<<endl;
            dchi2=Chi2Residuals(X0new[k],k);
        }
    }
    
    cout<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<" Original and Updated Optimal Error Shifts (Nuisance Parameters)"<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<endl;
    ePout<<endl;
    ePout<<"-----------------------------------------------------------------------"<<endl;
    ePout<<" Original and Updated Optimal Error Shifts (Nuisance Parameters)"<<endl;
    ePout<<"-----------------------------------------------------------------------"<<endl;
    ePout<<endl;
    
    for (int k=0;k<Ndata;k++) {
        if (dataIncluded[k]&&(error_type[k]==1)) {
            double dchi2;
            cout<<" For Data Set: "<<dataset[k]<<endl;
            cout<<" Old:"<<endl;
            ePout<<" For Data Set: "<<dataset[k]<<endl;
            ePout<<" Old:"<<endl;
            dchi2=NuisanceParameters(X0[k],k);
            cout<<" New:"<<endl;
            ePout<<" New:"<<endl;
            dchi2=NuisanceParameters(X0new[k],k);
        }
    }
    
    cout<<endl;
    cout<<"---------------------------------------------------------------"<<endl;
    cout<<"      Original and Updated Predictions and Uncertainties"<<endl;
    cout<<"------------------------------------------------------------"<<endl;
    cout<<endl;
    ePout<<endl;
    ePout<<"---------------------------------------------------------------"<<endl;
    ePout<<"      Original and Updated Predictions and Uncertainties"<<endl;
    ePout<<"---------------------------------------------------------------"<<endl;
    ePout<<endl;
    
    for (int k=0;k<Ndata;k++) {
        cout<<" For Observables X[a] in Data Set: "<<dataset[k]<<endl<<endl;
        cout<<"      a      X[a]      DXsym[a]                 DXasym[a]"<<endl;
        ePout<<" For Observables X[a] in Data Set: "<<dataset[k]<<endl<<endl;
        ePout<<"      a      X[a]      DXsym[a]                 DXasym[a]"<<endl;
        for (int a=0;a<Na[k];a++) {
            UpdateErrors(Xmat[k][a],Dsym,Dup,Ddown,DnewSym,DnewUp,DnewDown,D0new);
            cout<<"Old"<<setw(4)<<a+1<<"  "<<setw(11)<<setprecision(4)<<scientific<<X0[k][a]<<"  "<<Dsym<<"       -"<<Ddown<<"  +"<<Dup<<endl;
            cout<<"New"<<setw(4)<<a+1<<"  "<<setw(11)<<setprecision(4)<<scientific<<D0new<<"  "<<DnewSym<<"       -"<<DnewDown<<"  +"<<DnewUp<<endl<<endl;
            ePout<<"Old"<<setw(4)<<a+1<<"  "<<setw(11)<<setprecision(4)<<scientific<<X0[k][a]<<"  "<<Dsym<<"       -"<<Ddown<<"  +"<<Dup<<endl;
            ePout<<"New"<<setw(4)<<a+1<<"  "<<setw(11)<<setprecision(4)<<scientific<<D0new<<"  "<<DnewSym<<"       -"<<DnewDown<<"  +"<<DnewUp<<endl<<endl;
        }
        cout<<endl;
        ePout<<endl;
    }
    
    if (reportCCs) {
        cout<<endl;
        cout<<"---------------------------------------------------------------"<<endl;
        cout<<"          Original and Updated Correlation Cosines"<<endl;
        cout<<"------------------------------------------------------------"<<endl;
        cout<<endl;
        ePout<<endl;
        ePout<<"---------------------------------------------------------------"<<endl;
        ePout<<"          Original and Updated Correlation Cosines"<<endl;
        ePout<<"---------------------------------------------------------------"<<endl;
        ePout<<endl;

        for (int k=0;k<Ndata;k++) {
            for (int j=k;j<Ndata;j++) {
                cout<<" Original (Old) matrix of Correlation Cosines between X[a] ("<<dataset[k]<<") and Y[b] ("<<dataset[j]<<")"<<endl<<endl;
                ePout<<" Original (Old) matrix of Correlation Cosines between X[a] ("<<dataset[k]<<") and Y[b] ("<<dataset[j]<<")"<<endl<<endl;
                cout<<"     ";
                ePout<<"     ";
                for (int a=0;a<Na[j];a++) {
                    cout<<setw(7)<<"Y["<<setfill('0')<<setw(2)<<a+1<<"]"<<setfill(' ');
                    ePout<<setw(7)<<"Y["<<setfill('0')<<setw(2)<<a+1<<"]"<<setfill(' ');
                }
                cout<<endl;
                ePout<<endl;
                for (int a=0;a<Na[k];a++) {
                    cout<<" X["<<setfill('0')<<setw(2)<<a+1<<"]  "<<setfill(' ');
                    ePout<<" X["<<setfill('0')<<setw(2)<<a+1<<"]  "<<setfill(' ');
                    for (int b=0;b<Na[j];b++) {
                        cout<<setw(8)<<fixed<<setprecision(5);
                        ePout<<setw(8)<<fixed<<setprecision(5);
                        UpdateCorrCosine(Xmat[j][b], Xmat[k][a], DY1, DY2,Cos12, DY1new, DY2new,Cos12new);
                        cout<<Cos12<<"  ";
                        ePout<<Cos12<<"  ";
                    }
                    cout<<endl;
                    ePout<<endl;
                }
                cout<<endl;
                ePout<<endl;

                cout<<" Updated (New) matrix of Correlation Cosines between X[a] ("<<dataset[k]<<") and Y[b] ("<<dataset[j]<<")"<<endl<<endl;
                ePout<<" Updated (New) matrix of Correlation Cosines between X[a] ("<<dataset[k]<<") and Y[b] ("<<dataset[j]<<")"<<endl<<endl;
                cout<<"     ";
                ePout<<"     ";
                for (int a=0;a<Na[j];a++) {
                    cout<<setw(7)<<"Y["<<setfill('0')<<setw(2)<<a+1<<"]"<<setfill(' ');
                    ePout<<setw(7)<<"Y["<<setfill('0')<<setw(2)<<a+1<<"]"<<setfill(' ');
                }
                cout<<endl;
                ePout<<endl;
                for (int a=0;a<Na[k];a++) {
                    cout<<" X["<<setfill('0')<<setw(2)<<a+1<<"]  "<<setfill(' ');
                    ePout<<" X["<<setfill('0')<<setw(2)<<a+1<<"]  "<<setfill(' ');
                    for (int b=0;b<Na[j];b++) {
                        cout<<setw(8)<<fixed<<setprecision(5);
                        ePout<<setw(8)<<fixed<<setprecision(5);
                        UpdateCorrCosine(Xmat[j][b], Xmat[k][a], DY1, DY2,Cos12, DY1new, DY2new,Cos12new);
                        cout<<Cos12new<<"  ";
                        ePout<<Cos12new<<"  ";
                    }
                    cout<<endl;
                    ePout<<endl;
     
                }
            }
        }
    }
    cout<<endl;
    ePout<<endl;
    if (pdftype!=NONE) {
        UpdatePDFs();
    }
}


void ePump::OptimizeObservables(void) {
    OptimizedEVContributions();
    OptimizePDFs();
}





double MaxComp(const double &a, const double &b){
    if ((a<0.0)&&(b>0.0)) {return 0.0;}
    else if (a>-b) {return a;}
    else {return b;}
}

void TaylorCoeffs(double &dx, double &ddx, const double &X0, const double &Xp, const double &Xm, const double &Tp, const double &Tm, const double &Tbar) {
//    dx=(Xp-Xm)/sqrt(Tp*Tp+Tm*Tm)/sqrt(2);
    dx=((Tm/Tbar)*(Xp-X0)/Tp+(Tp/Tbar)*(X0-Xm)/Tm)/2.0;
    ddx=((Xp-X0)/Tp-(X0-Xm)/Tm)/Tbar;
}

