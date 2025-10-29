#include "ePump.h"
#include "Cholesky.h"
#include "random.h"

// This routine needs to be updated in order to use
// the current correct formats for ".theory" and ".data" files.

void ePump::MakePseudoData(const char* Xfile, const char* Datafile) {

    ReadInXi(Xfile);
    SetQ();
    
    ofstream outfile;
    outfile.open(Datafile);
    if ( !outfile.is_open() ) {
        cout<<"Error opening data file."<<endl;
    } else {
        cout<<"Writing PseudoData arrays."<<endl;
    }
    
    for (int a=0;a<Na;a++) {
        for (int b=0;b<Na;b++) {
            Cm[a][b]=0.0;
        }
    }

    outfile<<"Data for Experiments"<<endl;
    outfile<<endl;
    outfile<<"a           data_value"<<endl;
    
    Random r;
    r.seed(-11);
    for (int a=0;a<Na;a++) {
        XE[a]=r.normal(X0[a],sqrt(Q[a][a]/1.645));
        double fac=4.0;
        Cm[a][a]=fac*fac*(1.645*1.645)/Q[a][a];
        outfile<<"XE["<<a<<"]      "<<XE[a]<<endl;
    }
    
    outfile<<endl;
    outfile<<"Inverse Error Matrix, diagonal elements:"<<endl;
    for (int a=0;a<Na;a++) {
        outfile<<"Cm["<<a<<"]["<<a<<"]     "<<Cm[a][a]<<endl;
    }

    outfile<<endl;
    outfile<<"Inverse Error Matrix, off-diagonal elements:"<<endl;

    for (int a=0;a<Na-1;a++) {
        outfile<<"a="<<a<<",b="<<a+1<<"-"<<Na-1<<endl;
        for (int b=a+1;b<Na;b++) {
            outfile<<fixed<<setprecision(1)<<Cm[a][b]<<" ";
        }
        outfile<<endl;
    }

    outfile<<endl;
    outfile.close();

}

int main()
{
    int Ni=28;    // Number of eigenvector pairs
    int Na=16;    // Number of data points
    ePump EU(Ni,Na);
    EU.MakePseudoData("bin/cosinethetastar_dXi","bin/cosinethetastar_PseudoData");
    return 0;
}

