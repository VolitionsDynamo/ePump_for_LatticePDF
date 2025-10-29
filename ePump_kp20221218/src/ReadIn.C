#include "ePump.h"
#include "Cholesky.h"
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>


using namespace std;

void ePump::ReadInData1(char* Datafile, int k) {
    
    ifstream infile;
    string rawstring, inputstring;
    string dname;
    
    infile.open(Datafile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening data file: "<<Datafile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"Initializing Data arrays for Data Set "<<k+1<<": "<<Datafile<<" , error_type "<<error_type[k]<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"Initializing Data arrays for Data Set "<<k+1<<": "<<Datafile<<" , error_type "<<error_type[k]<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
    }
    cout<<endl;
    ePout<<endl;
    
    XE[k].resize(Na[k]);
    Cm[k].resize(Na[k]);
    for (int a=0;a<Na[k];a++) {
        Cm[k][a].resize(Na[k]);
    }
    
    do {
        getline(infile, rawstring);
				inputstring=trim(rawstring);
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*');

    double dum;
    double StatErr, UncSys, epercent;
    int Dcol, SEcol, UScol, CEcol;
    infile>>Nlam[k]>>Dcol>>SEcol>>UScol;
    cout<<"             "<<Nlam[k]<<"                "<<Dcol<<"                "<<SEcol<<"                "<<UScol;
    ePout<<"             "<<Nlam[k]<<"                "<<Dcol<<"                "<<SEcol<<"                "<<UScol;
        
    int dSE=SEcol-Dcol;
    int dUS=UScol-SEcol;
    int dCE=0;
    s[k].resize(Na[k]);
    beta[k].resize(Na[k]);
    calAinvBeta[k].resize(Na[k]);
    if (Nlam[k]!=0) {
        infile>>CEcol;
        cout<<"               "<<CEcol;
        ePout<<"               "<<CEcol;
        dCE=CEcol-UScol;
        for (int a=0;a<Na[k];a++) {
            beta[k][a].resize(Nlam[k]);
            calAinvBeta[k][a].resize(Nlam[k]);
        }
    }
    cout<<endl<<endl;
    ePout<<endl<<endl;
        
    cout<<"  Inverse covariance matrix to be constructed from "<<endl;
    cout<<"  uncorrelated errors and correlated systematic errors."<<endl;
    cout<<"  Nlam = "<<Nlam[k]<<"  weight (including global weight factor) = "<<weight[k]<<endl;
    ePout<<"  Inverse covariance matrix to be constructed from "<<endl;
    ePout<<"  uncorrelated errors and correlated systematic errors."<<endl;
    ePout<<"  Nlam = "<<Nlam[k]<<"  weight (including global weight factor) = "<<weight[k]<<endl;
    getline(infile, inputstring);
    getline(infile, rawstring);
    inputstring=trim(rawstring);
			
    cout<<endl;
    ePout<<endl;
    cout<<setw(2)<<right<<"a"<<setw(9)<<"      XE[a] "<<setw(12)<<"     s[a] "<<"      ";
    ePout<<setw(2)<<right<<"a"<<setw(9)<<"      XE[a] "<<setw(12)<<"     s[a] "<<"      ";
    for (int j=0;j<Nlam[k];j++) {
        cout<<setw(10)<<"beta[a]["<<setfill('0')<<setw(2)<<j+1<<"] "<<setfill(' ');
        ePout<<setw(10)<<"beta[a]["<<setfill('0')<<setw(2)<<j+1<<"] "<<setfill(' ');
    }
    cout<<endl;
    ePout<<endl;
    for (int a=0;a<Na[k];a++) {
        for (int i=1;i<Dcol;i++) {
            infile>>dum;
        }
        infile>>XE[k][a];
        for (int i=1;i<dSE;i++) {
            infile>>dum;
        }
        infile>>StatErr;
        for (int i=1;i<dUS;i++) {
            infile>>dum;
        }
        infile>>UncSys;
        for (int i=1;i<dCE;i++) {
            infile>>dum;
        }
        s[k][a]=sqrt(StatErr*StatErr+UncSys*UncSys);
        cout<<setw(2)<<right<<a+1<<"   "<<setw(9)<<XE[k][a]<<" "<<setw(12)<<s[k][a]<<"  ";
        ePout<<setw(2)<<right<<a+1<<"   "<<setw(9)<<XE[k][a]<<" "<<setw(12)<<s[k][a]<<"  ";
        for (int j=0;j<Nlam[k];j++) {
            infile>>epercent;
            if (error_type[k]==1) {
                beta[k][a][j]=epercent*X0[k][a]/100.0;
            } else if (error_type[k]==4) {
                beta[k][a][j]=epercent*XE[k][a]/100.0;
            } else {
                cerr<<"  It should be impossible to get here!!!"<<endl;
                exit(1);
            }
            cout<<setw(14)<<beta[k][a][j];
            ePout<<setw(14)<<beta[k][a][j];
        }
        cout<<endl;
        ePout<<endl;
        getline(infile, rawstring);
        inputstring=trim(rawstring);
				
    }
    cout<<endl;
    ePout<<endl;
    ConstructCm1(k);

    dChi2[k]=DeltaChiSquare(X0[k],k);
//    cout<<"  Delta Chi^2 for Data Set "<<dataset[k]<<" :   "<<dChi2[k]<<endl;
//    ePout<<"  Delta Chi^2 for Data Set "<<dataset[k]<<" :   "<<dChi2[k]<<endl;
    cout<<endl;
    ePout<<endl;
    
    infile.close();

}

void ePump::ReadInData2(char* Datafile, int k) {
    
    ifstream infile;
    string rawstring, inputstring;
    string dname;
    
    infile.open(Datafile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening data file: "<<Datafile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"Initializing Data arrays for Data Set "<<k+1<<": "<<Datafile<<" , error_type "<<error_type[k]<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"Initializing Data arrays for Data Set "<<k+1<<": "<<Datafile<<" , error_type "<<error_type[k]<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
    }
    cout<<endl;
    ePout<<endl;
    
    XE[k].resize(Na[k]);
    Cm[k].resize(Na[k]);
    rhocc[k].resize(Na[k]);
    for (int a=0;a<Na[k];a++) {
        Cm[k][a].resize(Na[k]);
        rhocc[k][a].resize(Na[k]);
    }
    
    do {
        getline(infile, rawstring);
        inputstring=trim(rawstring);
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*');

    double dum;
    double StatErr, UncSys, epercent;
    int Dcol, SEcol, UScol, CEcol;
    infile>>Nlam[k]>>Dcol>>SEcol>>UScol;
    cout<<"             "<<Nlam[k]<<"                "<<Dcol<<"                "<<SEcol<<"                "<<UScol;
    ePout<<"             "<<Nlam[k]<<"                "<<Dcol<<"                "<<SEcol<<"                "<<UScol;
    int dSE=SEcol-Dcol;
    int dUS=UScol-SEcol;
    int dCE=0;
    s_sys[k].resize(Na[k]);
    s_stat[k].resize(Na[k]);
    beta[k].resize(Na[k]);
    if (Nlam[k]!=0) {
        infile>>CEcol;
        cout<<"               "<<CEcol;
        ePout<<"               "<<CEcol;
        dCE=CEcol-UScol;
        for (int a=0;a<Na[k];a++) {
            beta[k][a].resize(Nlam[k]);
        }
    }
    cout<<endl<<endl;
    ePout<<endl<<endl;
        
    cout<<"  Inverse covariance matrix to be constructed from "<<endl;
    cout<<"  uncorrelated errors, correlation coefficients, and correlated systematic errors."<<endl;
    cout<<"  Nlam = "<<Nlam[k]<<"  weight (including global weight factor) = "<<weight[k]<<endl;
    ePout<<"  Inverse covariance matrix to be constructed from "<<endl;
    ePout<<"  uncorrelated errors, correlation coefficients, and correlated systematic errors."<<endl;
    ePout<<"  Nlam = "<<Nlam[k]<<"  weight (including global weight factor) = "<<weight[k]<<endl;
    getline(infile, inputstring);
    getline(infile, rawstring);
    inputstring=trim(rawstring);
		
    cout<<endl;
    ePout<<endl;
    cout<<setw(2)<<right<<"a"<<setw(9)<<"      XE[a] "<<setw(12)<<"s_stat[a] "<<setw(12)<<" s_sys[a] "<<"      ";
    ePout<<setw(2)<<right<<"a"<<setw(9)<<"      XE[a] "<<setw(12)<<"s_stat[a] "<<setw(12)<<" s_sys[a] "<<"      ";
    for (int j=0;j<Nlam[k];j++) {
        cout<<setw(10)<<"beta[a]["<<setfill('0')<<setw(2)<<j+1<<"] "<<setfill(' ');
        ePout<<setw(10)<<"beta[a]["<<setfill('0')<<setw(2)<<j+1<<"] "<<setfill(' ');
    }
    cout<<endl;
    ePout<<endl;
    for (int a=0;a<Na[k];a++) {
        for (int i=1;i<Dcol;i++) {
            infile>>dum;
        }
        infile>>XE[k][a];
        for (int i=1;i<dSE;i++) {
            infile>>dum;
        }
        infile>>StatErr;
        for (int i=1;i<dUS;i++) {
            infile>>dum;
        }
        infile>>UncSys;
        for (int i=1;i<dCE;i++) {
            infile>>dum;
        }
        s_sys[k][a]=UncSys;
        s_stat[k][a]=StatErr;
        cout<<setw(2)<<right<<a+1<<"   "<<setw(9)<<XE[k][a]<<" "<<setw(12)<<s_stat[k][a]<<"  "<<setw(12)<<s_sys[k][a]<<"  ";
        ePout<<setw(2)<<right<<a+1<<"   "<<setw(9)<<XE[k][a]<<" "<<setw(12)<<s_stat[k][a]<<"  "<<setw(12)<<s_sys[k][a]<<"  ";
        for (int j=0;j<Nlam[k];j++) {
//
//  Use this if correlated errors given as percents (multiplicative):
//            infile>>epercent;
//            beta[k][a][j]=epercent*X0[k][a]/100.0;
//
// Use this if correlated errors given as un-normalized and additive:
           infile>>beta[k][a][j];
//
// Use this to change correlated errors from additive to multiplicative:
//            beta[k][a][j]*=(X0[k][a]/XE[k][a]);
//
            cout<<setw(14)<<beta[k][a][j];
            ePout<<setw(14)<<beta[k][a][j];
        }
        cout<<endl;
        ePout<<endl;
        getline(infile, rawstring);
        inputstring=trim(rawstring);
				
    }
    cout<<endl;
    ePout<<endl;
    do {
        getline(infile, rawstring);
        inputstring=trim(rawstring);
				
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*');
    int aa,bb;
    double cc;
		if(inputstring.find("Matrix") != string::npos){
			for(int a=0;a<Na[k];a++){
				for(int b=0;b<Na[k];b++){
					infile>>rhocc[k][a][b];
				}
			}
		}
		else{
    	for (int a=0;a<Na[k];a++) {
				// change to read all the correlation matrix, rather half of them
        for (int b=0;b<Na[k];b++) {
            int aaa, bbb;
            infile>>aa>>bb;
            aaa=aa-1;
            bbb=bb-1;
            infile>>rhocc[k][aaa][bbb];
            getline(infile,inputstring);
            //cout<<setw(10)<<right<<aa<<"         "<<right<<bb<<"    "<<setw(4)<<rhocc[k][aaa][bbb]<<endl;
            //ePout<<setw(10)<<right<<aa<<"         "<<right<<bb<<"    "<<setw(4)<<rhocc[k][aaa][bbb]<<endl;
            if (bbb<aaa) {
                if(rhocc[k][bbb][aaa]!=rhocc[k][aaa][bbb]){
									cout<<bbb<<","<<aaa<<" "<<rhocc[k][bbb][aaa]
										<<" Not equal to "<< aaa<<","<<bbb<<" "<<rhocc[k][aaa][bbb]<<endl;
									//deal with the case that corr[b][a]!=corr[a][b]
									//It become necessary when the data table only contains upper triangle elements.
									rhocc[k][bbb][aaa]=rhocc[k][aaa][bbb];
								}		
            }// deal with lower triangle elements
        }//loop b
    	}//loop a
		}// default ET2 input
    ConstructCm2(k);
    
    dChi2[k]=DeltaChiSquare(X0[k],k);
    cout<<"  Delta Chi^2 for Data Set "<<dataset[k]<<" :   "<<dChi2[k]<<endl;
    ePout<<"  Delta Chi^2 for Data Set "<<dataset[k]<<" :   "<<dChi2[k]<<endl;
    cout<<endl;
    ePout<<endl;
    
    infile.close();
    
}

void ePump::ReadInData3(char* Datafile, int k) {
    
    ifstream infile;
    string rawstring, inputstring;
    string dname;
    
    infile.open(Datafile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening data file: "<<Datafile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"Initializing Data arrays for Data Set "<<k+1<<": "<<Datafile<<" , error_type "<<error_type[k]<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"Initializing Data arrays for Data Set "<<k+1<<": "<<Datafile<<" , error_type "<<error_type[k]<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
    }
    cout<<endl;
    ePout<<endl;

    //Na[k]=N(Observables) is read from .in file
    XE[k].resize(Na[k]);
    Cm[k].resize(Na[k]);
    for (int a=0;a<Na[k];a++) {
        Cm[k][a].resize(Na[k]);
    }
    
    do {
        getline(infile, rawstring);
        inputstring=trim(rawstring);		
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*'); // comments
    for (int a=0;a<Na[k];a++) {
        infile>>XE[k][a];
        getline(infile, rawstring);
        inputstring=trim(rawstring);				
        cout<<" XE["<<a<<"]  "<<XE[k][a]<<endl;
        ePout<<" XE["<<a<<"]  "<<XE[k][a]<<endl;
    }// read data
    do {
        getline(infile, rawstring);
        inputstring=trim(rawstring);
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*'); //comments
		if(inputstring.find("Inverse") != string::npos) {
			//default inverse covariance matrix
    	for (int a=0;a<Na[k];a++) {
        getline(infile, rawstring);
        inputstring=trim(rawstring);				
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
        for (int b=a;b<Na[k];b++) { //b=a..Na[k], the upper triangle elements
            infile>>Cm[k][a][b];
            if (a!=b) {
                Cm[k][b][a]=Cm[k][a][b]; // read ICov
            }
            cout<<Cm[k][a][b]<<"  ";
            ePout<<Cm[k][a][b]<<"  ";
        }
        getline(infile, rawstring); // skip line
        inputstring=trim(rawstring);				
        cout<<endl;
        ePout<<endl;
			}
		}
		else{//read covariance matrix rather than inverse
     for (int a=0;a<Na[k];a++) {
			 for(int b=0;b<Na[k];b++){// We read the whole N*N elements, rather than half
			 	infile>>Cm[k][a][b];
				if(b<a){
					if(Cm[k][a][b]!=Cm[k][b][a]){
						cout<<b+1<<","<<a+1<<" "<<Cm[k][b][a]
            <<" Not equal to "<<a+1<<","<<b+1<<" "<<Cm[k][a][b]<<endl;
                  //deal with the case that cov[b][a]!=cov[a][b]
                  //It become necessary when the data table only contains upper triangle elements.
                  Cm[k][b][a]=Cm[k][a][b];
					}
				}
			 }
		 	}
			//Convert into Cm
			ConstructCm3(k);
		}

    dChi2[k]=DeltaChiSquare(X0[k],k);
    cout<<"  Delta Chi^2 for Data Set "<<dataset[k]<<" :   "<<dChi2[k]<<endl;
    ePout<<"  Delta Chi^2 for Data Set "<<dataset[k]<<" :   "<<dChi2[k]<<endl;
    cout<<endl;
    ePout<<endl;
    
    infile.close();
    
}


void ePump::ReadInTheory(char* Theoryfile, int k) {
    ifstream infile;
    string rawstring, inputstring;
    stringstream line;
    double dum, dXmin;
    int Ncol;
    
    infile.open(Theoryfile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening theory file: "<<Theoryfile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"Initializing Theory arrays for Data Set "<<k+1<<": "<<Theoryfile<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
        ePout<<"Initializing Theory arrays for Data Set "<<k+1<<": "<<Theoryfile<<endl;
        ePout<<"-----------------------------------------------------------------"<<endl;
    }
    cout<<endl;
    ePout<<endl;
    
    X0[k].resize(Na[k]);
    X0new[k].resize(Na[k]);
    Xmat[k].resize(Na[k]);
		Xset[k].resize(2*Ni+1);
    dX[k].resize(Na[k]);
    Q[k].resize(Na[k]);

    for (int a=0;a<Na[k];a++) {
        Xmat[k][a].resize(2*Ni+1);
        dX[k][a].resize(Ni);
        Q[k][a].resize(Na[k]);
    }
		for (int i=0;i<2*Ni+1;i++) Xset[k][i].resize(Na[k]);
    
    do {
        getline(infile, rawstring);
        inputstring=trim(rawstring);			
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*');
    infile>>Ncol;
    getline(infile, inputstring);
    cout<<"      "<<Ncol<<endl;
    ePout<<"      "<<Ncol<<endl;
    
// Set up to output the theory file in with PDF number along columns
// and Observable along rows.
    
    cout<<"    D1         D2         D3";
    ePout<<"    D1         D2         D3";
    for (int i=0;i<2*Ni+1;i++) {
        cout<<"        X["<<setfill('0')<<setw(2)<<i<<"]";
        ePout<<"        X["<<setfill('0')<<setw(2)<<i<<"]";
    }
    cout<<endl;
    ePout<<endl;
    
    int Ndum;
    if (Ncol<4) {
        Ndum=3;
    } else {
        Ndum=Ncol-1;
    }
    d_mat D(Ndum);
    for (int b=0;b<Ndum;b++) {
        D[b].resize(Na[k]);
    }

		getline(infile, rawstring);
		inputstring=trim(rawstring);
		
    for (int a=0;a<Na[k];a++) {
        getline(infile, rawstring);
        inputstring=trim(rawstring);
				
        line.str(inputstring);
        for (int b=0;b<Ncol-1;b++) {
            line>>D[b][a];
        }
        if (Ncol<4) {
            for (int b=Ncol-1;b<Ndum;b++) {
                D[b][a]=0.0;
            }
        }
        line>>X0[k][a];
        Xmat[k][a][0]=X0[k][a];
				Xset[k][0][a]=Xmat[k][a][0];
        line.clear();
    }
    
    for (int i=0;i<Ni;i++) {
        getline(infile, rawstring); inputstring=trim(rawstring);
				
        for (int a=0;a<Na[k];a++) {
            getline(infile,inputstring);  line.str(inputstring);
            for (int b=0;b<Ncol-1;b++) {
                line>>dum;
            }
            line>>dX[k][a][i];
            Xmat[k][a][2*i+1]=dX[k][a][i];
						Xset[k][2*i+1][a]=Xmat[k][a][2*i+1];
            line.clear();
        }
        getline(infile, rawstring);  inputstring=trim(rawstring);				
        for (int a=0;a<Na[k];a++) {
						getline(infile, rawstring);  inputstring=trim(rawstring);
            line.str(inputstring);
            for (int b=0;b<Ncol-1;b++) {
                line>>dum;
            }
            line>>dXmin;
            Xmat[k][a][2*i+2]=dXmin;
						Xset[k][2*i+2][a]=Xmat[k][a][2*i+2];
            line.clear();
            dX[k][a][i]-=dXmin;
            dX[k][a][i]/=2.0;
        }
    }
    
    cout<<setfill(' ');
    ePout<<setfill(' ');
    for (int a=0;a<Na[k];a++) {
        cout<<setw(9)<<scientific<<uppercase<<setprecision(3)<<D[0][a];
        cout<<setw(11)<<scientific<<uppercase<<setprecision(3)<<D[1][a];
        cout<<setw(11)<<scientific<<uppercase<<setprecision(3)<<D[2][a];
        ePout<<setw(9)<<scientific<<uppercase<<setprecision(3)<<D[0][a];
        ePout<<setw(11)<<scientific<<uppercase<<setprecision(3)<<D[1][a];
        ePout<<setw(11)<<scientific<<uppercase<<setprecision(3)<<D[2][a];
        for (int i=0;i<2*Ni+1;i++) {
            cout<<setw(13)<<scientific<<uppercase<<setprecision(5)<<Xmat[k][a][i];
            ePout<<setw(13)<<scientific<<uppercase<<setprecision(5)<<Xmat[k][a][i];
        }
        cout<<endl;
        ePout<<endl;
    }
    
    cout<<endl;
    ePout<<endl;
    
    infile.close();
}



void ePump::SetCinv_Max(int k) {
    
    Cm[k].resize(Na[k]);
    for (int a=0;a<Na[k];a++) {
        Cm[k][a].resize(Na[k]);
        for (int b=0;b<Na[k];b++) {
            if (a==b) {
                Cm[k][a][b]=1.0/Q[k][a][b];
            } else {
                Cm[k][a][b]=0.0;
            }
        }
    }

}


void ePump::ReadInTheoryAndData(void) {
    char mainInfile[80];
    char Theoryfile[80];
    char Datafile[80];
    ifstream infile;
    string rawstring, inputstring;
    string pdfstyle;
    string PDFin;
    string PDFout;
    string dataYN;
    
    strcpy(mainInfile,mainfile);
    strcat(mainInfile,".in");
    infile.open(mainInfile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening input file: "<<mainInfile<<" "<<mainfile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------"<<endl;
        cout<<"Initializing Main input from file: "<<mainInfile<<endl;
        cout<<"-----------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------"<<endl;
        ePout<<"Initializing Main input from file: "<<mainInfile<<endl;
        ePout<<"-----------------------------------------------------------"<<endl;
    }
    
    cout<<endl;
    ePout<<endl;
    if (Pweight!=1.0) {
        cout<<"Contribution of Original data sets suppressed by factor p="<<Pweight<<" in update."<<endl<<endl;
        ePout<<"Contribution of Original data sets suppressed by factor p="<<Pweight<<" in update."<<endl<<endl;
    }
    if (DiagonalQuad) {
        cout<<"Including Diagonal Quadratic terms in X(z)"<<endl<<endl;
        ePout<<"Including Diagonal Quadratic terms in X(z)"<<endl<<endl;
    }
    if ((T_flag==1)&&(df_flag==0)&&(!DiagonalQuad)) {
        cout<<"Using old defaults for DX[i]"<<endl<<endl;
        ePout<<"Using old defaults for DX[i]"<<endl<<endl;
    }

    cout<<"Global weight factor for New data sets w="<<Gweight<<endl;
    ePout<<"Global weight factor for New data sets w="<<Gweight<<endl;
    cout<<endl;
    ePout<<endl;
    
    getline(infile, inputstring);
    cout<<inputstring<<endl;
    ePout<<inputstring<<endl;
    
    string DTYN, QUAD, SCAL;
    bool DynTol;
    
    infile>>Ni>>Ndata>>pdfstyle;
		infile >> QUAD ; // include diagonal quadratic
		if((QUAD=="Y")||(QUAD=="y")){
			DiagonalQuad=true;
		}
		else{
			DiagonalQuad=false;
		}

		infile>> DTYN; // dynamical tolerance
    if ((DTYN=="Y")||(DTYN=="y")) {
        DynTol=true;
        Tsq=100.0;
    } else {
        DynTol=false;
        infile>>Tsq;
    }
    getline(infile, inputstring);
    cout<<"        "<<Ni<<"             "<<Ndata<<"              "<<pdfstyle<<"              "<<DTYN<<"              "<<Tsq<<endl;
    ePout<<"        "<<Ni<<"             "<<Ndata<<"              "<<pdfstyle<<"              "<<DTYN<<"              "<<Tsq<<endl;
    if ((pdfstyle=="C")||(pdfstyle=="c")) {
        pdftype=CTEQ;
    } else if ((pdfstyle=="L")||(pdfstyle=="l")) {
        pdftype=LHAPDF;
    } else if ((pdfstyle=="N")||(pdfstyle=="n")) {
        pdftype=NONE;
    } else {
        cerr<<"PDFtype must be C, L, or N.  You input "<<pdfstyle<<endl;
        exit(1);
    }
//
    Tfix=sqrt(Tsq);
    Na.resize(Ndata);
    error_type.resize(Ndata);
    dataIncluded.resize(Ndata);
    weight.resize(Ndata);
    Xmat.resize(Ndata);
		Xset.resize(Ndata);
    dX.resize(Ndata);
    Xdyn.resize(Ndata);
    XXdyn.resize(Ndata);
    X0.resize(Ndata);
    X0new.resize(Ndata);
    Q.resize(Ndata);
    dChi2.resize(Ndata);
    dChi2new.resize(Ndata);
    XE.resize(Ndata);
    Cm.resize(Ndata);
    rhocc.resize(Ndata);
    dataset.resize(Ndata);
		theoryset.resize(Ndata);

//
// Re-diagonalization Objects:
//    LinvdX.resize(Ndata);
    M.resize(Ni);
    A.resize(Ni);
    L.resize(Ni);
    LinvA.resize(Ni);
    TdynP.resize(Ni);
    TdynM.resize(Ni);
    TdynBar.resize(Ni);
    TdynRMS.resize(Ni);
//
    for (int i=0;i<Ni;i++) {
        M[i].resize(Ni);
        L[i].resize(Ni);
        TdynP[i]=1.0;
        TdynM[i]=1.0;
        TdynBar[i]=1.0;
        TdynRMS[i]=1.0;
    }
    
    s.resize(Ndata);
    s_sys.resize(Ndata);
    s_stat.resize(Ndata);
    beta.resize(Ndata);
    calAinvBeta.resize(Ndata);
    Nlam.resize(Ndata);
    
    getline(infile, inputstring);
    cout<<inputstring<<endl;
    ePout<<inputstring<<endl;

    for (int i=0;i<Ndata;i++) {
			getline(infile, rawstring);
			inputstring=trim(rawstring); //cout<<inputstring<<endl;
			stringstream ss;
			ss<<inputstring; cout<<"stringstream: "<<ss.str()<<endl;
      ss>>dataset[i]>>Na[i]>>dataYN;
        if ((dataYN=="Y")||(dataYN=="y")) {
            dataIncluded[i]=true;
            ss>>error_type[i]>>weight[i];
            cout<<" Data set: "<<dataset[i]<<" Npt: "<<Na[i]<<" Y/N: "<<dataYN<<" error_type: "<<error_type[i]<<" weight: "<<weight[i]<<endl;
            ePout<<" "<<setw(15)<<left<<dataset[i]<<right<<setw(15)<<Na[i]<<"                "<<dataYN<<"           "<<error_type[i]<<"          "<<weight[i]<<endl;
            weight[i]*=Gweight;
						if(!ss.eof()) ss>> theoryset[i];
         } else {
            dataIncluded[i]=false;
						ss>>error_type[i]>>weight[i];
            cout<<"      "<<dataset[i]<<"            "<<Na[i]<<"                "<<dataYN<<endl;
            ePout<<" "<<setw(15)<<left<<dataset[i]<<right<<setw(15)<<Na[i]<<"                "<<dataYN<<endl;
						if(!ss.eof()) ss>> theoryset[i];
        }
				//getline(infile, inputstring);
				//cout<<"After ss: "<<inputstring<<endl;
    }//end loop Ndata
    
    if ((pdftype!=NONE)||DynTol) {
        getline(infile, inputstring);
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
        infile>>PDFin;
        cout<<"    "<<PDFin;
        ePout<<"    "<<PDFin;
        strcpy(PDFinfile,PDFin.c_str());
    }
        
    if (pdftype!=NONE) {
        infile>>PDFout;
        cout<<"        "<<PDFout<<endl;
        ePout<<"        "<<PDFout<<endl;
        strcpy(PDFoutfile,PDFout.c_str());
    }
        
    cout<<endl;
    ePout<<endl;
    infile.close();
    
    for (int i=0;i<Ndata;i++) {
        //strcpy(Theoryfile,dataset[i].c_str());
				if(theoryset[i].empty()) theoryset[i]=dataset[i];
				cout<<"Theoryfile: "<<theoryset[i]<<endl;
        if(theoryset[i].find(".theory")==string::npos)
					theoryset[i]=theoryset[i]+".theory";
        ReadInTheory(const_cast<char*>(theoryset[i].c_str()),i);
        SetQ(i);
//        X00new[i].resize(Na[i]);
//        Qnew[i].resize(Na[i]);
//        LinvdX[i].resize(Na[i]);
//        for (int a=0;a<Na[i];a++) {
//            Qnew[i][a].resize(Na[i]);
//            LinvdX[i][a].resize(Ni);
//        }
        if (dataIncluded[i]) {
            strcpy(Datafile,dataset[i].c_str());
            strcat(Datafile,".data");
            switch (error_type[i]) {
                case 1: ReadInData1(Datafile,i);
                        break;
                case 2: ReadInData2(Datafile,i);
                        break;
                case 3: ReadInData3(Datafile,i);
                        break;
                case 4: ReadInData1(Datafile,i);
                //error_types 1 and 4 are identical, except
                //1: the correlated systematic errors are normalized to the theory
                //4: the correlated systematic errors are normalized to the data
                        break;
                default: { cerr<<"  Error Type must be 1, 2, 3, or 4."<<endl;
                           exit(1);
                }
            }//switch error_type
        }//data Included
			for(int j=0;j<2*Ni+1; j++){
				cout<<"data set: "<<i<<" PDF set: "<<j<<" chi^2: "<<setprecision(3)<<DeltaChiSquare(Xset[i][j],i)<<endl;
				ePout<<"data set: "<<i<<" PDF set: "<<j<<" chi^2: "<<setprecision(3)<<DeltaChiSquare(Xset[i][j],i)<<endl;
			}
    }//Ndata

    if (DynTol) {
        ReadInTolerances();
    }
    for (int i=0;i<Ndata;i++) {
        SetXdyn(i);
    }
    
}

void ePump::ReadInTheory_Optimize(void) {
    char mainInfile[80];
    char Theoryfile[80];
    
    ifstream infile;
    string inputstring;
    string pdfstyle;
    string PDFin;
    string PDFout;
    
    strcpy(mainInfile,mainfile);
		string mainInfilestring=mainInfile;
		if(mainInfilestring.find(".in")==string::npos) strcat(mainInfile,".in");
    infile.open(mainInfile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening input file: "<<mainInfile<<" "<<mainfile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------"<<endl;
        cout<<"Initializing Main input from file: "<<mainInfile<<endl;
        cout<<"-----------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------"<<endl;
        ePout<<"Initializing Main input from file: "<<mainInfile<<endl;
        ePout<<"-----------------------------------------------------------"<<endl;
    }
    cout<<endl;
    ePout<<endl;
    
    getline(infile, inputstring);
    cout<<inputstring<<endl;
    ePout<<inputstring<<endl;
    
    infile>>Ni>>Ndata>>pdfstyle;
    getline(infile, inputstring);
    cout<<"        "<<Ni<<"             "<<Ndata<<"              "<<pdfstyle<<endl;
    ePout<<"        "<<Ni<<"             "<<Ndata<<"              "<<pdfstyle<<endl;
    if ((pdfstyle=="C")||(pdfstyle=="c")) {
        pdftype=CTEQ;
    } else if ((pdfstyle=="L")||(pdfstyle=="l")) {
        pdftype=LHAPDF;
    } else {
        cerr<<"PDFtype must be C or L.  You input "<<pdfstyle<<endl;
        exit(1);
    }
    
    Na.resize(Ndata);
    Xmat.resize(Ndata);
		Xset.resize(Ndata);
    dX.resize(Ndata);
    X0.resize(Ndata);
    X0new.resize(Ndata);
    Q.resize(Ndata);
    Cm.resize(Ndata);
    dataset.resize(Ndata);

//
// Re-diagonalization Object:
//
    M.resize(Ni);
    for (int i=0;i<Ni;i++) {
        M[i].resize(Ni);
    }
    
    getline(infile, inputstring);
    cout<<inputstring<<endl;
    ePout<<inputstring<<endl;
    for (int i=0;i<Ndata;i++) {
        infile>>dataset[i]>>Na[i];
        cout<<"      "<<dataset[i]<<"            "<<Na[i]<<endl;
        ePout<<"      "<<dataset[i]<<"            "<<Na[i]<<endl;
        getline(infile, inputstring);
    }
    
    getline(infile, inputstring);
    cout<<inputstring<<endl;
    ePout<<inputstring<<endl;
    infile>>PDFin>>PDFout;
    cout<<"    "<<PDFin<<"        "<<PDFout<<endl;
    ePout<<"    "<<PDFin<<"        "<<PDFout<<endl;
    strcpy(PDFinfile,PDFin.c_str());
    strcpy(PDFoutfile,PDFout.c_str());
    cout<<endl;
    ePout<<endl;
    infile.close();
    
    for (int i=0;i<Ndata;i++) {
        strcpy(Theoryfile,dataset[i].c_str());
        strcat(Theoryfile,".theory");
        ReadInTheory(Theoryfile,i);
        SetQ(i);
        SetCinv_Max(i);
    }

}

void ePump::ReadInTolerances(void) {
    char Tolfile[80];
    ifstream infile;
    string inputstring;
    double Ni2;
    
    strcpy(Tolfile,PDFinfile);
    strcat(Tolfile,".tol");
    infile.open(Tolfile);
    if ( !infile.is_open() ) {
        cerr<<"Error opening tolerance file: "<<Tolfile<<endl;
        exit(1);
    } else {
        cout<<"-----------------------------------------------------------"<<endl;
        cout<<"Initializing Tolerance array for file: "<<Tolfile<<endl;
        cout<<"-----------------------------------------------------------"<<endl;
        ePout<<"-----------------------------------------------------------"<<endl;
        ePout<<"Initializing Tolerance array for file: "<<Tolfile<<endl;
        ePout<<"-----------------------------------------------------------"<<endl;
    }
    cout<<endl;
    ePout<<endl;
    
    do {
        getline(infile, inputstring);
        cout<<inputstring<<endl;
        ePout<<inputstring<<endl;
    } while ( inputstring[0]=='*');
    infile>>Ni2;
    if (Ni2!=Ni) {
        cerr<<"Number of EVs not consistent in *.in and *.tol files."<<endl;
        exit(1);
    }
 
    getline(infile, inputstring);
    cout<<"    "<<Ni<<endl;
    ePout<<"    "<<Ni<<endl;
    
    double TolSum=0;
//    double TdynSq1, TdynSq2;
    double temp;
		cout<<" i TsqP[i] TsqM[i] Tsqbar[i]"<<endl;
    for (int i=0;i<Ni;i++) {
        infile>>TdynP[i];
        getline(infile, inputstring);
        //cout<<fixed<<setprecision(2)<<TdynP[i]<<endl;
        //ePout<<fixed<<setprecision(2)<<TdynP[i]<<endl;
        infile>>TdynM[i];
        getline(infile, inputstring);
        //cout<<fixed<<setprecision(2)<<TdynM[i]<<endl;
        //ePout<<fixed<<setprecision(2)<<TdynM[i]<<endl;
        TdynRMS[i]=(TdynP[i]+TdynM[i])/2.0;
				cout<<setw(2)<<i<<"   "<<TdynP[i]<<"   "<<TdynM[i]<<"   "<<TdynRMS[i]<<endl;
				ePout<<setw(2)<<i<<"   "<<TdynP[i]<<"   "<<TdynM[i]<<"   "<<TdynRMS[i]<<endl;
        TolSum+=TdynRMS[i];
        TdynP[i]=sqrt(TdynP[i]);
        TdynM[i]=sqrt(TdynM[i]);
        TdynBar[i]=(TdynP[i]+TdynM[i])/2.0;
        TdynRMS[i]=sqrt(TdynRMS[i]);
    }
    cout<<" TolSum= "<<TolSum<<endl;
    ePout<<endl;
    
    TsqAvg=TolSum/double(Ni);
    cout<<" Average Squared Dynamical Tolerance, T^2 = "<<fixed<<setprecision(2)<<TsqAvg<<endl<<endl;
    ePout<<" Average Squared Dynamical Tolerance, T^2 = "<<fixed<<setprecision(2)<<TsqAvg<<endl<<endl;
    Tol=sqrt(TsqAvg);
//    cout<<"   i      (T_dyn[i])^2"<<endl;
//    ePout<<"   i      (T_dyn[i])^2"<<endl;
    cout<<" Dynamical tolerances "<<endl;
    ePout<<" Dynamical tolerances "<<endl;
    cout<<"   i      TdynP[i]    TdynM[i]    TdynBar[i]  TdynRMS[i]"<<endl;
    ePout<<"   i      TdynP[i]    TdynM[i]    TdynBar[i]  TdynRMS[i]"<<endl;
    for (int i=0;i<Ni;i++) {
        cout<<"  "<<setw(2)<<i+1<<setw(12)<<fixed<<setprecision(2)<<TdynP[i]<<setw(12)<<TdynM[i]<<setw(12)<<TdynBar[i]<<setw(12)<<TdynRMS[i]<<endl;
        ePout<<"  "<<setw(2)<<i+1<<setw(12)<<fixed<<setprecision(2)<<TdynP[i]<<setw(12)<<TdynM[i]<<setw(12)<<TdynBar[i]<<setw(12)<<TdynRMS[i]<<endl;
        TdynP[i]/=Tfix;
        TdynM[i]/=Tfix;
        TdynBar[i]/=Tfix;
        TdynRMS[i]/=Tfix;
    }
    cout<<endl;
    ePout<<endl;
    infile.close();
}

void ePump::ConvertTolerance(const char* chi2file, const char* tolfile) {
    
 // Converts file with Chi2f values to file with dynamic tolerance values for each error PDF
    
    ifstream infile;
    ofstream outfile;
    string inputstring;
    double Chi2f00;
    double Chi2f;
    
    infile.open(chi2file);
    if ( !infile.is_open() ) {
        cerr<<"Error opening Chi2F file: "<<chi2file<<endl;
        exit(1);
    } else {
        cout<<"Converting Chi2F file to Tolerance file."<<endl;
    }
    outfile.open(tolfile);
    
    
    do {
        getline(infile, inputstring);
        cout<<inputstring<<endl;
    } while ( inputstring[0]=='*');
    
    infile>>Ni;
    cout<<"    "<<Ni<<endl;
    getline(infile, inputstring);
    infile>>Chi2f00;
    cout<<Chi2f00<<endl;
    getline(infile, inputstring);
    for (int i=1;i<2*Ni+1;i++) {
        infile>>Chi2f;
        getline(infile, inputstring);
        outfile<<setprecision(12)<<Chi2f-Chi2f00<<endl;
        cout<<setprecision(12)<<Chi2f<<"     "<<Chi2f-Chi2f00<<endl;
    }
    getline(infile, inputstring);
    infile.close();
    outfile.close();
}







