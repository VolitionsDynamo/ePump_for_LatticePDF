
#include "ePump.h"
#include <stdio.h>
#include <string.h>

static void show_usage(void)
{
    cerr<<"Usage: InterpolatePDFs [PDFfile] \n"
        <<"where\n"
        <<"[PDFfile] is the name of the PDFs to be interpolated\n"
        <<"That is, PDFfile0.pds and PDFfile1.pds are given,\n"
        <<"and new interpolated PDFs, PDFfileZ.pds are produced,\n"
        <<"where f(Z) = f(0.0) + Z*(f(1.0)-f(0.0)) for\n"
        <<"Z=0.0 to 2.0 by 0.1 increments."<<endl;
}



int main(int argc, char* argv[])
{
    char PDFfile[80];
    
    if (!(argc==2)) {
        show_usage();
        exit(1);
    }
    
    strcpy(PDFfile,argv[1]);
    
// Set up the Error Update class:
    
    ePump EU("dum");
    EU.Interpolate(PDFfile);
    
    return 0;
}

