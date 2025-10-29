
#include "ePump.h"
#include <stdio.h>
#include <string.h>

static void show_usage(void)
{
    cerr<<"Usage: ConvertToleranceFile [Infile] [Outfile]\n"
        <<"where\n"
        <<"[Infile] is the file with the list of Chi2F for each error PDF\n"
        <<"[Outfile] is the file with the list of dynamical tolerances"<<endl;
}

// Converts the CHI2F file to a Tolerance file used by epump.


int main(int argc, char* argv[])
{
    char Infile[80];
    char Outfile[80];
    
    if (!(argc==3)) {
        show_usage();
        exit(1);
    }
    
    strcpy(Infile,argv[1]);
    strcpy(Outfile,argv[2]);
    
// Set up the Error Update class:
    
    ePump EU("dum",100.0);
    EU.ConvertTolerance(Infile,Outfile);
    
    return 0;
}

