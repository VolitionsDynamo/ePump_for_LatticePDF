
#include "ePump.h"

static void show_usage(void)
{
    cerr<<"Usage: OptimizePDFs [Mainfile]\n"
        <<"where\n"
        <<"[Mainfile] indicates the name of the file with the input data, \"Mainfile.in\", \n"
        <<"and file to which the output is written, \"Mainfile.out\".\n"
        <<"See file 02OptimizePDFs for more info."<<endl;
}

int main(int argc, char* argv[])
{
    char mainfile[80];
    
    if (!(argc==2)) {
        show_usage();
        exit(1);
    }
    
    strcpy(mainfile,argv[1]);
    
// Set up the Error Update class:
    
    ePump EU(mainfile);
//    EU.verbose();
// Uncomment the previous line to print out more information and checks on the numerical calculations
    EU.ReadInTheory_Optimize();
    EU.ConstructOptimizeMatrix();
    EU.OptimizeObservables();    
    
    return 0;
}

