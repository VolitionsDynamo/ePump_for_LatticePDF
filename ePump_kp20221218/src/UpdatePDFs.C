
#include "ePump.h"
#include <stdio.h>
#include <string.h>

static void show_usage(void)
{
    cerr<<"Usage: UpdatePDFs [Mainfile] [weight (optional)] \n"
        <<"where\n"
        <<"[Mainfile] indicates the name of the file with the input data, \"Mainfile.in\", \n"
        <<"and file to which the output is written, \"Mainfile.out\".\n"
        <<"[weight] is a universal weight factor by which to multiply the chi^2 from each of the data sets. \n"
        <<"(This is on top of the individual weight factors given to each of the data sets, as set in \"Mainfile.in\".)\n"
        <<"The default is weight=1.0.\n"
        <<"See file 01UpdatePDFs for more info."<<endl;
}

int main(int argc, char* argv[])
{
    PDF_format_type pdftype;
    double weight;
    char mainfile[80];
    
    if (!(argc==3||argc==2)) {
        show_usage();
        exit(1);
    } else if (argc==3) {
        weight=atof(argv[2]);
    } else {
        weight=1.0;
    }
    
    strcpy(mainfile,argv[1]);
		cout<<"mainfile: "<<mainfile<<endl;
    
// Set up the Error Update class:
    
    ePump EU(mainfile,weight);

//    EU.set_linear();
//    EU.set_DiagonalQuad();  // Include Diagonal Quadratic terms for Best-fit predictions.
//    EU.set_oldDefaults();   // Linear terms only for best fits and use TdyRMS[i].
//    EU.set_Pweight(0.0);
    EU.suppressCCs();
// If ".out" file is too large, one can suppress the reporting of Correlation Cosines by uncommenting the previous line, if desired.
//    EU.verbose();
// Uncomment the previous line to print out more information and checks on the numerical calculations
    EU.ReadInTheoryAndData();
    EU.ConstructUpdateMatrix();
    EU.UpdateObservables();
    
    return 0;
}

