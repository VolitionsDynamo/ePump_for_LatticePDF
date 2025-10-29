
#CPP=xlC
#CFLAGS=-O -lcomplex -lm
#CPP=clang++
CPP=g++
F77=f77
CFLAGS= -O3 
CFLAGS=
LIBS= -lm -lstdc++ 
CCFLAGS= -O3 -std=c++11 -c -Wformat-overflow=0 -fpermissive
#CCFLAGS= -c

.SUFFIXES: .o .C .f

BASE=ePump.o ReadIn.o Cholesky.o Eigen.o trim.o
BASEHEADERS=ePump.h VectorDef.h Cholesky.h Eigen.h

EXES= UpdatePDFs OptimizePDFs ConvertToleranceFile InterpolatePDFs

install: $(EXES)

UpdatePDFs: UpdatePDFs.o $(BASE) $(BASEHEADERS) 
	$(CPP) $(CFLAGS) -o UpdatePDFs UpdatePDFs.o $(BASE) $(LIBS)

OptimizePDFs: OptimizePDFs.o $(BASE) $(BASEHEADERS)
	$(CPP) $(CFLAGS) -o OptimizePDFs OptimizePDFs.o $(BASE) $(LIBS)

PseudoData: PseudoData.o $(BASE)random.o $(BASEHEADERS)
	$(CPP) $(CFLAGS) -o PseudoData PseudoData.o $(BASE)random.o $(LIBS)

ConvertToleranceFile: ConvertToleranceFile.o $(BASE) $(BASEHEADERS)
	$(CPP) $(CFLAGS) -o ConvertToleranceFile ConvertToleranceFile.o $(BASE) $(LIBS)

InterpolatePDFs: InterpolatePDFs.o $(BASE) $(BASEHEADERS)
	$(CPP) $(CFLAGS) -o InterpolatePDFs InterpolatePDFs.o $(BASE) $(LIBS)

.C.o:
	$(CPP) $(CCFLAGS) $< 
.f.o:
	$(F77) $(CCFLAGS) $< 

clean:
	-rm *.o $(EXES)


