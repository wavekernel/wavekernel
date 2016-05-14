FC = mpifrtpx
CPPFLAG        =  -Cpp
MACRO_D_PREFIX =  -D
FCFLAGS = -Kfast,openmp
LIBS = -I$(HOME)/lib/EigenExa2.3c/ -L$(HOME)/lib/EigenExa2.3c/ -lEigenExa -I$(HOME)/lib/elpa_openmp-2015.05.001/include/elpa_openmp-2015.05.001/modules/ -L$(HOME)/lib/elpa_openmp-2015.05.001/lib/ -lelpa_openmp -SCALAPACK -SSL2BLAMP
WNO_UNUSED_DUMMY_ARGUMENT_FLAG = 
