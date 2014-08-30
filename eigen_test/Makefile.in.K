FC = mpifrtpx
FFLAGS = -Kfast,openmp
LIBS = -I $(HOME)/lib -L $(HOME)/lib -lEigenExa_timer -I $(HOME)/include/elpa_openmp-2014.06.001/modules -lelpa_openmp -SCALAPACK -SSL2BLAMP