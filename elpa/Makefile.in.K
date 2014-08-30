FC = mpifrtpx
FFLAGS = -Kfast,openmp
LIBS = -I $(HOME)/include/elpa_openmp-2014.06.001/modules -L $(HOME)/lib -lelpa_openmp -SCALAPACK -SSL2BLAMP
