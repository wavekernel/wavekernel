FC = mpifrtpx
FFLAGS = -Kfast,parallel,openmp
LIBS = -I $(HOME)/lib -L $(HOME)/lib -lEigenExa -SCALAPACK -SSL2BLAMP
