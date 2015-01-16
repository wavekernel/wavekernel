Usage: Calculating LDOS by using shifted COCG method.

Files:
		Makefile-elses	gnu make format, include "Makefile.inc" file
						 at the top directory of "ELSES" package
		Green_?.f (?=C,R)
						Main program for calculating Green's function of given Hamiltonian matrix,
						by using shifted COCG method.
						   (See R.Takayama, T.Hoshi, T.Sogabe, S.-L.Zhang, T.Fujiwara, PRB 73, 165108)
						?=R : Real version. Seed (reference energy) is restricted to be real quantity.
						?=C : Complex version.
						Module sCOCG is attached in Green_?.f

						Output files are "LDOS.dat" and "Green.dat". In the LDOS.dat,
						$\omega$'s and LDOS=$-\frac{1}{\pi}Im(G_{ii}(\omega))$'s, norms of the
						residual vectors, and the integrated LDOS are stored.
						In the "Green.dat", $\omega$'s, the real and imaginary part of
						$G(\omega)_{ii}$'s are stored.

		Green_SSW.f
						Main program for calculating Green's function of given Hamiltonian matrix,
						by using shifted COCG method with seed switching technique.

		Green_diag.f
						Main program for calculating Green's function of given Hamiltonian matrix,
						by using LAPACK. Hamiltonian is exactly diagonalized and calculate
						Green's function. This program is expected to give an exact solution
						(within the acccuracy of the LAPACK).
						The leading two lines of the file "input.txt" are ommitted. The following
						three lines define the energy mesh and indicate the element to be calculated.
						(See sample input file and following explanations bellow.)

						Output files are "LDOS_diag.dat", "TDOS_diag.dat" and "Green_diag.dat".
						In the "LDOS_diag.dat", $\omega$'s and LDOS=$-\frac{1}{\pi}Im(G_{ii}(\omega))$'s,
						9e-99, and the integrated LDOS are stored. The dummy value of
						9e-99 is due to keep the compatibility to Green_? (?=C,R).
						In the "TDOS_diag.dat", $\omega$'s and 
						total DOS=$-\frac{1}{\pi} trace G(\omega)$'s are stored.
						In the "Green.dat", $\omega$'s, the real and imaginary part of
						$G(\omega)_{ii}$'s are stored.
						
		Green_Lanc.f
						Main program for calculating Green's function of given Hamiltonian matrix,
						by using shifted Lanczos method.


		subroutine.f	Definitions of some modules used in "Green.f".
		diag.f			Wrappers of the LAPACK routines.


Perl script		
		residue.pl		Parsing standard output of the program Green and making the list of
						the norm of residual vector.
		residue.gnu		Command file to display the norm of the residual vectors.
						Try as follows.
						Green_C | residue.pl > log
						gnuplot
						(for the prompt of the gnuplot) iter=0
						(for the prompt of the gnuplot) load 'residue.gnu'

						!!Nov 7, 2007!! The above two scripts does not work currently, because the
						output format of Green_[C|R] is slightly changed.

Modules
		common_constants:		Physical constants, parameters related to the precision, etc.
		parallel_operations:	Basic tools for parallel operations, including parallelized dot_product.
				Notice! These routines are assuming OMP_DYNAMIC=FALSE.
		sort:					Merge sort routines. 
								input: N, value(N), output: index(N), value(index(:)) is sorted list
								(ascending order)
						imsort	sort list of integers (parallelized, but not so efficient)
						imsorts sort list of integers (serial version)
						lmsort	sort list of long int's (parallelized, but not so efficient)

		file_io:				Basic tools for file I/O.
						integer function vacant_unit() returns the vacant unit number.
						Example
						IUNIT =vacant_unit()
						IUNIT2=vacant_unit()
						OPEN(IUNIT, ... )
						OPEN(IUNIT2, ...)
						READ(IUNIT, ... )
						WRITE(IUNIT2,...)
						CLOSE(IUNIT)
						CLOSE(IUNIT2)

		sparse_matrix:			An implementation for treating sparse_matrix operation on a vector.
								Including type definition of "sparsemat"
						allocate_matrix(N,M,X)   allocate N x N sparse matrix X.
												 M denotes the number of non zero elements of X.
						parallel_setup(X)		 parallel setup for operating sparse matrix X on a vector.
						apply_matrix(X,vec,nvec) nvec = X vec

		control:				Control parameter set defined in Dr. Takayama's documents.
								And related routines.

		hamiltonian:			Module for treating Hamiltonian.
								Currently using the above "sparse_matrix" module.
								If you need to change the mechanism for treating the Hamiltonian,
								you are suggested to change this module.
								In that case, you must implement apply_hamiltonian and the definitions
								of the Hamiltonian.

						apply_hamiltonian(N,vec,nvec)
								Real(8) :: vec(N), nvec(N)
								nvec = H vec, where H denotes Hamiltonian.


		sCOCG:					Shifted COCG part of the program.
                                Depending on "hamiltonian", "parallel_operation", and "common_constants".

						apply_cocg		Doing COCG process and return alpa, beta, x(solution),
								   and more (dependent on the given arguments).
								   Calls apply_hamiltonian of the module "hamiltonian".
						apply_cg		Doing CG process and return alpa, beta, x(solution),
								   and more (dependent on the given arguments).
								   Calls apply_hamiltonian of the module "hamiltonian".

							Note:  "call apply_cg(nd,n,e,b,v)" is the most simple form of calling apply_cg,
								   and solve (e-H)x=b, where "H" denotes Hamiltoninan.
                            subroutine apply_cg(nd,n,e,b,v,alpha,beta,br,yr,y,yk,res,rnorm,xnorm,interval)
								On entry
								   nd: dimension of the matrix.
								   n:  iteration number
								   e:  reference energy
								   b:  initial vector
								   v:  work area v(nd,4)
								On exit
								   n:  actual number of iteration
								   v(:,1) : x (solution)
								   Other contents of v are destroyed.
								   Optional arguments
								   Imput
										res	  : Tolerance of the solution.
										   If the norm of the residual vector becomes lower than res before
										   iteration number reaches n, the iteration is discontinued.
										Interval: The interval of printing ||r_k|| on the screen.
								   Output
										alpha(n), beta(n) : store respective variables (alpha(k)=\alpha_{k-1})
										br(n)  : \bm{b}^T \cdot \bm{r}_{k-1}
										y(:,:) : Adjoint vectors
										yk(:)  : Alternative of y=0;y(yk(k),k)=1.
										rnorm(n): rnorm(k)  = ||\bm{r}_{k-1}||
										xnorm(n): xnorm     = ||\bm{x}_{m-1}||, (m=iteration number)

							function Gjk(sigma,m,j,k,a,b,yr,rnorm,er,interval,res); complex(16) :: Gjk
									 returns the \bm{y}_j^T (sigma+(Eref-H))^{-1} \bm{y}_k
									 sigma : energy shift (z-Eref)
									 j,k   : indeces of the elements of Green function to be solved.
									 a     : alpha
									 b     : beta
									 yr(j,i,k): y_j \cdot r^{(i)}_{k-1}
									       (initial vector = b_k, (i-1)-th iter.)
									 Optional
									    Imputs
										   res      : Tolerance of the solution.
										   interval : The interval of printing ||r^\sigma_k|| on the screen.
										   rnorm    : rnorm calculated in apply_cg or apply_cocg
										   renew    : if renew=T, \alpha, \beta, \bm{r} are renewed.
													(For seed switching)
										Outputs
									       er       : The error of approximate solution.
										         If "rnorm" is present in the list of argument,
												    er = ||r^\sigma_{m-1}||/pi_{m-1}
												      ( Not equals to ||r^\sigma_{m}||/pi_{m} ),
												 otherwise, er = |y^T r^\sigma_{m-1}|.
		*** Important notice *** "Gjk" in Green_R.f and that in Green_C.f are different.


						switch_seed(nd,m,j,k,z,sigma,alpha,beta,yr,res,rnorm,rvecs,yk,y,interval)
										On entry,
										   nd		: dimension of Hamiltonian
										   m		: length of alpha, beta, r. In other word, CGSTEP.
										   j, k, yk : 
												 input vector = (0,...,0,1[yk(k)-th column],0,...,0)^T
											   adjoint vector = (0,...,0,1[yk(j)-th column],0,...,0)
										   sigma    : energy shift
										   alpha    : \alpha_{i-1}
										   beta     :  \beta_{i-1}
										   y		: An alternative of the array "yk".
													 It gives adjoint vector, directly.
										   yr       : \bm{y}^T \bm{r}_{i-1}
										   res		: A criterion to stop recurrence Eq. in switch_seed.
										   rnorm    : ||\bm{r}_{i-1}||
										   rvecs    : rvecs(:,1)=\bm{r}_{m-2}
													: rvecs(:,2)=\bm{r}_{m-1}
													: rvecs(:,3)=\bm{r}_{m}
										On exit,
										   z		: new seed (reference energy) z+sigma
										   alpha	: new alpha
										   beta		: new beta
										   yr		: new yr

For calculating Green's function,
 call apply_cg or apply_cocg and calculate all the \alpha's, \beta's, and yr's (\bm{y}^T \bm{r}_k),
then call Gjk for every energy mesh point (sigma=(energy)-(reference energy)).

If reference energy is complex, you must use apply_cocg and complex version of Gjk.
You can use apply_cg, if and only if the reference energy is real.
The calculation time of apply_cg is nearly half of that of apply_cocg.
Usage of apply_cg may cause the more singular behavior of the result.
Therefore, take care when you call apply_cg.


Sample input file

  Essentially, the format of the input file is same as that
of the program written by Dr Takayama.
(H13_A04_soft01.tar.gz in 
http://act.jst.go.jp/content/h13/material/M04/PageMain.html
<the pages are written in Japanese>)
	  Differences between this version and Dr. Takayama's version are as follows.
				 1. This package does not treat energy scale conversion.
				 2. Input file is named as "input.txt".
				 3. Output file is named as "LDOS.dat" and "Green.dat".
					LDOS means the imaginary part of the Green's function with prefactor (-1/pi).
				 4. Format of LDOS.dat is as follows. First column shows energy, 
					second LDOS, third the norm of the residual vector, fourth integrated LDOS.
				 5. Format of Green.dat is "energy, real part, imaginary part" in each line.
			!!   For the program Green_diag, the string "_diag" is appended to the base name
			!! of the respective output file.
-------From here <input.txt>------
-16.0d0, -16.0d0, 500,    eps_ref,eps_shift,nitemax,
 0,-0.14d0,               iswEscale,Esample,
-0.60d0, 0.400d0, 1001,   enmin,enmax,nenemax,
 1, 1.00d-3, 2.0d0,       iswgamma,gamma,gammaXmesh,
   1,   1,   1,           target_i,target_j,target_k,
-------To here <input.txt>------
eps_ref	  : 10^{eps_ref} is a criterion of convergence at the seed (the reference energy).
eps_shift : Originally, 10^{eps_shift} is a criterion of convergence in shifted systems.
		  But now (Nov 11, 2007), this parameter is not used. "eps_ref" is used instead.
nitemax	  : The maximum number of COCG iterations at the seed (the reference energy).
iswEscale : (Reserved) Put 0, for the current version.
Esample	  : The seed (the reference energy).
enmin	  : The minimum end point of the energy mesh.
enmax	  : The maximum end point of the energy mesh.
nenemax	  : The number of the energy mesh points.
		  (Therefore, the interval of the energy mesh, denergy=(enmax-enmin)/(nenemax-1))
iswgamma,gamma,gammaXmesh  : Controls the way to determin the imaginary part of the energies
		  which smear the density of states.
		  If iswgamma=0, (the imaginary part of energies)=gamma.
		  If iswgamma=1, (the imaginary part of energies)=denergy*gammaXmesh.
target_i : Choose the elements of the Green's function to be calculated.
		 G_{target_i, target_i} (\omega) is calculated.
		 !! If target_i is equal to zero, the random vector v is generated
		 !! and <v|G(\omega)|v> is calculated.  This is useful to estimate
		 !! the total density of state.
target_j,target_k : (Reserved)



Format of the "Hamiltonian.txt"
----From here <a part of Hamiltonian.txt>----
    2048
    139264
    1  1  1   -0.15158976318922790000D-01              0.00000000000000000000D+00
    2  1  2   -0.59257816519425539000D-01              0.00000000000000000000D+00
----To here   <a part of Hamiltonian.txt> (following lines are omitted)----
	   The number in the first line shows the dimension of the Hamiltonian matrix.
	   The number in the second line shows the number of the non zero elements
		   in the Hamiltonian matrix.
	   The following lines show
		   n, i, j, Re(H_{ij}) and Im(H_{ij}), where "n" denotes the serial number of the elements,
		   and the H_{ij} denotes the (i,j) element of the Hamiltonian matrix and
		   Re(.) and Im(.) stands for the real and the imaginary part of the argument, respectively.


Sample output (on screen) of Green_diag
> cp input_ref input
> cp Hamiltonian_ref.txt Hamiltonian.txt
> Green_diag
----From here <sample output on screen>----
 #### gamma is replaced by denergy*gammaXmesh ####
 gamma=  2.000000000000000E-003
eps_ref,eps_shift,nitemax=   0.100000E-15  0.100000E-15  500
iswEscale,Esample=             0 -0.140000    
enmin,enmax,nenemax=        -0.600000      0.400000     1001
iswgamma,gamma,gammaXmesh=     1  0.200000E-02   2.00000    
target_i,target_j,target_k=    1    1    1
DIM=    2048
# of non-zero elements=    139264
check_hamiltonian:
   Dim   =                            2048
   Nelem =                          139264
   # of Ham%id1(i) > Ham%id2(i)=     68608
   (Ham%Nelem-Ham%dim)/2       =     68608
 set_NPROC: NPROC=           4
   symmetric=T
----To here   <sample output on screen>----



Update History:

*Aug. 4, 2008
The format of the files "Hamiltonian.txt" and "input.txt"
is extended in order to improve the usability of LDOS package.
Now we can use the atom and atomic-orbital pair in order to specify
the orbital that the density of state is projected onto.
A negative "target_i" triggers the extended mode to specify the
orbital. A sample input.txt is listed bellow.

-------From here, new way to specify the orbital in <input.txt>------
-16.0d0, -16.0d0,   500,   eps_ref,   eps_shift, nitemax
      0, -0.14d0,          iswEscale, Esample
  -50d0,    50d0,  1001,   enmin,     enmax,     nenemax
      0, 1.00d+0, 2.0d0,   iswgamma,  gamma,     gammaXmesh
     -1,       1,     1,   target_i,  target_j,  target_k
      4,                   Number of orbitals (for extended mode only)
      1,       1           (atom, atomic-orbital) pair 1
      1,       2           (atom, atomic-orbital) pair 2
      1,       3           (atom, atomic-orbital) pair 3
      1,       4           (atom, atomic-orbital) pair 4
-------To here, new way to specify the orbital in <input.txt>------

The local (projected) densities of state with respect to respective
orbitals are wrriten into "LDOS.dat", separately.
Blank line separates the respective densities of states, in "LDOS.dat".
