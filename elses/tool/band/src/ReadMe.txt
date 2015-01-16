On the BAND-structure calculation with GENO parameters

GENO : Work flow for band calculation
      1. Write Si_unit.xml
      2. Write Si_bulk_512atom.xml (only contains the size of resultant rectangular cell)
      3. Execute 
	../../tool/Perl/mkSupercell.pl Si_unit.xml Si_bulk_512atom_new.xml > mkSupercell.log
	(output : Si_bulk_512atom_new.xml & ForBandCalc.txt)
	(Si_bulk_512atom_new.xml contains every coordinates of atom in recutangular supercell)
      4. Execute
        ../../bin/elses -band > elses-band.log
        (output : OverlapAndHamiltonian.txt)
      5. ../../tool/band/src/band

Program 
      band.f90
      diag.f90 
      Makefile
      Lapack 

Input files : Output files of ELSES package 
    (1) ForBandCalc.txt
    Structure data file for the band structure calculation (Output files of ELSES package)
      # translation vectors of Bravais lattice (a.u.)
        6.741033718396681e+00   0.000000000000000e+00   0.000000000000000e+00
        0.000000000000000e+00   6.741033718396681e+00   0.000000000000000e+00
        0.000000000000000e+00   0.000000000000000e+00   6.741033718396681e+00
      # number of positions in the unitcell
         2
      # positions in the unitcell (internal co-ordinate)
        0.000000000000000e+00   0.000000000000000e+00   0.000000000000000e+00
        2.500000000000000e-01   2.500000000000000e-01   2.500000000000000e-01
      # number of unitcell
        32
      # translation vectors of the origin of unitcell
      # (internal co-ordinate, with shifted center)
        ..........................
      # k-th unitcell includes origin.
      # k is wrtten in the following line.
      29
    (2) OverlapAndHamiltonian.txt (Output files of ELSES package)
    Files of Hamiltonian, Overlap matrices, neighboring list, etc
      # size of Overlap & Hamiltonian storage S(int1,int2,int3,int4), H(int1,int2,int3,int4)
      # List of neighbor atoms (int3-th neighbor of int4)
      # Number of valence orbitals of int4-th atom
      # S(int1,int2,int3,int4) = double1, H(int1,int2,int3,int4) = double2
    (3) SymLine.txt  (Prepared by the user)
    Define the symmetric lines etc.
    -------<from here> a sample SymLine.txt-------
    6
    40   0.5    0.5    0.5    0.0    0.0    0.0       #L-G 
    40   0.0    0.0    0.0    0.0    1.0    0.0       #G-X
    20   0.0    1.0    0.0    0.5    1.0    0.0       #X-W
    30   0.5    1.0    0.0    0.5    0.5    0.5       #W-L
    20   0.5    0.5    0.5    0.0    0.75   0.75      #L-K
    40   0.0    0.75   0.75   0.0    0.0    0.0       #K-G
    -------<to here> a sample SymLine.txt-------
   # The first line contains the number of symmetric lines to be drawn.
   # Following 6 lines contain the specification of respective symmetry lines.
   # The format of specification is as follows.
   #	 number of points between b_start and b_end, b_start(x,y,z), b_end(x,y,z) 
   # All the first and the last points of symmetric lines should be given in a unit of 2*pi/aLat,
   # where aLat=|a_1| (a_i are the translation vectors of Bravais lattice (a.u.)

 
Output files:
    (1) EigenEnergy.txt
    #  kx ky kz K E1 E2 E3 ........
    # of kpoint, vector component of k-vector, k-distance from the first k-point, corresponding eigen enerfies (eV unit)   
    (2) GNU-Ekcurve-**** is one example for gnuplot routine. Here the file "EigenEnergy.txt" should be the input.
    K : horizontal coordinate 
    E1, E2, ... : vertical coordinmate
    The size of figures, the name of symmetric point, K-value of the symmetric point can be found in EigenEnergy.txt.  
    The resultant figure can be found in this directry.
