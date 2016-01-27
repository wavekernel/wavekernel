!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
      program Green_Lanc
      use common_constants
      use control
      use hamiltonian, only: read_hamiltonian, Ham, apply_hamiltonian
      use parallel_operation, only: init_workshare, pdprod
      implicit none
      integer :: CGSTEP, NE1, NE2, j, m
      type(control_parameters) :: input

      real(knd), allocatable :: b(:,:)
      real(knd) :: eps, err, dtmp
      real(knd),allocatable:: vec(:,:),alpha(:),beta(:)
      integer, allocatable :: yk(:)
      integer   :: ND, i, niter, np, id, lb, ub

*     reading parameters from file
      call read_parameters("input.txt",input)
      NE1=1
      NE2=1                     ! G(1:NE1, 1:NE2)
      CGSTEP=input%nitemax

*     reading Hamiltonian from file
c      call read_hamiltonian("Hamiltonian.dat")
                                ! binary data (Private format)
      call read_hamiltonian("Hamiltonian.txt") 
                   ! ascii data (See Dr. Takayama's documents)
*
      call set_list_of_orbitals_to_be_shown(input)
*     allocating prerequisites
      ND=Ham%dim
      allocate(b(ND,NE2))
      allocate(vec(ND,CGSTEP+1))
      allocate(alpha(CGSTEP-1))
      allocate(beta (CGSTEP))

      eps =input%eps_ref
      eps =sqrt(eps)            ! to keep backward compatibility
                                ! (With Takayama version)

      if(input%target_i /= 0) then
                                ! Compatible to Takayama         
         allocate(yk(NE1))
         do j=1, input%number_of_orbitals_to_be_shown
            m = input%list_of_orbitals_to_be_shown(j)
            yk(1)=m
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
            call init_workshare(ND,np,id,lb,ub)
            b(lb:ub,NE2)=0d0
C$OMP END PARALLEL
            b(m,NE2)=1d0
            call calc_Lanczos_coeff(b(:,1),CGSTEP,vec,alpha,beta,eps)
            if(CGSTEP.ne.input%nitemax) write(*,'("CGSTEP=",I6)')CGSTEP
            write(*,*)'min alpha=',minval(alpha(1:CGSTEP-1))
            write(*,*)'max alpha=',maxval(alpha(1:CGSTEP-1))
            write(*,*)'min beta =',minval(beta (1:CGSTEP))
            write(*,*)'max beta =',maxval(beta (1:CGSTEP))
            call write_green(j/=1) !Append_mode=(j/=1)
         end do
      else
                                ! Total DOS calc.
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(ND,np,id,lb,ub)
         call random_number(b(lb:ub,NE2)) ! it does not work on PGF90
         b(lb:ub,NE2)=b(lb:ub,NE2)-5d-1
         !because random_number is not thread safe on PGF90
C$OMP END PARALLEL
         dtmp=pdprod(b(:,NE2),b(:,NE2))
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(ND,np,id,lb,ub)
         b(lb:ub,NE2)=b(lb:ub,NE2)/dsqrt(dtmp)
C$OMP END PARALLEL
         call calc_Lanczos_coeff(b(:,1),CGSTEP,vec,alpha,beta,eps)
         if(CGSTEP.ne.input%nitemax) write(*,'("CGSTEP=",I6)') CGSTEP
         write(*,*)'min alpha=',minval(alpha(1:CGSTEP-1))
         write(*,*)'max alpha=',maxval(alpha(1:CGSTEP-1))
         write(*,*)'min beta =',minval(beta (1:CGSTEP))
         write(*,*)'max beta =',maxval(beta (1:CGSTEP))
         call write_green(.false.)
      end if


      contains

      subroutine calc_Lanczos_coeff(v1,M,v,a,b,eps)
      implicit none
      real(knd) :: v1(:),v(:,:),a(:),b(:),eps
      integer,intent(inout)::M
      intent(in) :: v1,eps
      intent(out):: v,a,b
      integer :: i, ND, np, id, lb, ub

      ND=size(v,1)

      if(size(v1) .ne. ND)
     $     stop 'calc_Lanczos : size mismatch (v1,v)'
      if(size(b).ne.size(a)+1)
     $     stop 'calc_Lanczos : size mismatch (a,b)'
      if(size(v,2).ne.size(b)+1)
     $     stop 'calc_Lanczos : size mismatch (v,b)'


C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(ND,np,id,lb,ub)
      v(lb:ub,1)=v1(lb:ub)
C$OMP END PARALLEL
      do i=1, M-1
         call apply_hamiltonian(ND,v(:,i),v(:,i+1))
         b(i)=pdprod(v(:,i),v(:,i+1))
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(ND,np,id,lb,ub)
         if(i.eq.1) then
            v(lb:ub,i+1)=v(lb:ub,i+1)-b(i)*v(lb:ub,i)
         else
            v(lb:ub,i+1)=v(lb:ub,i+1)-b(i)*v(lb:ub,i)
     $           -a(i-1)*v(lb:ub,i-1)
         end if
C$OMP END PARALLEL
         a(i)=sqrt(pdprod(v(:,i+1),v(:,i+1)))
         if( a(i) .le. eps ) exit
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(ND,np,id,lb,ub)
         v(lb:ub,i+1)=v(lb:ub,i+1)/a(i)
C$OMP END PARALLEL
      end do
      if(i.ne.M) then
         M=M-1
         write(*,'("calc_Lanczos_coeff: M=",I8)') M
      else
         call apply_hamiltonian(ND,v(:,i),v(:,i+1))
         b(i)=pdprod(v(:,i),v(:,i+1))
      end if
      end subroutine calc_Lanczos_coeff

      function continued_frac(z,M,b,a)
      complex(knd) :: continued_frac,result
      complex(knd),intent(in)::z
      integer, intent(in)::M
      real(knd),intent(in)::b(M),a(M-1)
      integer :: i
      result=z-b(M)
      do i=M-1, 1, -1
         result=z-b(i)-a(i)*a(i)/result
      end do
      continued_frac=1/result
      end function continued_frac
      

      subroutine write_green(append_mode)
*     write g_{jk}(\omega) to the file (for all \omega)
      use common_constants
      use file_io
      logical, intent(in):: append_mode
      integer :: IUNIT, IUNT2
      character(len=*), parameter :: FNAME="Green_Lanc.dat"
      character(len=*), parameter :: FNAM2="LDOS_Lanc.dat"
      character(len=8) :: MODE_STRING
      integer :: i, nmesh, j
      real(knd)    :: gamma, denergy, emin
      real(knd)    :: pi, integrated_ldos, ldos
      complex(knd) :: omega, resolvent
      pi=atan(1d0)*4

      nmesh   = input%nenemax
      denergy = input%denergy
      emin    = input%enmin
      gamma   = input%gamma

      if(append_mode)then
         MODE_STRING='OLD'
      else
         MODE_STRING='REPLACE'
      end if

      IUNIT=vacant_unit()
      open(IUNIT,FILE=FNAME,STATUS=trim(MODE_STRING),POSITION='APPEND')
      IUNT2=vacant_unit()
      open(IUNT2,FILE=FNAM2,STATUS=trim(MODE_STRING),POSITION='APPEND')

      integrated_ldos=0d0        ! initialize integrated_ldos

      do i=1, nmesh
         omega=cmplx(emin+(i-1)*denergy,gamma,knd)

         resolvent=0d0
         err      =0d0

         resolvent=continued_frac(omega,CGSTEP,beta,alpha)

         write(IUNIT,'(3ES21.13)') dble(omega),
     $        dble(resolvent),imag(resolvent)

         ldos= -imag(resolvent)/pi
         if(i==1 .or. i==nmesh)then
            integrated_ldos=integrated_ldos+ldos/2
         else
            integrated_ldos=integrated_ldos+ldos
         end if
         write(IUNT2,'(4ES21.13)') dble(omega), ldos,
     $        max(err/pi,9d-99), integrated_ldos*denergy
                  !max() here is a quick hack for gnuplot.
                  !gnuplot cannot read floting point without E or D,
                  !for example,  "1.0-108".
      end do
      write(IUNIT,*)
      write(IUNT2,*)
      close(IUNIT)
      close(IUNT2)
      end subroutine write_green

      end program Green_Lanc
