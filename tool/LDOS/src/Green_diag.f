!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
      program Green_diag
      use common_constants
      use control
      use hamiltonian, only: read_hamiltonian, Ham
      use sparse_matrix, only: apply_matrix
      use parallel_operation, only: init_workshare, pdprod
      use diag
      implicit none
      integer :: CGSTEP, NE1, NE2
      type(control_parameters) :: input

      real(knd), allocatable :: b(:,:)
      real(knd) :: eps, err, dtmp
      real(knd),allocatable:: dham(:,:),eigv(:), wght(:)
      integer, allocatable :: yk(:)
      integer   :: ND, niter, np, id, lb, ub, j, k, m

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

      eps =input%eps_ref
      eps =sqrt(eps)            ! to keep backward compatibility
                                ! (With Takayama version)
      allocate(dham(ND,ND))
      allocate(eigv(ND))
      allocate(wght(ND))
      call sparse_to_dense
      call diag_sym(ND,dham,eigv)
      ! Now, 'dham' and 'eigv' store the eigenvectors and eigenvalues,
      ! respectively.

      if(input%target_i == 0) then
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
         call dgemv('T',ND,ND,1d0,dham,ND,b,1,0d0,wght,1)
         wght(:)=wght(:)*wght(:)
         call write_green(.false.) !Append_mode=.false.
      else
         do j=1, input%number_of_orbitals_to_be_shown
            m = input%list_of_orbitals_to_be_shown(j)
            do k=1, ND
               wght(k)=dham(m,k)*dham(m,k)
            end do
            call write_green(j/=1) ! Only for the first iteration,
                                   ! append_mode=.false. .
         end do
      end if

      contains

      subroutine write_green(append_mode)
*     write g_{jk}(\omega) to the file (for all \omega)
      use common_constants
      use file_io
      logical, intent(in) :: append_mode
      integer :: IUNIT, IUNT2, IUNT3, IUNT4
      character(len=*), parameter :: FNAME="Green_diag.dat"
      character(len=*), parameter :: FNAM2="LDOS_diag.dat"
      character(len=*), parameter :: FNAM3="TDOS_diag.dat"
      character(len=*), parameter :: FNAM4="EIGN_diag.dat"
      character(len=8) :: MODE_STRING
      integer :: i, nmesh, j
      real(knd)    :: gamma, denergy, emin
      real(knd)    :: pi, ldos, tdos, integrated_ldos, integrated_tdos
      real(knd)    :: ldos_p, tdos_p,  trapezoid_ldos,  trapezoid_tdos
      real(knd)    :: int_emin, int_a_x, w, or
      allocatable  :: int_emin(:)
      complex(knd) :: omega, resolvent, ctmp
      pi=atan(1d0)*4

      nmesh   = input%nenemax
      denergy = input%denergy
      emin    = input%enmin
      gamma   = input%gamma

      if(append_mode)then
         MODE_STRING='OLD'
      else
         MODE_STRING='REPLACE'
         IUNT4=vacant_unit()
         open(IUNT4,FILE=FNAM4,STATUS=trim(MODE_STRING))
         write(IUNT4,'(ES21.13)') eigv
         close(IUNT4)
      end if
      
      IUNIT=vacant_unit()
      open(IUNIT,FILE=FNAME,STATUS=trim(MODE_STRING),POSITION='APPEND')
      IUNT2=vacant_unit()
      open(IUNT2,FILE=FNAM2,STATUS=trim(MODE_STRING),POSITION='APPEND')
      IUNT3=vacant_unit()
      open(IUNT3,FILE=FNAM3,STATUS=trim(MODE_STRING),POSITION='APPEND')

      allocate(int_emin(nd))
      int_emin=0d0
      err     =0d0
      do j=1, nd
         err=err+wght(j)
         omega=cmplx(emin,gamma,knd)
         int_emin(j)=atan((emin-eigv(j))/gamma)
      end do
      err=err-1d0
      if(abs(err) > 1d-9)then
         write(*,'("write_green: wrong normalization")')
         stop
      end if

      trapezoid_ldos=0d0
      trapezoid_tdos=0d0
      do i=1, nmesh
         or   =emin+(i-1)*denergy
         omega=cmplx(or,gamma,knd)
         resolvent      =0d0
         integrated_ldos=0d0
         integrated_tdos=0d0
         ctmp           =0d0
         do j=1, nd
            w        =wght(j)
            resolvent=resolvent+w/(omega-eigv(j))
            ctmp     =ctmp+   1d0/(omega-eigv(j))
            int_a_x=atan((or-eigv(j))/gamma)-int_emin(j)
            int_a_x=int_a_x/pi
            integrated_tdos=integrated_tdos+int_a_x   ! exact
            integrated_ldos=integrated_ldos+int_a_x*w
         end do
         write(IUNIT,'(3ES21.13)') or,dble(resolvent),imag(resolvent)

         ldos= -imag(resolvent)/pi
         tdos= -imag(ctmp)/pi
         if(i .ne. 1) then
            trapezoid_ldos=trapezoid_ldos+(ldos+ldos_p)*denergy/2
            trapezoid_tdos=trapezoid_tdos+(tdos+tdos_p)*denergy/2
         end if
         write(IUNT2,'(4ES21.13)') or,ldos,
     $        abs(trapezoid_ldos-integrated_ldos),
     $        integrated_ldos
                  !gnuplot cannot read floting point without E or D,
                  !for example,  "1.0-108".
         write(IUNT3,'(4ES21.13)') or, tdos,
     $        abs(trapezoid_tdos-integrated_tdos),
     $        integrated_tdos
         ldos_p=ldos
         tdos_p=tdos
      end do
      write(IUNT2,*)
      write(IUNT3,*)
      close(IUNIT)
      close(IUNT2)
      close(IUNT3)
      deallocate(int_emin)
      end subroutine write_green

      subroutine sparse_to_dense
      integer i, j, k
      dham=0d0
      do i=1, Ham%nelem
         j=Ham%id1(i)
         k=Ham%id2(i)
         dham(j,k)=dham(j,k)+Ham%val(i)
      end do
      end subroutine sparse_to_dense

      subroutine test_ham
      use sparse_matrix
      integer :: i,j
      real(knd),allocatable   ::A(:,:)
      !complex(knd),allocatable::EIG(:)
      real(knd),allocatable::EIG(:),vec(:)
      allocate(A(nd,nd))
      allocate(EIG(nd))
      allocate(vec(nd))
      A=0d0
      do i=1, Ham%nelem
         A(Ham%id1(i),Ham%id2(i))=A(Ham%id1(i),Ham%id2(i))+Ham%val(i)
c         A(Ham%id1(i),Ham%id2(i))=Ham%val(i)
      end do
c      call diag_gen(nd,A,eig)
      call diag_sym(nd,A,eig)
c      write(*,*) eig
      do i=1, Ham%dim
         call apply_matrix(Ham,A(:,i),vec)
         write(*,*) dot_product(A(:,i),vec)
         vec=vec-eig(i)*A(:,i)
         if( sqrt(dot_product(vec,vec)) .gt. 1d-10 ) then
            write(*,*) i, eig(i), A(:,i)
         end if
      end do
      deallocate(A)
      deallocate(EIG)
      deallocate(vec)
      end subroutine test_ham
      end program Green_diag
