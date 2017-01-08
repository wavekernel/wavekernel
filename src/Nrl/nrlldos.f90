!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! Module for calculating LDOS in a scheme with overlap matrix.
!
module MNrlLDos
  implicit none

  private

  ! public types
  public TNrlLDos
  public TNrlLDosCalculator

  ! public procedures
  public Initialize
  public Release
  public AssociateMatrices
  public SetSigma
  public GetLDos
  public GetPartialDos
  public GetCohp
  public WriteTo

  !!
  ! TNrlLDos stores the LDOS and related data.
  !
  type TNrlLDos
    integer :: count
    real(8) :: interval
    real(8) :: minimum_energy
    real(8) :: maximum_energy
    real(8), pointer :: nos(:)
    real(8), pointer :: dos(:)
  end type

  !!
  ! TNrlLDosCalculator stores the parameters required in LDOS calculations.
  !
  type TNrlLDosCalculator
    integer :: mode
    integer :: dimension
    real(8) :: sigma
    real(8), pointer :: non_orthogonal(:,:)
    real(8), pointer :: orthogonal(:,:)
    real(8), pointer :: energies(:)
    real(8), pointer :: overlap(:,:)
    real(8), pointer :: weights(:)
  end type

  ! Followings are public interfaces. 

  interface Initialize
    module procedure LDInitLDosCalculator
    module procedure LDInitLDos
  end interface

  interface Release
    module procedure LDReleaseLDosCalculator
    module procedure LDReleaseLDos
  end interface

  interface SetSigma
    module procedure LDSetSigma
  end interface

  interface WriteTo
    module procedure LDWriteTo
  end interface

  interface AssociateMatrices
    module procedure LDAssociateMatrices
  end interface

  interface GetLDos
    module procedure LDGetLDos
  end interface

  interface GetPartialDos
    module procedure LDGetPartialDos
  end interface

  interface GetCohp
    module procedure LDGetCohp
  end interface

  real(8), parameter :: M_PI = 3.141592653589793238462d0
    
contains

  !!
  ! Initialize procedure for TNrlLDos
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDInitLDos(io_ldos)
    type(TNrlLDos), intent(inout) :: io_ldos
    io_ldos%count = 0
    io_ldos%interval = 0d0
    io_ldos%minimum_energy = 0d0
    io_ldos%maximum_energy = 0d0
    nullify(io_ldos%nos)
    nullify(io_ldos%dos)
  end subroutine

  !!
  ! Finalize procedure for TNrlLDos
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDReleaseLDos(io_ldos)
    type(TNrlLDos), intent(inout) :: io_ldos
    call DeallocateLDosRange(io_ldos)
  end subroutine

  !!
  ! Release memory for TNrlLDos
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine DeallocateLDosRange(io_ldos)
    type(TNrlLDos), intent(inout) :: io_ldos
    if(associated(io_ldos%nos)) deallocate(io_ldos%nos)
    if(associated(io_ldos%dos)) deallocate(io_ldos%dos)
  end subroutine

  !!
  ! Set LDOS parameter
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetLDosRange(io_ldos, i_emin, i_emax, i_count, o_err)
    type(TNrlLDos), intent(inout) :: io_ldos
    real(8), intent(in) :: i_emin, i_emax
    integer, intent(in) :: i_count
    integer, intent(out) :: o_err

    call DeallocateLDosRange(io_ldos)
    allocate(io_ldos%nos(i_count), &
             io_ldos%dos(i_count), &
             stat=o_err)
    if(o_err /= 0) return

    io_ldos%count = i_count
    io_ldos%minimum_energy = i_emin
    io_ldos%maximum_energy = i_emax
  end subroutine     

  real(8) function GetStep(i_ldos) result(o_result)
    type(TNrlLDos), intent(in) :: i_ldos
    if (i_ldos%count > 1) then
      o_result = (i_ldos%maximum_energy - i_ldos%minimum_energy) &
                 / (i_ldos%count - 1)
    else
      o_result= 0d0
    end if
  end function

  !!
  ! Constructs DOS value by differentiate NOS.numerically.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ConstructDosFromNos(io_ldos)
    type(TNrlLDos), intent(inout) :: io_ldos
    integer :: count, e_pos
    real(8) :: de

    ! DOS calculation 
    count = io_ldos%count
    de = GetStep(io_ldos)
    io_ldos%dos(1) = 0d0
    io_ldos%dos(count) = 0d0
    do e_pos = 2, count - 1
      io_ldos%dos(e_pos) = (io_ldos%nos(e_pos + 1) - io_ldos%nos(e_pos - 1)) &
                           / (2 * de)
    end do
  end subroutine

  !!
  ! Write LDOS component in the indicated file.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDWriteTo(i_ldos, i_unit)
    type(TNrlLDos), intent(in) :: i_ldos
    integer, intent(in) :: i_unit

    integer :: i
    real(8) :: e, step

    step = GetStep(i_ldos)
    do i = 1, i_ldos%count
      e = i_ldos%minimum_energy + (i - 1) * step
      write(i_unit, "(3F18.12)") e, i_ldos%nos(i), i_ldos%dos(i) 
    end do
  end subroutine

  !!
  ! Initialize procedure for TNrlLDosCalculator
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDInitLDosCalculator(io_calculator)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator

    io_calculator%sigma = 1d-3
    io_calculator%dimension = 0

    nullify(io_calculator%non_orthogonal)
    nullify(io_calculator%orthogonal)
    nullify(io_calculator%overlap)
    nullify(io_calculator%energies)
    nullify(io_calculator%weights)
  end subroutine

  !!
  ! Finalize procedure for TTNrlLDosCalculator
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDReleaseLDosCalculator(io_calculator)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    call DeallocateLDosCalculator(io_calculator)
  end subroutine

  !!
  ! Release memory for TNrlLDosCalculator
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine DeallocateLDosCalculator(io_calculator)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    if(associated(io_calculator%orthogonal)) &
      deallocate(io_calculator%orthogonal)
    if(associated(io_calculator%weights)) &
      deallocate(io_calculator%weights)
  end subroutine

  !!
  ! Sets sigma of weight calculation
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDSetSigma(io_calculator, i_sigma)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    real(8), intent(in) :: i_sigma

    io_calculator%sigma = i_sigma
  end subroutine

  !!
  ! Set working matrices
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDAssociateMatrices(io_calculator, i_eigen, i_energies, &
      i_overlap, i_dimension, o_err)
    type(TNrlLdosCalculator), intent(inout) :: io_calculator
    real(8), intent(in), target :: i_eigen(:,:)
    real(8), intent(in), target :: i_energies(:)
    real(8), intent(in), target :: i_overlap(:,:)
    integer, intent(in) :: i_dimension
    integer, intent(out) :: o_err

    call DeallocateLDosCalculator(io_calculator)

    allocate(io_calculator%weights(i_dimension), &
             stat=o_err)
    if(o_err /= 0) return

    io_calculator%non_orthogonal => i_eigen
    io_calculator%energies => i_energies
    io_calculator%overlap => i_overlap
    io_calculator%dimension = i_dimension

    o_err = 0
  end subroutine

  !!
  ! Calculate integrated Lorentzian broadening.
  !
  real(8) function GetOccupation(i_energy, i_ecenter, i_sigma) &
      result(oresult)
    real(8) :: i_energy, i_ecenter, i_sigma
    oresult = 1.0d0 / M_PI * atan((i_energy - i_ecenter) / i_sigma) + 0.5d0
  end function

  !! 
  !  Calculates S^(1/2)
  !    where S is a overlap matrix
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SqrtMatrix(io_calculator, o_err)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    integer, intent(out) :: o_err

    integer :: dim, lwork, n, i, j
    real(8) :: e
    real(8), pointer :: matrix(:,:)
    real(8), allocatable :: buffer1(:,:), buffer2(:,:), eigen(:), work(:)

    dim = io_calculator%dimension
    matrix => io_calculator%overlap

    lwork = 3 * dim - 1
    allocate(buffer1(dim, dim), &
             buffer2(dim, dim), &
             eigen(dim), &
             work(lwork), stat=o_err)
    if(o_err /= 0) return

    ! get eigen value
    call dsyev('V', 'U',  dim, matrix(1, 1), dim, eigen(1), work(1), lwork, &
               o_err)

    ! b2 = e L_in
    if(o_err /= 0) return
    do n = 1, dim
      e = sqrt(eigen(n))
      do i = 1, dim
        buffer2(i, n) = matrix(i, n) * e
      end do
    end do
    
    ! b1 = L_in e L_jn
    call dgemm('N', 'T',  dim, dim, dim, &
               1.0d0, matrix(1,1), dim,  & 
                      buffer2(1,1), dim, &
               0.0d0, buffer1(1,1), dim)
    ! substitution
    matrix(1:dim, 1:dim) = buffer1(1:dim, 1:dim)
    deallocate(buffer1)
    deallocate(buffer2)
    deallocate(eigen)
    deallocate(work)

    o_err = 0
  end subroutine

  !! 
  !  Changes the basis from non-orthogonal to orthogonal
  !  Before calling this procedure, you must set S^(1/2) by calling
  !  FactorizeByCholesky.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine MultiplyMatrix(io_calculator)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    integer :: i, dim

    dim = io_calculator%dimension
    ! now, S^(1/2)
    call dsymm('L', 'U',  dim, dim, &
               1.0d0, io_calculator%overlap(1, 1), dim, & 
                      io_calculator%non_orthogonal(1, 1), dim, &
               0.0d0, io_calculator%orthogonal(1, 1), dim)
  end subroutine

  !!
  ! Converts basis from non-orthogonal to orthogonal
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ConvertBasisAsOrthogonal(io_calculator, o_err)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    integer, intent(out) :: o_err

    ! Calculates
    call SqrtMatrix(io_calculator, o_err)
    if(o_err /= 0) return
    call MultiplyMatrix(io_calculator)
    o_err = 0
  end subroutine

  !!
  ! Create the row list which has a weight for integration.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ConstructRowList(io_row_list, i_weights, i_dimension)
    integer, intent(inout) :: io_row_list(:)
    real(8), intent(in) :: i_weights(:)
    integer, intent(in) :: i_dimension
    integer :: i, pos

    pos = 1
    do i = 1, i_dimension
      if (abs(i_weights(i))  > 1d-12) then
        io_row_list(pos) = i
        pos = pos + 1
      end if
    end do
    io_row_list(pos) = i_dimension + 1
  end subroutine

  !!
  ! Calculate local density of states.
  ! After calculating LDOS, associated working matrix becomes undefined.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDGetLDos(io_calculator, i_emin, i_emax, i_count, &
      i_weights, io_ldos, o_err)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    integer, intent(in) :: i_count
    real(8), intent(in) :: i_emin
    real(8), intent(in) :: i_emax
    real(8), intent(in) :: i_weights(:)
    integer, intent(out) :: o_err
    type(TNrlLDos), intent(inout) :: io_ldos
 
    real(8) :: de, e, sum, occupation
    integer :: e_pos, n, i, j, dim, pos
    integer, allocatable :: row_list(:)
    real(8), allocatable :: dots(:)

    call SetLDosRange(io_ldos,i_emin, i_emax, i_count, o_err)
    if(o_err /= 0) return
 
    call DeallocateLDosCalculator(io_calculator)
    dim = io_calculator%dimension
    allocate(io_calculator%orthogonal(dim, dim), &
             row_list(dim + 1), &
             dots(dim),&
             stat=o_err)
    if(o_err /= 0) return

    call ConstructRowList(row_list, i_weights(:), dim)

    de  = GetStep(io_ldos)
    call ConvertBasisAsOrthogonal(io_calculator, o_err)
    if(o_err /= 0) return

    ! Calculates elements of local wavefunction
    do n = 1, dim
      sum = 0d0
      pos = 1
      do while(row_list(pos) <= dim) 
        j = row_list(pos)
        sum = sum + i_weights(j) * io_calculator%orthogonal(j, n) ** 2
        pos = pos + 1
      end do
      dots(n) = sum
    end do
    
    ! NOS calculation
    e = i_emin
    do e_pos = 1, i_count
      e = e + de
      sum = 0.0d0
      do n = 1, dim
        occupation = GetOccupation(e, io_calculator%energies(n), io_calculator%sigma) 
        if(occupation < 1d-12) cycle
        sum = sum + occupation * dots(n)
      end do
      ! "two" represents spin degeneracy
      io_ldos%nos(e_pos) = 2.0d0 * sum 
    end do

    ! DOS calculation 
    call ConstructDosFromNos(io_ldos)

    deallocate(row_list)
    deallocate(dots)
    o_err = 0
  end subroutine

  !!
  ! Calculate local density of states.
  ! After calculating LDOS, associated working matrix becomes undefined.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDGetPartialDos(io_calculator, i_emin, i_emax, i_count, &
      i_weights, io_ldos, o_err)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    integer, intent(in) :: i_count
    real(8), intent(in) :: i_emin
    real(8), intent(in) :: i_emax
    real(8), intent(in) :: i_weights(:)
    integer, intent(out) :: o_err
    type(TNrlLDos), intent(inout) :: io_ldos
 
    real(8) :: de, e, sum, occupation, value
    integer :: e_pos, n, i, j, dim, pos
    integer, allocatable :: row_list(:)
    real(8), allocatable :: dots(:)

    call SetLDosRange(io_ldos,i_emin, i_emax, i_count, o_err)
    if(o_err /= 0) return
 
    call DeallocateLDosCalculator(io_calculator)
    dim = io_calculator%dimension
    allocate(io_calculator%orthogonal(dim, dim), &
             row_list(dim + 1), &
             dots(dim), &
             stat=o_err)
    if(o_err /= 0) return

    call ConstructRowList(row_list, i_weights, dim)

    de  = GetStep(io_ldos)
    call ConvertBasisAsOrthogonal(io_calculator, o_err)
    if(o_err /= 0) return

    ! Calculates elements of local wavefunction
    do n = 1, dim
      dots(n) = 0d0
      sum = 0d0
      pos = 1
      do while(row_list(pos) <= dim) 
        j = row_list(pos)
        sum = sum + i_weights(j) * io_calculator%orthogonal(j, n)
        pos = pos + 1
      end do
      dots(n) = sum**2
    end do

    ! NOS calculation
    e = i_emin
    do e_pos = 1, i_count
      e = e + de
      sum = 0d0
      do i = 1, io_calculator%dimension
        occupation = GetOccupation(e, io_calculator%energies(i), io_calculator%sigma) 
        if(occupation < 1d-12) cycle

        sum = sum  + dots(i) * occupation
      end do
      ! "two" represents spin degeneracy
      io_ldos%nos(e_pos) = 2.0d0 * sum 
    end do
    ! DOS calculation 
    call ConstructDosFromNos(io_ldos)

    o_err = 0
    deallocate(row_list)
    deallocate(dots)
  end subroutine
  
  !!
  ! Calculate local density of states.
  ! After calculating LDOS, associated working matrix becomes undefined.
  ! This procedure is non-orthogonal basis based.
  ! 
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine LDGetCohp(io_calculator, i_emin, i_emax, i_count, &
      i_hamiltonian, i_weights1, i_weights2, io_ldos, o_err)
    type(TNrlLDosCalculator), intent(inout) :: io_calculator
    integer, intent(in) :: i_count
    real(8), intent(in) :: i_emin
    real(8), intent(in) :: i_emax
    real(8), intent(in) :: i_hamiltonian(:,:)
    real(8), intent(in) :: i_weights1(:)
    real(8), intent(in) :: i_weights2(:)
    integer, intent(out) :: o_err
    type(TNrlLDos), intent(inout) :: io_ldos
 
    real(8) :: de, e, sum, sum1, sum2, occupation, v1, v2, h
    real(8) :: weighted_value1, weighted_value2
    integer :: e_pos, n, i, j, dim, pos1, pos2
    integer, allocatable :: row_list1(:), row_list2(:)

    call SetLDosRange(io_ldos,i_emin, i_emax, i_count, o_err)
    if(o_err /= 0) return
 
    call DeallocateLDosCalculator(io_calculator)
    dim = io_calculator%dimension
    allocate(io_calculator%orthogonal(dim, dim), &
             row_list1(dim + 1), &
             row_list2(dim + 1), &
             stat=o_err)
    if(o_err /= 0) return

    call ConstructRowList(row_list1, i_weights1, dim)
    call ConstructRowList(row_list2, i_weights2, dim)

    de  = GetStep(io_ldos)

    ! COHP calculation
    e = i_emin
    do e_pos = 1, i_count
      e = e + de
      sum = 0d0
      do n = 1, dim
        occupation = &
          GetOccupation(e, io_calculator%energies(n), io_calculator%sigma) 
        if(occupation < 1d-12) cycle

        sum1 = 0d0
        pos1 = 1        
        do while(row_list1(pos1) <= dim)
          i = row_list1(pos1)
          v1 = io_calculator%non_orthogonal(i, n)

          sum2 = 0d0
          pos2 = 1
          do while(row_list2(pos2) <= dim)
            j = row_list2(pos2)

            v2 = io_calculator%non_orthogonal(j, n)
            h = i_hamiltonian(j, i)
            sum2 = sum2 + v1 * v2 * h * i_weights2(j)**2
            pos2 = pos2 + 1
          end do
          sum1 = sum1 + sum2 * i_weights1(i)**2
          pos1 = pos1 + 1
        end do
        sum = sum + sum1 * occupation
      end do
      ! "two" represents spin degeneracy
      io_ldos%nos(e_pos) = 2.0d0 * sum 
    end do
    ! DOS calculation 
    call ConstructDosFromNos(io_ldos)

    o_err = 0
    deallocate(row_list1)
    deallocate(row_list2)
  end subroutine
end module
