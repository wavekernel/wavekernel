module generalized_to_standard
  use descriptor_parameters
  use processes, only : check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time
  implicit none

  private
  public :: reduce_generalized, recovery_generalized

contains

  subroutine reduce_generalized(dim, A, desc_A, B, desc_B)
    integer, intent(in) :: dim, desc_A(9), desc_B(9)
    double precision, intent(inout) :: A(:, :), B(:, :)
    integer :: base_count
    double precision :: times(2)

    integer :: info
    double precision :: scale, work_pdlaprnt(desc_B(block_row_))

    call get_wall_clock_base_count(base_count)

    ! B = LL', overwritten to B
    call pdpotrf('L', dim, B, 1, 1, desc_B, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdpotrf): ", i0)', info
      if (info > 0) then
        info = min(info, 10)
        if (check_master()) print &
             '("The leading minor that is not positive definite (up to order 10) is:")'
        call eigentest_pdlaprnt(info, info, B, 1, 1, desc_B, 0, 0, '  B', 6, work_pdlaprnt)
      end if
      call terminate('[Error] reduce_generalized: pdpotrf failed')
    end if

    call get_wall_clock_time(base_count, times(1))

    ! Reduction to standard problem by A <- L^(-1) * A * L'^(-1)
    call pdsygst(1, 'L', dim, A, 1, 1, desc_A, B, 1, 1, desc_B, scale, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdsygst): ", i0)', info
      call terminate('[Error] reduce_generalized: pdsygst failed')
    end if

    call get_wall_clock_time(base_count, times(2))
    if (check_master()) then
      print *, '  reduce_generalized pdpotrf: ', times(1)
      print *, '  reduce_generalized pdsygst: ', times(2)
    end if
  end subroutine reduce_generalized


  subroutine recovery_generalized(dim, n_vec, B, desc_B, Vectors, desc_Vectors)
    integer, intent(in) :: dim, n_vec, desc_B(9), desc_Vectors(9)
    double precision, intent(in) :: B(:, :)
    double precision, intent(inout) :: Vectors(:, :)

    integer :: info

    ! Recovery eigenvectors by V <- L'^(-1) * V
    call pdtrtrs('L', 'T', 'N', dim, n_vec, B, 1, 1, desc_B, &
         Vectors, 1, 1, desc_Vectors, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdtrtrs): ", i0)', info
      call terminate('[Error] reduce_generalized: pdtrtrs failed')
    end if
  end subroutine recovery_generalized
end module generalized_to_standard
