module generalized_to_standard
  use mpi
  use descriptor_parameters
  use event_logger_m, only : add_event
  use processes, only : check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time
  implicit none

  private
  public :: reduce_generalized, recovery_generalized

contains

  subroutine reduce_generalized(dim, A, desc_A, B, desc_B)
    integer, intent(in) :: dim, desc_A(9), desc_B(9)
    double precision, intent(inout) :: A(:, :), B(:, :)

    integer :: info
    double precision :: scale, work_pdlaprnt(desc_B(block_row_))
    double precision :: time_start, time_end

    time_start = mpi_wtime()

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
      call terminate('reduce_generalized: pdpotrf failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('reduce_generalized:pdpotrf', time_end - time_start)
    time_start = time_end

    ! Reduction to standard problem by A <- L^(-1) * A * L'^(-1)
    call pdsygst(1, 'L', dim, A, 1, 1, desc_A, B, 1, 1, desc_B, scale, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdsygst): ", i0)', info
      call terminate('reduce_generalized: pdsygst failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('reduce_generalized:pdsygst', time_end - time_start)
    time_start = time_end
  end subroutine reduce_generalized


  subroutine recovery_generalized(dim, n_vec, B, desc_B, Vectors, desc_Vectors)
    integer, intent(in) :: dim, n_vec, desc_B(9), desc_Vectors(9)
    double precision, intent(in) :: B(:, :)
    double precision, intent(inout) :: Vectors(:, :)

    double precision :: time_start, time_end
    integer :: info

    time_start = mpi_wtime()

    ! Recovery eigenvectors by V <- L'^(-1) * V
    call pdtrtrs('L', 'T', 'N', dim, n_vec, B, 1, 1, desc_B, &
         Vectors, 1, 1, desc_Vectors, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdtrtrs): ", i0)', info
      call terminate('reduce_generalized: pdtrtrs failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('recovery_generalized:pdtrtrs', time_end - time_start)
  end subroutine recovery_generalized
end module generalized_to_standard
