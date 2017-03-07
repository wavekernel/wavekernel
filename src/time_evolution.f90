module wk_time_evolution_m
  use mpi
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_processes_m
  use wk_event_logger_m
  use wk_global_variables_m
  implicit none

  private
  public :: step_forward

contains

  subroutine add_timer_event(name, wtime_start)
    character(*), intent(in) :: name
    real(8), intent(inout) :: wtime_start

    real(8) :: wtime_end

    if (check_master()) then
      wtime_end = mpi_wtime()
      call add_event(name, wtime_end - wtime_start)
      wtime_start = wtime_end
    end if
  end subroutine add_timer_event


  ! 時間ステップを1つ進める.
  ! Complexity: O(n^3).
  subroutine step_forward_crank_nicolson(H, H_desc, delta_t, dv_alpha, dv_alpha_next)
    complex(kind(0d0)), intent(inout) :: H(:, :)
    integer, intent(in) :: H_desc(desc_size)
    real(8), intent(in) :: delta_t
    complex(kind(0d0)), intent(in) :: dv_alpha(H_desc(rows_))
    complex(kind(0d0)), intent(out) :: dv_alpha_next(H_desc(rows_))

    complex(kind(0d0)) :: elem
    complex(kind(0d0)), allocatable :: dv_alpha_next_dist(:, :)
    integer :: dim, mb, i, ipiv_size, dv_desc(desc_size), info
    integer :: nprow, npcol, myrow, mycol, local_rows
    integer, allocatable :: ipiv(:)
    real(8) :: wtime_start
    integer :: numroc

    wtime_start = mpi_wtime()

    dim = H_desc(rows_)

    ! alpha_next = alpha initially.
    dv_alpha_next(1 : dim) = dv_alpha(1 : dim)

    call blacs_gridinfo(H_desc(context_), nprow, npcol, myrow, mycol)
    mb = H_desc(block_row_)
    call descinit(dv_desc, dim, 1, mb, H_desc(block_col_), &
         H_desc(rsrc_), H_desc(csrc_), H_desc(context_), dim, info)
    if (info /= 0) then
      call terminate('descinit failed in step_forward_crank_nicolson', info)
    end if
    local_rows = max(1, numroc(dim, mb, myrow, 0, nprow))
    allocate(dv_alpha_next_dist(local_rows, 1))
    call distribute_matrix_complex(H_desc(context_), 0, 0, dim, dv_alpha_next, &
       dv_desc, dv_alpha_next_dist)

    ! For pzgetrf.
    ipiv_size = numroc(dim, mb, myrow, H_desc(rsrc_), nprow) + mb
    allocate(ipiv(ipiv_size), stat=info)
    ipiv(:) = 0
    if (info /= 0) then
      call terminate('step_forward_crank_nicolson: allocation failed', info)
    end if

    call add_timer_event('step_forward_crank_nicolson:preparation', wtime_start)

    ! H' = L + U s.t. LU = (I + i H delta_t / 2) / 2, H <- H'
    H(:, :) = kImagUnit * delta_t / 4.0d0 * H(:, :)
    do i = 1, dim
      call pzelget('Self', ' ', elem, H, i, i, H_desc)
      elem = elem + (0.5d0, 0.0d0)
      call pzelset(H, i, i, H_desc, elem)
    end do

    call add_timer_event('step_forward_crank_nicolson:transform_A', wtime_start)

    call pzgetrf(dim, dim, H, 1, 1, H_desc, ipiv, info)
    if (info /= 0) then
      call terminate('pzgetrf failed', info)
    end if

    call add_timer_event('step_forward_crank_nicolson:pzgetrf', wtime_start)

    ! \alpha' = \frac{I - i H delta_t / 2}{I + i H delta_t / 2} \alpha
    !       = ((I + i H delta_t / 2) / 2)^{-1} \alpha - \alpha
    ! , alpha_next <- \alpha'
    call pzgetrs('No', dim, 1, &
         H, 1, 1, H_desc, ipiv, &
         dv_alpha_next_dist, 1, 1, dv_desc, info)
    deallocate(ipiv)
    if (info /= 0) then
      call terminate('step_forward_crank_nicolson: pzgetrs failed', info)
    end if
    call gather_matrix_complex(dv_desc(context_), dv_desc, dv_alpha_next_dist, 0, 0, dim, dv_alpha_next)
    call mpi_bcast(dv_alpha_next, dim, mpi_double_complex, 0, mpi_comm_world, info)
    if (info /= 0) then
      call terminate('step_forward_crank_nicolson: mpi_bcast failed', info)
    end if

    call add_timer_event('step_forward_crank_nicolson:pzgetrs', wtime_start)

    !call pzaxpy(dim, -kOne, &  ! alpha_next = alpha_next - alpha
    !     filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
    !     filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
    dv_alpha_next(1 : dim) = dv_alpha_next(1 : dim) - dv_alpha(1 : dim)

    call add_timer_event('step_forward_crank_nicolson:subtract_alpha', wtime_start)
  end subroutine step_forward_crank_nicolson


  !subroutine step_forward_taylor1(H, H_desc, delta_t, filtered_vecs, filtered_vecs_desc)
  !  complex(kind(0d0)), intent(inout) :: H(:, :), filtered_vecs(:, :)
  !  integer, intent(in) :: H_desc(desc_size), filtered_vecs_desc(desc_size)
  !  real(8), intent(in) :: delta_t
  !
  !  integer :: dim
  !
  !  dim = H_desc(rows_)
  !  ! alpha_next <- (I - i delta_t H1_eve) alpha.
  !  ! alpha_next = alpha initially.
  !  call pzcopy(dim, &
  !       filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !  call pzgemv('No', dim, dim, - kImagUnit * delta_t, &
  !      H, 1, 1, H_desc, &
  !      filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !      kOne, &
  !      filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !end subroutine step_forward_taylor1
  !
  !
  !subroutine step_forward_taylor2(H, H_desc, delta_t, filtered_vecs, filtered_vecs_desc)
  !  complex(kind(0d0)), intent(inout) :: H(:, :), filtered_vecs(:, :)
  !  integer, intent(in) :: H_desc(desc_size), filtered_vecs_desc(desc_size)
  !  real(8), intent(in) :: delta_t
  !
  !  integer :: dim
  !
  !  dim = H_desc(rows_)
  !
  !  ! alpha_next <- (I - i delta_t H1_eve - delta_t ^ 2 / 2 H1_eve ^ 2) alpha
  !  call pzcopy(dim, &
  !       filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !  call pzgemv('N', dim, dim, - kImagUnit * delta_t, &
  !       H, 1, 1, H_desc, &
  !       filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !       kZero, &
  !       filtered_vecs, 1, col_filtered_work, filtered_vecs_desc, 1)
  !  call pzgemv('N', dim, dim, - kImagUnit * delta_t / 2d0, &
  !       H, 1, 1, H_desc, &
  !       filtered_vecs, 1, col_filtered_work, filtered_vecs_desc, 1, &
  !       kZero, &
  !       filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1)
  !  call pzaxpy(dim, kOne, &
  !       filtered_vecs, 1, col_filtered_work, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !  call pzaxpy(dim, kOne, &
  !       filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !end subroutine step_forward_taylor2
  !
  !
  !subroutine step_forward_taylor3(H, H_desc, delta_t, filtered_vecs, filtered_vecs_desc)
  !  complex(kind(0d0)), intent(inout) :: H(:, :), filtered_vecs(:, :)
  !  integer, intent(in) :: H_desc(desc_size), filtered_vecs_desc(desc_size)
  !  real(8), intent(in) :: delta_t
  !
  !  integer :: dim
  !
  !  dim = H_desc(rows_)
  !
  !  call pzgemv('N', dim, dim, - kImagUnit * delta_t, &
  !       H, 1, 1, H_desc, &
  !       filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !       kZero, &
  !       filtered_vecs, 1, col_filtered_work, filtered_vecs_desc, 1)
  !  call pzgemv('N', dim, dim, - kImagUnit * delta_t / 2d0, &
  !       H, 1, 1, H_desc, &
  !       filtered_vecs, 1, col_filtered_work, filtered_vecs_desc, 1, &
  !       kZero, &
  !       filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1)
  !  call pzgemv('N', dim, dim, - kImagUnit * delta_t / 3d0, &
  !       H, 1, 1, H_desc, &
  !       filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1, &
  !       kZero, &
  !       filtered_vecs, 1, col_filtered_work2, filtered_vecs_desc, 1)
  !  call pzcopy(dim, &
  !       filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !  call pzaxpy(dim, kOne, &
  !       filtered_vecs, 1, col_filtered_work, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !  call pzaxpy(dim, kOne, &
  !       filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !  call pzaxpy(dim, kOne, &
  !       filtered_vecs, 1, col_filtered_work2, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_alpha_next, filtered_vecs_desc, 1)
  !end subroutine step_forward_taylor3


  subroutine step_forward(mode, H, H_desc, delta_t, dv_alpha, dv_alpha_next)
    character(*), intent(in) :: mode
    complex(kind(0d0)), intent(in) :: dv_alpha(:)
    complex(kind(0d0)), intent(inout) :: H(:, :)
    complex(kind(0d0)), intent(out) :: dv_alpha_next(:)
    integer, intent(in) :: H_desc(desc_size)
    real(8), intent(in) :: delta_t

    if (mode == 'crank_nicolson') then
      call step_forward_crank_nicolson(H, H_desc, delta_t, dv_alpha, dv_alpha_next)
    !else if (mode == 'taylor1') then
    !  call step_forward_taylor1(H, H_desc, delta_t, filtered_vecs, filtered_vecs_desc)
    !else if (mode == 'taylor2') then
    !  call step_forward_taylor2(H, H_desc, delta_t, filtered_vecs, filtered_vecs_desc)
    !else if (mode == 'taylor3') then
    !  call step_forward_taylor3(H, H_desc, delta_t, filtered_vecs, filtered_vecs_desc)
    else
      call terminate('not implemented', 61)
    end if
  end subroutine step_forward
end module wk_time_evolution_m
