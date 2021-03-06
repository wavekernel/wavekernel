module wk_util_m
  implicit none

  private
  public :: wk_energy_t, wk_error_t, index_to_coordinate, truncate_imag, &
       add_number_to_filename, add_numbers_to_filename, add_numbers_to_filename_multiple_initial, &
       add_postfix_to_filename, &
       remove_directory_from_filename, read_vector, wk_random_number, comb_sort, &
       find_range_index

  type wk_energy_t
    real(8) :: tightbinding, nonlinear, total, tightbinding_deviation
  end type wk_energy_t

  type wk_error_t
    real(8) :: absolute, relative
  end type wk_error_t

contains

  ! i = 1 .. dim を x \in [0, 1] に等分割
  double precision function index_to_coordinate(dim, i)
    integer :: dim, i

    index_to_coordinate = dble(i - 1) / dble(dim - 1)
  end function index_to_coordinate


  double precision function truncate_imag(z)
    complex(kind(0d0)) :: z
    if (abs(imag(z)) > 1.0d-10) then
      write(0, *) 'Warn: non-negligible imaginary part is truncated'
      write(0, *) 'z: ', z
    end if
    truncate_imag = dble(z)
  end function truncate_imag


  function add_numbers_to_filename(filename_base, n_begin, n_end) result(filename)
    character(*), intent(in) :: filename_base
    integer, intent(in) :: n_begin, n_end
    character(len=len_trim(filename_base)+14) :: filename
    character(len=6) :: n_begin_str, n_end_str
    integer :: index_last_dot

    write(n_begin_str, '(I6.6)') n_begin
    write(n_end_str, '(I6.6)') n_end
    index_last_dot = index(filename_base, '.', back=.true.)

    if (index_last_dot == 0) then
      filename = filename_base // '_' // n_begin_str // '-' // n_end_str
    else
      filename = filename_base(1 : index_last_dot - 1) // '_' // n_begin_str // &
           '-' // n_end_str // filename_base(index_last_dot : len_trim(filename_base))
    end if
  end function add_numbers_to_filename


  function add_number_to_filename(filename_base, n) result(filename)
    character(*), intent(in) :: filename_base
    integer, intent(in) :: n
    character(len=len_trim(filename_base)+7) :: filename
    character(len=6) :: n_str
    integer :: index_last_dot

    write(n_str, '(I6.6)') n
    index_last_dot = index(filename_base, '.', back=.true.)

    if (index_last_dot == 0) then
      filename = filename_base // '_' // n_str
    else
      filename = filename_base(1 : index_last_dot - 1) // '_' // n_str // &
           filename_base(index_last_dot : len_trim(filename_base))
    end if
  end function add_number_to_filename


  function add_numbers_to_filename_multiple_initial(filename_base, multiple_initial_index, n_begin, n_end) result(filename)
    character(*), intent(in) :: filename_base
    integer, intent(in) :: multiple_initial_index, n_begin, n_end
    character(len=len_trim(filename_base)+21) :: filename
    character(len=6) :: multiple_initial_index_str, n_begin_str, n_end_str
    integer :: index_last_dot

    write(multiple_initial_index_str, '(I6.6)') multiple_initial_index
    write(n_begin_str, '(I6.6)') n_begin
    write(n_end_str, '(I6.6)') n_end
    index_last_dot = index(filename_base, '.', back=.true.)

    if (index_last_dot == 0) then
      filename = filename_base // '_' // multiple_initial_index_str // '_' // n_begin_str // '-' // n_end_str
    else
      filename = filename_base(1 : index_last_dot - 1) // '_' // multiple_initial_index_str // '_' // n_begin_str // &
           '-' // n_end_str // filename_base(index_last_dot : len_trim(filename_base))
    end if
  end function add_numbers_to_filename_multiple_initial


  function add_postfix_to_filename(filename_base, postfix) result(filename)
    character(*), intent(in) :: filename_base, postfix
    character(len=len_trim(filename_base)+len_trim(postfix)) :: filename
    integer :: index_last_dot

    index_last_dot = index(filename_base, '.', back=.true.)

    if (index_last_dot == 0) then
      filename = filename_base // postfix
    else
      filename = filename_base(1 : index_last_dot - 1) // postfix // &
           filename_base(index_last_dot : len_trim(filename_base))
    end if
  end function add_postfix_to_filename


  function remove_directory_from_filename(filename_base) result(filename)
    character(*), intent(in) :: filename_base
    character(len=len_trim(filename_base)) :: filename
    integer :: index_last_slash

    index_last_slash = index(filename_base, '/', back=.true.)
    filename = filename_base(index_last_slash + 1 : )
  end function remove_directory_from_filename


!  subroutine read_matrix(filename, matrix)
!    character(*), intent(in) :: filename
!    complex(kind(0d0)), allocatable, intent(out) :: matrix(:, :)
!
!    integer, parameter :: iunit = 8
!    integer :: rows, cols, nnz, info, i
!    integer, allocatable :: indx(:), jndx(:)
!    double precision, allocatable :: rval(:)
!    character :: rep*10, field*7, symm*19
!
!    open(iunit, file = filename, status='old')
!    call mminfo(iunit, rep, field, symm, rows, cols, nnz)
!    !print *, rep, field, symm, rows, cols, nnz
!    if (trim(rep) /= 'coordinate' .or. trim(field) /= 'real' .or. &
!         trim(symm) /= 'symmetric' .or. rows /= cols) then
!      stop 'invalid matrix file'
!    end if
!    rewind(iunit)
!    allocate(indx(nnz), jndx(nnz), rval(nnz), matrix(rows, cols), stat=info)
!    if (info /= 0) then
!      stop 'memory allocation failed'
!    end if
!    call mmread(iunit, rep, field, symm, rows, cols, nnz, &
!         nnz, indx, jndx, 0, rval, 0)
!    close(iunit)
!
!    matrix(:, :) = (0.0d0, 0.0d0)
!    do i = 1, nnz
!      matrix(indx(i), jndx(i)) = rval(i)  ! Cast occurs
!      matrix(jndx(i), indx(i)) = rval(i)
!    end do
!  end subroutine read_matrix


  subroutine read_vector(dim, filename, vector)
    integer, intent(in) :: dim
    character(*), intent(in) :: filename
    complex(kind(0d0)), intent(out) :: vector(dim)

    integer :: i, j
    integer, parameter :: iunit = 13
    double precision :: re, im

    vector(:) = (0d0, 0d0)
    open(iunit, file=trim(filename), status='old')
    do i = 1, dim
      read(iunit, *) j, re, im
      if (i /= j) then
        write(0, *) 'Warn: malformed vector file'
        write(0, *) 'filename, line: ', filename, j
      end if
      vector(j) = cmplx(re, im, kind(0d0))
    end do
    close(iunit)
  end subroutine read_vector


  subroutine wk_random_number(seed, val)
    integer, intent(inout) :: seed
    real(8), intent(out) :: val
    integer, parameter :: a = 1029, c = 555, n = 65536

    seed = mod(seed * a + c, n)
    val = dble(seed) / dble(n)
  end subroutine wk_random_number


  real(8) function get_charge_factor(atom_i, atom_elements, charge_factor_common, charge_factor_H, charge_factor_C)
    integer, intent(in) :: atom_i
    character, intent(in) :: atom_elements(:)
    real(8), intent(in) :: charge_factor_common, charge_factor_H, charge_factor_C

    if (atom_elements(atom_i) == 'H' .and. charge_factor_H >= 0d0) then
      get_charge_factor = charge_factor_H
    else if (atom_elements(atom_i) == 'C' .and. charge_factor_C >= 0d0) then
      get_charge_factor = charge_factor_C
    else  ! TODO: add oxygen.
      get_charge_factor = charge_factor_common
    end if
  end function get_charge_factor


  subroutine comb_sort(n, xs, js)  ! Sort xs in ***descending*** order.
    integer, intent(in) :: n
    real(8), intent(inout) :: xs(n)
    integer, intent(out) :: js(n)

    integer :: h, i, j_temp
    real(8) :: x_temp
    logical :: is_swapped

    h = n
    is_swapped = .false.

    do i = 1, n
      js(i) = i
    end do

    do while (h > 1 .or. is_swapped)
      if (h > 1) then
        h = (h * 10) / 13
      end if
      is_swapped = .false.
      do i = 1, n - h
        if (xs(i) < xs(i + h)) then
          x_temp = xs(i)
          j_temp = js(i)
          xs(i) = xs(i + h)
          js(i) = js(i + h)
          xs(i + h) = x_temp
          js(i + h) = j_temp
          is_swapped = .true.
        end if
      end do
    end do
  end subroutine comb_sort


  ! Find index of the range that contains given element.
  ! when ranges = [1, 3, 6, 11], it expresses list of successive ranges
  ! 1 <= x < 3, 3 <= x < 6, 6 <= x < 11.
  ! Searched element x is 5, returned range index is 2 (3 <= x < 6).
  integer function find_range_index(ranges, x)
    integer, intent(in) :: ranges(:), x

    integer :: n, s, e, m, i

    n = size(ranges, 1) - 1
    if (n < 1) then
      stop 'find_range_index: input does not express ranges '
    end if

    ! Bisection method.
    s = 1
    e = n
    do while (e - s > 7)
      m = (s + e) / 2
      if (x < ranges(m)) then
        e = m - 1
      else
        s = m
      end if
    end do

    ! When the width of range reaches 8, do sequential search.
    do i = s, e
      if (ranges(i) <= x .and. x < ranges(i + 1)) then
        find_range_index = i
        return
      end if
    end do
    find_range_index = 0  ! Not found.
  end function find_range_index
end module wk_util_m
