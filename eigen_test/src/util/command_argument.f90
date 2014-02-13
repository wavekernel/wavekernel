module command_argument
  use processes, only : check_master, terminate
  implicit none

  ! Matrix Market format
  type matrix_info
    character(len=10) :: rep
    character(len=7) :: field
    character(len=19) :: symm
    integer :: rows, cols, entries
  end type matrix_info

  type argument
    character(len=256) :: matrix_A_filename = ''
    character(len=256) :: matrix_B_filename = '' ! Empty means standard eigenvalue problem
    type(matrix_info) :: matrix_A_info, matrix_B_info
    character(len=256) :: solver_type
    character(len=256) :: output_filename = 'eigenvalues.dat'
    logical :: is_generalized_problem
    logical :: is_printing_grid_mapping = .false.
    logical :: is_dry_run = .false.
    integer :: n_vec = -1, n_check_vec = -1 ! These default -1 mean 'all the vectors'
    ! When zero, orthogonality is not evaluated.
    integer :: ortho_check_index_start = 0
    integer :: ortho_check_index_end = 0
    character(len=256) :: eigenvector_dir = '.'
    integer :: printed_vecs_start = 0 ! Zero means do not print eigenvectors
    integer :: printed_vecs_end = 0
    integer :: verbose_level = 0
  end type argument

  private
  public :: matrix_info, argument, required_memory, &
       read_command_argument, print_command_argument

contains

  subroutine print_help()
    if (check_master()) then
      print *, 'Usage: eigen_test -s <solver_type> <options> <matrix_A> [<matrix_B>]'
      print *, 'Solver types are:'
      print *, '  scalapack_all (standard)'
      print *, '  scalapac_select (standard, selecting)'
      print *, '  general_scalapack_all (generalized)'
      print *, '  general_scalapack_select (generalized, selecting)'
      print *, 'Options are:'
      print *, '  -n <num>  (available with selecting solvers) Compute only &
           &<num> eigenpairs in ascending order of their eigenvalues'
      print *, '  -c <num>  Consider only <num> eigenvectors in residual norm checking'
      print *, '  -o <file>  Set output file name for eigenvalues to <file>'
      print *, '  -d <dir>  Set output files directory for eigenvectors to <dir>'
      print *, '  -p <num>  Specify the number of eigenvector to be output'
      print *, '  -p <num1>,<num2>  Specify range of the number of eigenvectors to be output'
      print *, '  -t <num1>,<num2>  Consider eigenvectors indexed <num1> to <num2>(included) in orthogonality checking'
      print *, '  -h  Print this help and exit'
      print *, '  --dry-run  Read command arguments and matrix files and instantly exit'
      print *, '  --print-grid-mapping  Print which process is assigned to each coordinate in BLACS grid'
      call flush(6)
    end if
  end subroutine print_help


  subroutine wrap_mminfo(filename, info)
    character(*), intent(in) :: filename
    type(matrix_info), intent(out) :: info

    integer, parameter :: iunit = 8

    open(unit=iunit, file=filename)
    call mminfo(iunit, info%rep, info%field, info%symm, &
         info%rows, info%cols, info%entries)
    close(iunit)
  end subroutine wrap_mminfo


  subroutine validate_argument(arg)
    type(argument), intent(in) :: arg

    integer :: dim
    logical :: is_size_valid, is_solver_valid, is_n_vec_valid

    ! Is matrix size appropriate?
    dim = arg%matrix_A_info%rows
    is_size_valid = dim == arg%matrix_A_info%cols
    if (arg%is_generalized_problem) then
      is_size_valid = is_size_valid &
           .and. (dim == arg%matrix_B_info%rows) &
           .and. (dim == arg%matrix_B_info%cols)
    end if
    if (.not. is_size_valid) then
      stop '[Error] validate_argument: Matrix dimension mismatch'
    end if

    ! Solver type and problem type matched?
    select case (trim(arg%solver_type))
    case ('lapack')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('scalapack_all')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('scalapack_select')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('general_scalapack_all')
      is_solver_valid = arg%is_generalized_problem
    case ('general_scalapack_select')
      is_solver_valid = arg%is_generalized_problem
    case ('eigenexa')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('general_eigenexa')
      is_solver_valid = arg%is_generalized_problem
    case default
      is_solver_valid = .false.
      stop '[Error] validate_argument: Unknown solver'
    end select
    if (.not. is_solver_valid) then
      if (arg%is_generalized_problem) then
        stop '[Error] validate_argument: This solver is not for generalized eigenvalue problem'
      else
        stop '[Error] validate_argument: This solver is not for standard eigenvalue problem'
      end if
    end if

    ! For all eigenpair computation, n_vec = n?
    select case (trim(arg%solver_type))
    case ('lapack')
      is_n_vec_valid = arg%n_vec == dim
    case ('scalapack_all')
      is_n_vec_valid = arg%n_vec == dim
    case ('general_scalapack_all')
      is_n_vec_valid = arg%n_vec == dim
    case default
      is_n_vec_valid = .true.
    end select
    if (.not. is_n_vec_valid) then
      stop '[Error] validate_argument: This solver does not support partial eigenvalue computation'
    end if

    if (arg%printed_vecs_start < 0 .or. arg%printed_vecs_end < 0 .or. &
         arg%printed_vecs_end > arg%n_vec .or. &
         arg%printed_vecs_start > arg%printed_vecs_end) then
      stop '[Error] validate_argument: Specified numbers with -p option are not valid'
    end if

    if (arg%ortho_check_index_start < 0 .or. arg%ortho_check_index_end < 0 .or. &
         arg%ortho_check_index_end > arg%n_vec .or. &
         arg%ortho_check_index_start > arg%ortho_check_index_end) then
      stop '[Error] validate_argument: Specified numbers with -t option are not valid'
    end if
  end subroutine validate_argument


  double precision function required_memory_lapack(arg)
    type(argument), intent(in) :: arg

    double precision :: num_double, dim

    num_double = real(arg%matrix_A_info%entries)
    dim = real(arg%matrix_A_info%rows)
    num_double = num_double + dim * dim
    required_memory_lapack = 8.0d0 * num_double
  end function required_memory_lapack


  double precision function required_memory_parallel_standard(arg)
    ! This is just an approximation and partial eigenvector computation
    ! (means reduced columns of eigenvector storage) is not supported yet
    ! Generalized version below has the same problem
    type(argument), intent(in) :: arg

    integer :: my_rank, n_procs
    double precision :: num_double, dim

    num_double = real(arg%matrix_A_info%entries)

    call blacs_pinfo(my_rank, n_procs)
    dim = real(arg%matrix_A_info%rows)
    ! 2 is for the input matrix and eigenvectors
    num_double = num_double + dim * dim * 2.0d0 / real(n_procs)

    required_memory_parallel_standard = 8.0d0 * num_double
  end function required_memory_parallel_standard


  double precision function required_memory_parallel_generalized(arg)
    type(argument), intent(in) :: arg

    integer :: my_rank, n_procs
    double precision :: num_double, dim

    num_double = real(arg%matrix_A_info%entries + arg%matrix_B_info%entries)

    call blacs_pinfo(my_rank, n_procs)
    dim = real(arg%matrix_A_info%rows)
    ! 3 is for the input matrices (A and B) and eigenvectors
    num_double = num_double + dim * dim * 3.0d0 / real(n_procs)

    required_memory_parallel_generalized = 8.0d0 * num_double
  end function required_memory_parallel_generalized


  double precision function required_memory(arg)
    type(argument), intent(in) :: arg

    select case (trim(arg%solver_type))
    case ('lapack')
      required_memory = required_memory_lapack(arg)
    case ('scalapack_all')
      required_memory = required_memory_parallel_standard(arg)
    case ('scalapack_select')
      required_memory = required_memory_parallel_standard(arg)
    case ('general_scalapack_all')
      required_memory = required_memory_parallel_generalized(arg)
    case ('general_scalapack_select')
      required_memory = required_memory_parallel_generalized(arg)
    case default
      required_memory = -1.0d0 ! Required memory unknown for this solver
    end select
  end function required_memory


  subroutine read_command_argument(arg)
    type(argument), intent(out) :: arg

    integer :: argi = 1
    character(len=256) :: arg_str
    logical :: exists
    integer :: i

    do while (argi <= command_argument_count())
      call get_command_argument(argi, arg_str)

      if (arg_str(1:1) == '-') then
        select case (trim(arg_str(2:)))
        case ('s')
          call get_command_argument(argi + 1, arg_str)
          arg%solver_type = trim(arg_str)
          argi = argi + 1
        case ('n')
          call get_command_argument(argi + 1, arg_str)
          read (arg_str, *) arg%n_vec
          argi = argi + 1
        case ('c')
          call get_command_argument(argi + 1, arg_str)
          read (arg_str, *) arg%n_check_vec
          argi = argi + 1
        case ('o')
          call get_command_argument(argi + 1, arg_str)
          arg%output_filename = trim(arg_str)
          argi = argi + 1
        case ('d')
          call get_command_argument(argi + 1, arg_str)
          arg%eigenvector_dir = trim(arg_str)
          argi = argi + 1
        case ('p')
          call get_command_argument(argi + 1, arg_str)
          i = index(arg_str, ',')
          if (i == 0) then
            read (arg_str, *) arg%printed_vecs_start
            arg%printed_vecs_end = arg%printed_vecs_start
          else
            read (arg_str(1 : i - 1), *) arg%printed_vecs_start
            read (arg_str(i + 1 :), *) arg%printed_vecs_end
          end if
          argi = argi + 1
        case ('t')
          call get_command_argument(argi + 1, arg_str)
          i = index(arg_str, ',')
          if (i == 0) then
            stop '[Error] read_command_argument: wrong format for -t option'
          else
            read (arg_str(1 : i - 1), *) arg%ortho_check_index_start
            read (arg_str(i + 1 :), *) arg%ortho_check_index_end
          end if
          argi = argi + 1
        case ('v')
          arg%verbose_level = 1
        case ('h')
          call print_help()
          stop ''
        case ('-dry-run')
          arg%is_dry_run = .true.
        case ('-print-grid-mapping')
          arg%is_printing_grid_mapping = .true.
        case default
          call print_help()
          stop '[Error] read_command_argument: unknown option'
        end select
      else if (len_trim(arg%matrix_A_filename) == 0) then
        ! The first non-option argument specifies the (left) input matrix
        arg%matrix_A_filename = trim(arg_str)
        ! Check whether the file exists
        inquire(file = trim(arg%matrix_A_filename), exist = exists)
        if (.not. exists) then
          stop '[Error] read_command_argument: Matrix A file not found'
        end if
      else
        arg%matrix_B_filename = trim(arg_str)
        inquire(file = arg%matrix_B_filename, exist = exists)
        if (.not. exists) then
          stop '[Error] read_command_argument: Matrix B file not found'
        end if
      end if
      argi = argi + 1
    enddo

    arg%is_generalized_problem = (len_trim(arg%matrix_B_filename) /= 0)

    call wrap_mminfo(arg%matrix_A_filename, arg%matrix_A_info)
    if (arg%is_generalized_problem) then
      call wrap_mminfo(arg%matrix_B_filename, arg%matrix_B_info)
    end if

    if (arg%n_vec == -1) then ! unspecified in command line arguments
      arg%n_vec = arg%matrix_A_info%rows
    end if

    if (arg%n_check_vec == -1) then
      arg%n_check_vec = arg%n_vec
    end if

    call validate_argument(arg)
  end subroutine read_command_argument


  subroutine print_matrix_info(name, info)
    character(*), intent(in) :: name
    type(matrix_info), intent(in) :: info

    print "('matrix ', A, ' field: ', A)", name, trim(info%field)
    print "('matrix ', A, ' symm: ', A)", name, trim(info%symm)
    print "('matrix ', A, ' rows: ', I0)", name, info%rows
    print "('matrix ', A, ' cols: ', I0)", name, info%cols
    print "('matrix ', A, ' entries: ', I0)", name, info%entries
  end subroutine print_matrix_info


  subroutine print_command_argument(arg)
    type(argument), intent(in) :: arg

    if (arg%is_generalized_problem) then
      print '("problem type: generalized")'
    else
      print '("problem type: standard")'
    end if

    print '("matrix A file: ", a)', trim(arg%matrix_A_filename)
    call print_matrix_info('A', arg%matrix_A_info)

    if (arg%is_generalized_problem) then
      print '("matrix B file: ", a)', trim(arg%matrix_B_filename)
      call print_matrix_info('B', arg%matrix_B_info)
    end if

    print '("solver: ", a)', trim(arg%solver_type)
    print '("eigenvalues output file: ", a)', trim(arg%output_filename)
    print '("required eigenpairs: ", i0)', arg%n_vec
    print '("verified eigenpairs: ", i0)', arg%n_check_vec
  end subroutine print_command_argument
end module command_argument
