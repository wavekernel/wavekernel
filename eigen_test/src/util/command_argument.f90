module command_argument
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
    integer :: n_vec = -1, n_check_vec = -1 ! These default -1 mean 'all the vectors'
    integer :: verbose_level = 0
  end type argument

  private
  public :: matrix_info, argument, required_memory, &
       read_command_argument, print_command_argument

contains
  subroutine print_help()
    print *, 'help not written'
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


  integer function validate_argument(arg)
    type(argument), intent(in) :: arg

    ! for all eigenpair computation, n_vec = n?
    ! file exists?
    ! square matrix?
    validate_argument = 0 ! Not implemented yet
  end function validate_argument


  integer function required_memory_lapack(arg)
    type(argument), intent(in) :: arg

    integer :: dim, num_double = 0

    num_double = arg%matrix_A_info%entries
    dim = arg%matrix_A_info%rows
    num_double = num_double + dim * dim
    required_memory_lapack = 8 * num_double
  end function required_memory_lapack


  integer function required_memory_parallel_standard(arg)
    ! This is just an approximation and partial eigenvector computation
    ! (means reduced columns of eigenvector storage) is not supported yet
    ! Generalized version below has the same problem
    type(argument), intent(in) :: arg

    integer :: my_rank, n_procs, dim, num_double = 0

    num_double = arg%matrix_A_info%entries

    call blacs_pinfo(my_rank, n_procs)
    dim = arg%matrix_A_info%rows
    ! 2 is for the input matrix and eigenvectors
    num_double = num_double + dim * dim * 2 / n_procs

    required_memory_parallel_standard = 8 * num_double
  end function required_memory_parallel_standard


  integer function required_memory_parallel_generalized(arg)
    type(argument), intent(in) :: arg

    integer :: my_rank, n_procs, dim, num_double = 0

    num_double = arg%matrix_A_info%entries + arg%matrix_B_info%entries

    call blacs_pinfo(my_rank, n_procs)
    dim = arg%matrix_A_info%rows
    ! 3 is for the input matrices (A and B) and eigenvectors
    num_double = num_double + dim * dim * 3 / n_procs

    required_memory_parallel_generalized = 8 * num_double
  end function required_memory_parallel_generalized


  integer function required_memory(arg)
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
    case ('eigenexa')
      stop 'Eigen Exa is not supported yet'
    case default
      print *, 'Error(command_argument), solver type: ', trim(arg%solver_type)
      stop
    end select
  end function required_memory


  subroutine read_command_argument(arg)
    implicit none

    type(argument), intent(out) :: arg
    integer :: argi = 1
    character(len=256) :: arg_str

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
        case ('v')
          arg%verbose_level = 1
        case default
          call print_help()
          stop 'unknown option'
        end select
      else if (len_trim(arg%matrix_A_filename) == 0) then
        arg%matrix_A_filename = trim(arg_str)
      else
        arg%matrix_B_filename = trim(arg_str)
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

    if (validate_argument(arg) /= 0) then
      stop 'Illegal argument'
    end if
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

    print *, '----- Settings -----'

    if (arg%is_generalized_problem) then
      print *, 'problem type: generalized'
    else
      print *, 'problem type: standard'
    end if

    print *, 'matrix A file: ', trim(arg%matrix_A_filename)
    call print_matrix_info('A', arg%matrix_A_info)

    if (arg%is_generalized_problem) then
      print *, 'matrix B file: ', trim(arg%matrix_B_filename)
      call print_matrix_info('B', arg%matrix_B_info)
    end if

    print *, 'solver: ', trim(arg%solver_type)
    print *, 'output file: ', trim(arg%output_filename)
    print *, 'computed eigenpairs: ', arg%n_vec
    print *, 'verified eigenpairs: ', arg%n_check_vec
  end subroutine print_command_argument
end module command_argument
