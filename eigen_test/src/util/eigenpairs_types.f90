module eigenpairs_types
  type eigenpairs_local ! type number = 1
    double precision, pointer :: values(:)
    double precision, pointer :: vectors(:, :)
    integer :: vector_to_value_index_start, vector_to_value_index_end
  end type eigenpairs_local

  type eigenpairs_blacs ! type number = 2
    double precision, pointer :: values(:)
    integer :: desc(9)
    double precision, pointer :: Vectors(:, :)
    integer :: vector_to_value_index_start, vector_to_value_index_end
  end type eigenpairs_blacs

  type eigenpairs_types_union
    integer :: type_number
    type(eigenpairs_local) :: local
    type(eigenpairs_blacs) :: blacs
  end type eigenpairs_types_union
contains

!  subroutine initialize_eigenpairs_local(dim, n_vec, eigenpairs)
!    integer, intent(in) :: dim, n_vec
!    type(eigenpairs_types_union), intent(out) :: eigenpairs
!
!    eigenpairs%type_number = 1
!
!    if (allocated(eigenpairs%local%values)) then
!      deallocate(eigenpairs%local%values)
!    end if
!
!    if (allocated(eigenpairs%local%vectors)) then
!      deallocate(eigenpairs%local%vectors)
!    end if
!
!    allocate(eigenpairs%local%values(dim))
!    allocate(eigenpairs%local%vectors(dim, n_vec))
!
!    eigenpairs%local%values(:) = 0.0d0
!    eigenpairs%local%vectors(:, :) = 0.0d0
!  end subroutine initialize_eigenpairs_local

end module eigenpairs_types
