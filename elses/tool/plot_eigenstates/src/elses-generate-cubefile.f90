!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
program elses_generate_cubefile

  implicit none

  integer :: i_verbose
  integer :: level_index
  integer :: natom,n_tot_base, n_eigen_vectors
  integer, parameter :: nmo=1
  integer :: level_index_ini, level_index_fin
  logical :: sign_inversion

  real(8) :: unitvectorA(3), unitvectorB(3), unitvectorC(3)
! integer :: mesh_point=80
  integer :: mesh_point_x
  integer :: mesh_point_y
  integer :: mesh_point_z
!
  integer :: mesh_points(3)
  real(8) :: origin_of_grid_region(3)
  real(8) :: size_of_grid_region(3)
  character(len=50) :: unit_of_origin    !  'au', 'internal', or 'nodata'
  character(len=50) :: unit_of_size      !  'au', 'internal'  or 'nodata'
!
  integer :: i_pbc, i_resize
  integer :: i_pbc_x, i_pbc_y, i_pbc_z

  real(8), allocatable :: atmp(:,:)

  real(8) :: r_cut   ! cutoff for plot in au
!
  TYPE Tatom_info
     integer :: atomic_number
     integer :: num_vaL_orb
     integer :: principal(9)
     integer :: angular(9)
     real(8) :: zeta(9)
     real(8) :: zeta2(9)
     real(8) :: c1(9)
     real(8) :: c2(9)
     real(8) :: position(3)
  END TYPE Tatom_info

  type(Tatom_info), allocatable :: atom_info(:)

  character(len=64) :: filename_input
  character(len=64) :: filename_input_wrk
  character(len=50) :: filename_output
  character(len=64) :: cube_filename_header
  character(len=64) :: cube_filename_header_wrk
  character(len=50) :: chara_tmp
  character(len=1024) :: file_format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  filename_input=""
  i_verbose=0
  i_pbc=0    ! This will be replaced in 'read_eigenstates_file'
  i_resize=0
! file_format='current'
!
! file_format='upto_v0.03.09'   ! Activate this line, if needed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  write(*,*)'elses-generate-cubefile'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  filename_input_wrk=""
  cube_filename_header_wrk=""
  write(*,*)'Read the options'
  call read_options(level_index_ini, level_index_fin, sign_inversion, r_cut, filename_input_wrk, cube_filename_header_wrk)
!     -----> The values of  (level_index_ini) and (level_index_fin) are -1 (dummy value),
!                  if not specified
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (trim(adjustl(filename_input_wrk)) == "") then
    filename_input='output_wavefunction.txt'
    write(*,*)'LCAO coef. file name (default) :', trim(adjustl(filename_input))
  else
    filename_input=trim(adjustl(filename_input_wrk))
    write(*,*)'LCAO coef. file name (custom ) :', trim(adjustl(filename_input))
  endif
!
  if (trim(adjustl(cube_filename_header_wrk)) == "") then
    cube_filename_header='eigen_state_'
    write(*,*)'cube file name header (default) :', trim(adjustl(cube_filename_header))
  else
    cube_filename_header=trim(adjustl(cube_filename_header_wrk))
    write(*,*)'cube file name header (custom)  :', trim(adjustl(cube_filename_header))
  endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  call get_file_format_info(file_format)
!
  write(*,*)'File format =', trim(file_format)
!
  write(*,*)'Read the mesh size from the file (optional)'
  call read_mesh_grid(mesh_points, origin_of_grid_region, unit_of_origin, size_of_grid_region, unit_of_size)
  mesh_point_x=mesh_points(1)
  mesh_point_y=mesh_points(2)
  mesh_point_z=mesh_points(3)
!
! write(*,*)'Read the options'
! call read_options(level_index_ini, level_index_fin, sign_inversion, r_cut)
!     -----> The values of  (level_index_ini) and (level_index_fin) are -1 (dummy value),
!                  if not specified
!
  write(*,*)'Read eigenstate file'
  call read_eigenstates_file
!


!
! i_pbc=1                         ! Activate this line, if needed
!
! write(*,*)'file_format =',trim(file_format)
  write(*,*)'i_pbc =', i_pbc
! if (r_cut > 0.0d0) then
!   write(*,*)'OPTION: cutoff radius for plot [au] =', r_cut
! endif
!
  if (i_verbose >= 1) call plot_detailed_information
!
  if (i_pbc == 1) call reduction_with_periodic_bc
  if ((i_pbc == 0) .and. (i_resize ==1)) call resize_mesh_region
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  do level_index=level_index_ini,level_index_fin
!
    write(chara_tmp, '(i6.6)') level_index
    filename_output = trim(adjustl(cube_filename_header))//trim(chara_tmp)//".cube"
!
    call output_cubefile(sign_inversion, r_cut)
!
  enddo
!
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read mesh grid from the file
!
  subroutine read_mesh_grid(mesh_points, origin_of_grid_region, unit_of_origin, size_of_grid_region, unit_of_size)
    implicit none
    integer, intent(out) :: mesh_points(3)
    real(8), intent(out) :: origin_of_grid_region(3)
    real(8), intent(out) :: size_of_grid_region(3)
    character(len=50), intent(out) :: unit_of_origin !  'au', 'internal', or 'nodata'
    character(len=50), intent(out) :: unit_of_size   !  'au', 'internal'  or 'nodata'
    integer, parameter :: iunit=11
    integer :: ierr
    logical :: file_exist
    character(len=32) :: file_name_for_mesh_grid
!
    real(8)           :: data_wrk(3)
!   character(len=50) :: chara_wrk1
    character(len=50) :: chara_wrk2
    integer           :: itry, itry_max
!
    file_name_for_mesh_grid='input_mesh_grid.txt'
!
    unit_of_origin = 'dummy'
    unit_of_size   = 'dummy'
    origin_of_grid_region(:) = -1000.0d0 ! dummy value
    size_of_grid_region(:)   = -1000.0d0 ! dummy value
!
    write(*,*)'file_ name=',trim(file_name_for_mesh_grid)
!
    inquire (file=trim(file_name_for_mesh_grid), exist=file_exist)
    write(*,*)'file_ exist=',file_exist
    if (file_exist .eqv. .true.) then
      open(iunit,file=file_name_for_mesh_grid, status='old')
      read(iunit,*) mesh_points(1:3)
    else
      mesh_points(1:3)=80
    endif
!
    write(*,*)'mesh grid in the x direction=', mesh_points(1)
    write(*,*)'mesh grid in the y direction=', mesh_points(2)
    write(*,*)'mesh grid in the z direction=', mesh_points(3)
!
    if (.not. file_exist) then
      unit_of_origin = 'internal'
      unit_of_size   = 'internal'
      origin_of_grid_region(:) = -0.5d0
      size_of_grid_region(:)   = 1.0d0
      return
    endif
!
    itry_max=2
    do itry=1,itry_max
!     write(*,*)'itry=', itry
      read(iunit,*, iostat=ierr) data_wrk(1:3), chara_wrk2
      if (ierr /= 0) then
        stop 'ERROR in read_mesh_grid'
      endif
      if (itry == 1) then
        if (trim(chara_wrk2)=='au') unit_of_origin='au'
        if (trim(chara_wrk2)=='internal') unit_of_origin='internal'
        if (trim(unit_of_origin) /= 'dummy') then
          origin_of_grid_region(1:3)=data_wrk(1:3)
!         write(*,'(a,3f20.10,a,a)')'grid origin = ', data_wrk(1:3), '    unit=', trim(chara_wrk2)
        endif
      endif
      if (itry == 2) then
        if (trim(chara_wrk2)=='au') unit_of_size='au'
        if (trim(chara_wrk2)=='internal') unit_of_size='internal'
        if (trim(unit_of_size) /= 'dummy') then
          size_of_grid_region(1:3)=data_wrk(1:3)
!         write(*,'(a,3f20.10,a,a)')'grid sizes  = ', data_wrk(1:3), '    unit=', trim(chara_wrk2)
        endif
      endif
    enddo
!
    if (file_exist .eqv. .true.) close(iunit)
!
  end subroutine read_mesh_grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read the options for the level index and the sign value
!
  subroutine read_options(level_index_ini,level_index_fin,sign_inversion,r_cut_wrk, filename_wrk, cube_header_wrk)
!
    implicit none
    integer, intent(out) :: level_index_ini
    integer, intent(out) :: level_index_fin
    logical, intent(out) :: sign_inversion
    character(len=*), intent(inout) :: filename_wrk
    character(len=*), intent(inout) :: cube_header_wrk
    integer :: iargc
    integer :: i,count, n
    integer :: k, ierr
    logical :: acceptable_option
    real(8) :: r_cut_wrk
!
    character(len=80) :: argc
    count=iargc()
!
!   call get_number_of_bases(n)  ! the number of bases
!
    level_index_ini=-1     ! dummy value
    level_index_fin=-1     ! dummy value
    sign_inversion = .false.
    r_cut_wrk = -1.0d0 ! dummy value
!
!   write(*,*)'INFO:the number of the options  : ', count
!
    do i=1, count
      call getarg(i,argc)
!     write(*,*)'i, argc=',i, trim(argc)
      read(argc,*,iostat=ierr) k
      if (ierr /= 0) then
        acceptable_option = .false.
!       write(*,*)'INFO:an option found : ', trim(argc)
        if (argc(1:3)  == '-si') then
          acceptable_option = .true.
          sign_inversion = .true.
          cycle
        endif
        if (argc(1:15) == '-sign_inversion') then
          acceptable_option = .true.
          sign_inversion = .true.
          cycle
        endif
        if (argc(1:11) == '-cutoff_au=') then
          read(argc(12:),*,iostat=ierr) r_cut_wrk
          if (ierr /= 0) then
            write(*,*)'ERROR in option:', trim(argc)
            stop
          endif
!         write(*,*)'OPTION: cutoff [au] =', r_cut_wrk
          cycle
        endif
        if (argc(1:16) == '-lcao_coef_file=') then
          acceptable_option = .true.
          read(argc(17:),*,iostat=ierr) filename_wrk
          if (ierr /= 0) then
            write(*,*)'ERROR in option:', trim(argc)
            stop
          endif
          if (trim(adjustl(filename_wrk)) == '') then
            write(*,*)'ERROR in option :', trim(argc)
            stop
          endif
          filename_wrk=trim(adjustl(filename_wrk))
          cycle
        endif
        if (argc(1:22) == '-cube_filename_header=') then
          acceptable_option = .true.
          read(argc(23:),*,iostat=ierr) cube_header_wrk
          if (ierr /= 0) then
            write(*,*)'ERROR in option:', trim(argc)
            stop
          endif
          if (trim(adjustl(cube_header_wrk)) == '') then
            write(*,*)'ERROR in option :', trim(argc)
            stop
          endif
          cube_header_wrk=trim(adjustl(cube_header_wrk))
          cycle
        endif
        if (.not. acceptable_option) then
          write(*,*)'ERROR(generate cube) in option : ', trim(argc)
          stop
        endif
      else
        if (level_index_ini == -1) then
!         write(*,*)'level_index_ini = ', k
          level_index_ini = k
          cycle
        else
          if (level_index_fin == -1) then
!           write(*,*)'level_index_fin = ', k
            level_index_fin = k
            cycle
          endif
        endif
      endif
    enddo
!
    write(*,*)'INFO:sign inversion = ', sign_inversion
!
    if (level_index_ini /= -1) then
      write(*,*)'INFO:the lowest  level index = ', level_index_ini
    endif
!
    if (level_index_fin /= -1) then
      write(*,*)'INFO:the highest level index = ', level_index_fin
    endif
!
    return
!

  end subroutine read_options
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read the file format information, if available
!
  subroutine get_file_format_info(file_format_result)
    implicit none
    character(len=1024), intent(out) :: file_format_result
    character(len=1024)              :: str_wrk
    integer :: ierr
!
    open(10,file=filename_input,status='old')
      read(10,'(a)',iostat=ierr) str_wrk
      if (ierr /=0) then
        write(*,*)'READ ERROR(get_file_format_info)'
        stop
      endif
    close(10)
!
    if (index(str_wrk, 'file_format') == 0) then
      file_format_result='not_specified'
      return
    endif
!
    str_wrk=trim(adjustl(str_wrk))
    str_wrk=str_wrk(16:)
    file_format_result=trim(adjustl(str_wrk))
!
    ierr=1
    if (trim(adjustl(file_format_result)) == 'v0.04.05') ierr=0
!
    write(*,*)'file format=',trim(adjustl(file_format_result))
!
    if (ierr /=0) then
      write(*,*)'ERROR(get_file_format_info):file_format=',trim(adjustl(file_format_result))
      stop
    endif
!
  end subroutine get_file_format_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read the number of the bases from eigenstates.dat
!
  subroutine get_number_of_bases(number_of_bases)

    implicit none
    integer, intent(out) :: number_of_bases
    integer :: number_of_atoms
    integer :: ierr
    character(len=1024)              :: str_wrk
!
    open(10,file=filename_input,status='old')
!
    if (trim(file_format) /= 'not_specified') then
      read(10,'(a)',iostat=ierr) str_wrk
!       ---> skip the first line of the file
!     write(*,*)' first line=',trim(str_wrk)
    endif
!
!   read(10,'(a)',iostat=ierr) str_wrk
!   write(*,*)' secondt line=',trim(str_wrk)
!   stop 'stop manually'
!
    read(10,*,iostat=ierr) number_of_atoms, number_of_bases
    if (ierr /= 0) then
      write(*,*)'ERROR:get_number_of_bases'
      stop
    endif
!
!   write(*,*) ' number of bases =', number_of_bases
    close(10)
!
  end subroutine get_number_of_bases
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read system information and eigenstates
!from eigenstates.dat
!
  subroutine read_eigenstates_file

    implicit none

    integer :: i,j,index_j, index_i
    integer :: ierr, ierr2
    character(len=16) :: chara
    character(len=1024) :: str_wrk
!
    ierr2=0
!
    open(10,file=filename_input,status='old')
!
    if (trim(file_format) /= 'not_specified') then
      read(10,'(a)',iostat=ierr) str_wrk
!       ---> skip the first line of the file
!     write(*,*)' first line=',trim(str_wrk)
    endif
!
!
    read(10,*,iostat=ierr) natom, n_tot_base
    if (ierr /= 0) then
      write(*,*)'ERROR(read_eigenstates_file):wrong format'
      stop
    endif
!
    write(*,*) ' number of atoms         =', natom
    write(*,*) ' number of bases         =', n_tot_base
!   write(*,*) ' number of eigen vectors =', n_eigen_vectors
!
    write(*,*)' file_format = ', trim(file_format)
!
    i_pbc_x=0
    i_pbc_y=0
    i_pbc_z=0
!
    if (trim(file_format) /= 'upto_v0.03.09') then
      read(10,*,iostat=ierr) unitvectorA(1:3), i_pbc_x
      write(*,*)' i_pbc_x appears in the file'
    else
      read(10,*,iostat=ierr) unitvectorA(1:3)
      write(*,*)' i_pbc_x does not appear in the file'
      i_pbc_x=0
    endif
!
    if (unitvectorA(1)       < 0.0d0)   ierr2=1
    if (dabs(unitvectorA(2)) > 1.0d-10) ierr2=1
    if (dabs(unitvectorA(3)) > 1.0d-10) ierr2=1
    if (i_pbc_x*i_pbc_x /= i_pbc_x)     ierr2=1   ! i_pbc should be zero or one
!
    if ((ierr /=0) .or. (ierr2 /=0)) then
      write(*,*)' ERROR in IO (unitvectorA)'
      write(*,*)'   The file format may be old (upto_v0.03.09)'
      stop
    endif
!
    if (trim(file_format) /= 'upto_v0.03.09') then
      read(10,*,iostat=ierr) unitvectorB(1:3), i_pbc_y
      write(*,*)' i_pbc_y appears in the file'
    else
      read(10,*,iostat=ierr) unitvectorB(1:3)
      write(*,*)' i_pbc_y does not appear in the file'
      i_pbc_y=0
    endif
!
    if (unitvectorB(2)       < 0.0d0)   ierr2=1
    if (dabs(unitvectorB(1)) > 1.0d-10) ierr2=1
    if (dabs(unitvectorB(3)) > 1.0d-10) ierr2=1
    if (i_pbc_y*i_pbc_y /= i_pbc_y)     ierr2=1   ! i_pbc should be zero or one
!
    if ((ierr /=0) .or. (ierr2 /=0)) then
      write(*,*)' ERROR in IO (unitvectorB)'
      write(*,*)'   The file format may be old (upto_v0.03.09)'
      stop
    endif
!
    if (trim(file_format) /= 'upto_v0.03.09') then
      read(10,*,iostat=ierr) unitvectorC(1:3), i_pbc_z
      write(*,*)' i_pbc_z appears in the file'
    else
      read(10,*,iostat=ierr) unitvectorC(1:3)
      write(*,*)' i_pbc_z does not appear in the file'
      i_pbc_z=0
    endif
!
    if (unitvectorC(3)       < 0.0d0)   ierr2=1
    if (dabs(unitvectorC(1)) > 1.0d-10) ierr2=1
    if (dabs(unitvectorC(2)) > 1.0d-10) ierr2=1
    if (i_pbc_z*i_pbc_z /= i_pbc_z)     ierr2=1   ! i_pbc should be zero or one
!
    if ((ierr /=0) .or. (ierr2 /=0)) then
      write(*,*)' ERROR in IO (unitvectorC)'
      write(*,*)'   The file format may be old (upto_v0.03.09)'
      stop
    endif
!
    if (trim(file_format) /= 'upto_v0.03.09') then
      i_pbc=i_pbc_x*i_pbc_y*i_pbc_z
      write(*,*) 'FileRead:i_pbc_x, y, z=',i_pbc_x, i_pbc_y, i_pbc_z
      write(*,*) 'FileRead:i_pbc=',i_pbc
    endif
!
    if (i_pbc_x /= i_pbc_y) then
      write(*,*)'ERROR;Now not supported:Periodic in 1 or 2 direction(s)'
      stop
    endif
!
    if (i_pbc_y /= i_pbc_z) then
      write(*,*)'ERROR;Now not supported:Periodic in 1 or 2 direction(s)'
      stop
    endif
!
    if (i_pbc_z /= i_pbc_x) then
      write(*,*)'ERROR;Now not supported:Periodic in 1 or 2 direction(s)'
      stop
    endif
!

    write(*,*)'FileRead : atom information'
!
    allocate(atom_info(natom))

    do i=1, natom
       atom_info(i)%atomic_number=0
       atom_info(i)%num_val_orb=0
       atom_info(i)%principal=0
       atom_info(i)%angular=0
       atom_info(i)%zeta=0.0d0
       atom_info(i)%zeta2=0.0d0
       atom_info(i)%c1=0.0d0
       atom_info(i)%c2=0.0d0
       atom_info(i)%position=0.0d0
    end do

    do i=1, natom
!
       read(10,*) atom_info(i)%atomic_number
       if ((atom_info(i)%atomic_number > 100) .or. (atom_info(i)%atomic_number < 1)) then
          write(*,*)'ERROR:atomic number',i,atom_info(i)%atomic_number
          stop
       endif
!
       read(10,*) atom_info(i)%num_val_orb
       if ((atom_info(i)%num_val_orb > 9) .or. (atom_info(i)%num_val_orb < 1)) then
          write(*,*)'ERROR:number of orbitals',i,atom_info(i)%num_val_orb
          stop
       endif
!
       read(10,*) atom_info(i)%principal(1:atom_info(i)%num_val_orb)
       read(10,*) atom_info(i)%angular(1:atom_info(i)%num_val_orb)
       read(10,*) atom_info(i)%zeta(1:atom_info(i)%num_val_orb)
       read(10,*) atom_info(i)%zeta2(1:atom_info(i)%num_val_orb)
       read(10,*) atom_info(i)%c1(1:atom_info(i)%num_val_orb)
       read(10,*) atom_info(i)%c2(1:atom_info(i)%num_val_orb)
       read(10,*) atom_info(i)%position(1:3)
!
       if (i_verbose >=1) then
         write(*,'(I2)') atom_info(i)%atomic_number
         write(*,'(I8)') atom_info(i)%num_val_orb
         write(*,'(10I6)') atom_info(i)%principal(1:9)
         write(*,'(10I6)') atom_info(i)%angular(1:9)
         write(*,'(4ES10.3)') atom_info(i)%zeta(1:9)
         write(*,'(3F23.15)') atom_info(i)%position(1:3)
       endif

    end do
!
    write(*,*)'FileRead : coefficient information'
!
    n_eigen_vectors=n_tot_base
    if (trim(file_format) /= 'not_specified') then
      read(10,'(a)',iostat=ierr) str_wrk
      if (ierr /= 0) then
        write(*,*) 'ERROR:wrong separator'
        stop
      endif
      if (index(trim(str_wrk),'LCAO') == 0) then
        write(*,*) 'ERROR:wrong separator:',trim(str_wrk)
        stop
      endif

      write(*,*)'separator=',trim(str_wrk)
      read(10,*) n_eigen_vectors
    endif
    write(*,*) ' number of eigen vectors =', n_eigen_vectors
!
    allocate(atmp(n_tot_base,n_eigen_vectors))

    if (trim(file_format) == 'upto_v0.03.09') then
      do i=1,n_eigen_vectors
        do j=1,n_tot_base
          read(10,*) index_j, atmp(j,i)
        enddo
      enddo
    else
      do i=1,n_eigen_vectors
        do j=1,n_tot_base
          read(10,*) index_j, atmp(j,i), chara, index_i
!         write(*,'(a,2i10,f20.10)')' wfn=',i,j,atmp(j,i)
          if (i /= index_i) then
            write(*,*) 'ERROR:Unmatched index_i',index_i,i
            stop
          endif
          if (j /= index_j) then
            write(*,*) 'ERROR:Unmatched index_j',index_j,j
            stop
          endif
        end do
      end do
    endif
!
    close(10)
!
    if (level_index_ini == -1) then
      level_index_ini = 1
      level_index_fin = n_eigen_vectors
    else
      if (level_index_fin == -1) then
        level_index_fin = level_index_ini
      endif
    endif
!
    ierr=0
    if (level_index_ini < 1) ierr=1
    if (level_index_fin > n_eigen_vectors) ierr=1
    if (level_index_ini > level_index_fin) ierr=1
!
    if (ierr /= 0) then
      write(*,*)'ERROR!(generate cube) :the lowest  level index = ', level_index_ini
      write(*,*)'ERROR!(generate cube) :the highest level index = ', level_index_fin
      stop
    endif
!
    write(*,*)'INFO:the lowest  level index for output = ', level_index_ini
    write(*,*)'INFO:the highest level index for output = ', level_index_fin
!
  end subroutine read_eigenstates_file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculation of amptitudes of wave function at
!every mesh points, and output a cube file.
!
  subroutine output_cubefile(sign_inversion, r_cut_wrk)

    implicit none
    logical, intent(in) :: sign_inversion
    real(8), intent(in) :: r_cut_wrk ! cutoff for plot in au
!
    integer :: i,j,k,atom_index,orbital_index,wf_index,mesh_index
    real(8) :: x,y,z,l,m,n,r
    real(8) :: sto
    real(8), parameter :: pi=3.14159265358979323
    real(8) :: zeta,zeta2,c1,c2,f
    real(8) :: wf !!! amplitude of wave function
    real(8) :: dx, dy, dz
!   real(8) :: dd
    real(8) :: origin_of_grid_region_au(3)   ! origin of the grid region in au
    real(8) :: size_of_grid_region_au(3)     ! size   of the grid region in au
    real(8) :: size_of_cell_au(3)            ! size   of the cell
    real(8), parameter :: zero_value=0.0d0
    real(8) :: sign_value                    ! = +1.0d0 or -1.0d0
    logical :: cutoff_is_imposed
!
    if (sign_inversion) then
      sign_value = -1.0d0
    else
      sign_value =  1.0d0
    endif
!
    if (r_cut_wrk > 1.0d-10) then
      write(*,*)'OPTION:The cutoff radius for plot [au] =', r_cut_wrk
      cutoff_is_imposed = .true.
    else
      cutoff_is_imposed = .false.
    endif
!
    size_of_cell_au(1) = unitvectorA(1)   ! Re-definition
    size_of_cell_au(2) = unitvectorB(2)   ! Re-definition
    size_of_cell_au(3) = unitvectorC(3)   ! Re-definition
!
    origin_of_grid_region_au(1) = -unitvectorA(1)/2.0  ! Default setting
    origin_of_grid_region_au(2) = -unitvectorB(2)/2.0  ! Default setting
    origin_of_grid_region_au(3) = -unitvectorC(3)/2.0  ! Default setting
!
    size_of_grid_region_au(1)   =  unitvectorA(1)      ! Default setting
    size_of_grid_region_au(2)   =  unitvectorB(2)      ! Default setting
    size_of_grid_region_au(3)   =  unitvectorC(3)      ! Default setting
!
    if (trim(unit_of_origin) == 'internal') then
      origin_of_grid_region_au(1:3)=origin_of_grid_region(1:3)*size_of_cell_au(1:3)
    endif
!
    if (trim(unit_of_origin) == 'au') then
      origin_of_grid_region_au(1:3)=origin_of_grid_region(1:3)
    endif
!
    if (trim(unit_of_size) == 'internal') then
      size_of_grid_region_au(1:3)=size_of_grid_region(1:3)*size_of_cell_au(1:3)
    endif
!
    if (trim(unit_of_origin) == 'au') then
      size_of_grid_region_au(1:3)=size_of_grid_region(1:3)
    endif
!
    write(*,*)'Generate a cube file:level_index=',level_index
!
    write(*,'(a,3f20.10)')'grid origin [au] = ', origin_of_grid_region_au(1:3)
    write(*,'(a,3f20.10)')'grid size   [au] = ', size_of_grid_region_au(1:3)
!
    open(10,file=filename_output)

    write(10,'("ELSES MO")')
    write(10,'(I10,"-th MO")') level_index

    write(10,'(I6,3F16.8)') -natom, origin_of_grid_region_au(1:3)
    write(10,'(I6,3F16.8)') mesh_point_x,size_of_grid_region_au(1)/mesh_point_x, zero_value, zero_value
    write(10,'(I6,3F16.8)') mesh_point_y,zero_value, size_of_grid_region_au(2)/mesh_point_y, zero_value
    write(10,'(I6,3F16.8)') mesh_point_z, zero_value, zero_value, size_of_grid_region_au(3)/mesh_point_z

    do i=1, natom
       write(10,'(I2,1X,I2,1X,3F16.8)') atom_info(i)%atomic_number, &
            atom_info(i)%atomic_number,atom_info(i)%position(1:3)
    end do

    write(10,'(2I10)') nmo, level_index

    do i=0, mesh_point_x-1
       mesh_index=1
       do j=0,mesh_point_y-1
          mesh_index=1
          do k=0, mesh_point_z-1

                x = origin_of_grid_region_au(1) + size_of_grid_region_au(1)/mesh_point_x*i
                y = origin_of_grid_region_au(2) + size_of_grid_region_au(2)/mesh_point_y*j
                z = origin_of_grid_region_au(3) + size_of_grid_region_au(3)/mesh_point_z*k
!
!               x=-unitvectorA(1)/2.0+unitvectorA(1)/mesh_point_x*i
!               y=-unitvectorB(2)/2.0+unitvectorB(2)/mesh_point_y*j
!               z=-unitvectorC(3)/2.0+unitvectorC(3)/mesh_point_z*k

!               write(*,'("(x,y,z)",3F20.10)') x,y,z

                wf=0.0d0
                wf_index=1

             do atom_index=1,natom

 !               write(*,'("atom_index=",I)') atom_index
!
                 dx=(x-atom_info(atom_index)%position(1))/unitvectorA(1)
                 dy=(y-atom_info(atom_index)%position(2))/unitvectorB(2)
                 dz=(z-atom_info(atom_index)%position(3))/unitvectorC(3)
                 if (i_pbc == 1) then
                   dx = modulo(dx + 0.5d0, 1.0d0) - 0.5d0
                   dy = modulo(dy + 0.5d0, 1.0d0) - 0.5d0
                   dz = modulo(dz + 0.5d0, 1.0d0) - 0.5d0
                 endif
                 dx=dx*unitvectorA(1)
                 dy=dy*unitvectorB(2)
                 dz=dz*unitvectorC(3)
                 r=dsqrt(dx**2+dy**2+dz**2)
                 l=dx
                 m=dy
                 n=dz
!
                 if (cutoff_is_imposed) then
                   if (r > r_cut_wrk) then
                     wf_index = wf_index+atom_info(atom_index)%num_val_orb
                     cycle
                   endif
                 endif
!
!                if(r < 1e-10) then
!                  r=0.0d0
!                  l=0.0d0
!                  m=0.0d0
!                  n=0.0d0
!                else
!                  l=dx/r
!                  m=dy/r
!                  n=dz/r
!                end if
!
!                write(*,'("r=",F)') r
!                write(*,'("(l,m,n)=",3F)') l,m,n

!!!! Until now, the visualization of the 4d orbital and higher orbital
!!!! was not implemented

                do orbital_index=1, atom_info(atom_index)%num_val_orb

                   sto=0.0d0

                   zeta=atom_info(atom_index)%zeta(orbital_index)
                   zeta2=atom_info(atom_index)%zeta2(orbital_index)
                   c1=atom_info(atom_index)%c1(orbital_index)
                   c2=atom_info(atom_index)%c2(orbital_index)

                   if(c2 == 0d0) then
                      c1=1.0d0
                      f=1.0d0
                   else
                      f=sqrt(c1*c1+c2*c2+2.0d0*c1*c2* &
                           (4*zeta*zeta2/(zeta+zeta2)**2.0d0)**(7.0d0/2.0d0))
                   end if

                   if(atom_info(atom_index)%principal(orbital_index)==1) then
                       sto=1.0d0/sqrt(pi)/f* &
                            (c1*(zeta**(3.0d0/2.0d0))*exp(-zeta*r) &
                            +c2*(zeta2**(3.0d0/2.0d0))*exp(-zeta2*r))
!                       write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                   else if(atom_info(atom_index)%principal(orbital_index)==2) then
                      if(orbital_index==1) then
                         sto=1.0/sqrt(3.0d0*pi)/f*r* &
                              (c1*(zeta**(5.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(5.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2,"  sto=",F)') orbital_index,sto
                      else if(orbital_index==2) then
                         sto=1.0/sqrt(pi)/f*l* &
                              (c1*(zeta**(5.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(5.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2,"  sto=",F)') orbital_index,sto
                      else if(orbital_index==3) then
                         sto=1.0/sqrt(pi)/f*m* &
                              (c1*(zeta**(5.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(5.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2,"  sto=",F)') orbital_index,sto
                      else if(orbital_index==4) then
                         sto=1.0/sqrt(pi)/f*n* &
                              (c1*(zeta**(5.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(5.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2,"  sto=",F)') orbital_index,sto
                      end if
                   else if(atom_info(atom_index)%principal(orbital_index)==3) then
                      if(orbital_index==1) then
                         sto=2.0d0/3.0d0/sqrt(10.0d0*pi)/f*(r**2)* &
                              (c1*(zeta**(7.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(7.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==2) then
                         sto=1.0d0/sqrt(15.0d0*pi)/f*r*l* &
                              (c1*(zeta**(7.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(7.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==3) then
                         sto=1.0d0/sqrt(15.0d0*pi)/f*r*m* &
                              (c1*(zeta**(7.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(7.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==4) then
                         sto=1.0d0/sqrt(15.0d0*pi)/f*r*n* &
                              (c1*(zeta**(7.0d0/2.0d0))*exp(-zeta*r) &
                              +c2*(zeta2**(7.0d0/2.0d0))*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==5) then
                         sto=1.0d0/sqrt(6.0d0*pi)*2.0d0/f*l*m* &
                             (c1*zeta**(7.0d0/2.0d0)*exp(-zeta*r) &
                             +c2*zeta2**(7.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==6) then
                         sto=1.0d0/sqrt(6.0d0*pi)*2.0d0*m*n* &
                             (c1*zeta**(7.0d0/2.0d0)*exp(-zeta*r) &
                             +c2*zeta2**(7.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==7) then
                         sto=1.0d0/sqrt(6.0d0*pi)*2.0d0*n*l* &
                             (c1*zeta**(7.0d0/2.0d0)*exp(-zeta*r) &
                             +c2*zeta2**(7.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==8) then
                         sto=1.0d0/sqrt(6.0d0*pi)*(l*l-m*m)* &
                             (c1*zeta**(7.0d0/2.0d0)*exp(-zeta*r) &
                             +c2*zeta2**(7.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==9) then
                         sto=1.0d0/sqrt(6.0d0*pi)*(3*n*n-r*r)/sqrt(3.0d0)* &
                             (c1*zeta**(7.0d0/2.0d0)*exp(-zeta*r) &
                             +c2*zeta2**(7.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("orbital_index=",I2," sto=",F)') orbital_index,sto
                      end if
                   else if(atom_info(atom_index)%principal(orbital_index)==4) then
                      if(orbital_index==1) then
                         sto=1.0d0/sqrt(105.0d0)/f*r**3.0d0* &
                              (c1*zeta**(9.0d0/2.0d0)*exp(-zeta*r) &
                              +c2*zeta2**(9.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("n=4 orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==2) then
                         sto=1.0d0/sqrt(35.0d0*pi)/f*l*r**2.0d0* &
                              (c1*zeta**(9.0d0/2.0d0)*exp(-zeta*r) &
                              +c2*zeta2**(9.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("n=4 orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==3) then
                         sto=1.0d0/sqrt(35.0d0*pi)/f*m*r**2.0d0* &
                              (c1*zeta**(9.0d0/2.0d0)*exp(-zeta*r) &
                              +c2*zeta2**(9.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("n=4 orbital_index=",I2," sto=",F)') orbital_index,sto
                      else if(orbital_index==4) then
                         sto=1.0d0/sqrt(35.0d0*pi)/f*n*r**2.0d0* &
                              (c1*zeta**(9.0d0/2.0d0)*exp(-zeta*r) &
                              +c2*zeta2**(9.0d0/2.0d0)*exp(-zeta2*r))
!                         write(*,'("n=4 orbital_index=",I2," sto=",F)') orbital_index,sto
                      else
                         write(*,'("!!! ERROR !!! Unexpected orbital is included")')
                         stop
                      end if
                   else
                      write(*,'("!!! ERROR !!! Unexpected orbital is included")')
                      stop
                   end if

!                   write(*,'("wf_index=",I10)') wf_index
!                   write(*,'("atmp(",I10,",",I10,")=",F)') wf_index,level_index, &
!                   atmp(wf_index,level_index)

                   wf=wf+sto*atmp(wf_index,level_index)

!                   write(*,'(" wf=",F)') wf

                   wf_index=wf_index+1

                end do

!                write(*,'(" wf=",F\)') wf
!                write(*,'("")')

             end do
!
!            if (dabs(wf) <= 1.0d-80) wf=0.0d0
             if (dabs(wf) >= 1.0d80) then
                write(*,*)'ERROR:Too large amplitude'
                stop
             endif
!
!             write(18,'("wf=",E15.8)') wf
!
             wf = wf * sign_value
!
             ! Start a new line before writing 6n+1 th elements.
             if(mesh_index > 1 .and. mod(mesh_index-1,6)==0) then
                write(10, '()')
             end if
             write(10,'(1X,E18.8e3)',advance='no') wf

             mesh_index=mesh_index+1
          end do
          write(10,'()')
       end do
    end do

!   deallocate(atom_info)
!   deallocate(atmp)

    close(10)

  end subroutine output_cubefile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Reduce the atomic positions within the periodic boundary condition
!
  subroutine reduction_with_periodic_bc
!
    implicit none
    integer :: i,j
    real(8) :: pos(3)
    real(8) :: cell(3)
!
    cell(1)=unitvectorA(1)
    cell(2)=unitvectorB(2)
    cell(3)=unitvectorC(3)
!
    write(*,*)'@@ reduction_within_periodic_bc'
    write(*,*)'  cell(1)=',cell(1)
    write(*,*)'  cell(2)=',cell(2)
    write(*,*)'  cell(3)=',cell(3)
!
    do i=1, natom
      pos(1:3) = atom_info(i)%position(1:3)
      pos(1:3) = pos(1:3)/cell(1:3)
!     write(*,*)'i,tx0=',i,pos(1)
!     write(*,*)'i,ty0=',i,pos(2)
!     write(*,*)'i,tz0=',i,pos(3)
      do j=1,3
        pos(j)=pos(j)-dble(int(pos(j)))
        if (pos(j) < 0) then
          pos(j)=pos(j)+1.0d0
        endif
      enddo
!     write(*,*)'i,tx1=',i,pos(1)
!     write(*,*)'i,ty1=',i,pos(2)
!     write(*,*)'i,tz1=',i,pos(3)
      pos(1:3) = pos(1:3)*cell(1:3)
!     write(*,*)'i,tx=',i,pos(1), atom_info(i)%position(1)
!     write(*,*)'i,ty=',i,pos(2), atom_info(i)%position(2)
!     write(*,*)'i,tz=',i,pos(3), atom_info(i)%position(3)
      atom_info(i)%position(1:3)=pos(1:3)
    enddo
!   stop
  end subroutine reduction_with_periodic_bc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Resize the mesh region (in non PBC)
!
  subroutine resize_mesh_region
!
    implicit none
    integer :: i,j
    real(8) :: pos_max(3), pos_min(3), pos_cen(3), ddd
    real(8) :: cell(3)
    real(8) :: supercell_scale(3)
!
    cell(1)=unitvectorA(1)
    cell(2)=unitvectorB(2)
    cell(3)=unitvectorC(3)
!
    do j=1,3
      do i=1, natom
        ddd=atom_info(i)%position(j)/cell(j)
        if ( i == 1 ) then
          pos_cen(j)=0.0d0
          pos_max(j)=ddd
          pos_min(j)=ddd
        endif
        pos_cen(j)=pos_cen(j)+ddd/dble(natom)
        if (ddd > pos_max(j)) pos_max(j)=ddd
        if (ddd < pos_min(j)) pos_min(j)=ddd
      enddo
      write(*,*)'j,max,min,cen=', j, pos_max(j), pos_min(j), pos_cen(j)
    enddo
!
    do j=1,3
      ddd=pos_max(j)
      if (ddd > 1.0d0) then
         pos_max(j)=dble(int(ddd)+1)
         write(*,'(a,i10,2f20.10)')'j, Max_region, Max_pos=',j,pos_max(j), ddd
      endif
      ddd=pos_min(j)
      if (ddd < 0.0d0) then
         pos_min(j)=-dble(int(abs(ddd))+1)
         write(*,'(a,i10,2f20.10)')'j, Min_region, Min_pos=',j,pos_min(j), ddd
      endif
      supercell_scale(j)=pos_max(j)-pos_min(j)
      write(*,'(a,i10,f20.10)')'j, supercell_scale =',j, supercell_scale(j)
    enddo
!
!   Set the origin ((x,y,z)=(0,0,0)) to be the center of the visible region
    do j=1,3
      do i=1, natom
        ddd=atom_info(i)%position(j)/cell(j)
        ddd=ddd-(pos_max(j)+pos_min(j))/2.0d0
        atom_info(i)%position(j)=ddd*cell(j)
        write(*,*)'j,i,pos_new=',j,i,atom_info(i)%position(j)
      enddo
    enddo
!
    supercell_scale(1:3)=supercell_scale(1:3)*2.0d0
!
    unitvectorA(1)=cell(1)*supercell_scale(1)
    unitvectorB(2)=cell(2)*supercell_scale(2)
    unitvectorC(3)=cell(3)*supercell_scale(3)
!
!   stop
!
  end subroutine resize_mesh_region
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Plot detailed information (optional)
!
  subroutine plot_detailed_information

    implicit none
    integer, parameter :: iunit1=21
    integer, parameter :: iunit2=22
    integer :: j, ierr
    integer :: lvl_index, atm_index, orb_index , atm_num_index
    integer, parameter :: max_atomic_number=100
    integer :: max_nval
    real(8) :: coef
    real(8), allocatable :: weight(:,:)
    integer, allocatable :: weight_flag(:)
!
    character(len=*), parameter :: filename1="output_wavefunction_checked.txt"
    character(len=*), parameter :: filename2="output_wavefunction_analyzed.txt"
!
    write(*,*)'Check and plot detailed information'
    write(*,*)'n_tot_base=',n_tot_base
!
    max_nval=maxval(atom_info(:)%num_val_orb)
    write(*,*)'  max_nvaL=',max_nval
!
    allocate(weight_flag(max_atomic_number),stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error'
      stop
    endif
!
    allocate(weight(max_nval, max_atomic_number),stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error'
      stop
    endif
!
    open(iunit1,file=filename1, status='unknown')
    open(iunit2,file=filename2, status='unknown')
!
    write(iunit2,*)'Analysis on wavefunction weight (or |coefficient|^2) [experimental function]'
!
    do lvl_index=1,n_tot_base
      weight(:,:)=0.0d0
      weight_flag(:)=0
!
      j=0
      do atm_index=1,natom
!       write(*,*)'atm_index=',atm_index
        atm_num_index=atom_info(atm_index)%atomic_number
        weight_flag(atm_num_index)=1
        do orb_index=1, atom_info(atm_index)%num_val_orb
          j=j+1
          if (j > n_tot_base) then
            write(*,*)'ERROR(plot_detailed_information) '
          endif
          coef=atmp(j,lvl_index)
          write(iunit1,'(4i10,f20.10)') lvl_index, atm_index, atm_num_index, orb_index,coef
          weight(orb_index, atm_num_index)=weight(orb_index, atm_num_index)+coef*coef
        end do
      end do
!
!     write(iunit2,*)'level_index=',lvl_index
      do atm_num_index=1, max_atomic_number
        if (weight_flag(atm_num_index)  == 0) cycle
        if (max_nval == 4) then
          write(iunit2,'(2i10,f25.10,5f10.5)') lvl_index, atm_num_index, &
               sum(weight(:,atm_num_index)), weight(1, atm_num_index), &
               sum(weight(2:4, atm_num_index)), weight(2:4, atm_num_index)
        else
           write(iunit2,'(2i10,f25.10, 3f10.5)') lvl_index, atm_num_index, &
                sum(weight(:,atm_num_index)), weight(1, atm_num_index),&
                sum(weight(2:4, atm_num_index)), sum(weight(5:9, atm_num_index))
        endif
      enddo
    end do
!
    close(iunit1)
    close(iunit2)
!
    write(*,*)'... The data check was carried out without problem'
!
  end subroutine plot_detailed_information
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end program elses_generate_cubefile
