module M_gid_basic_routines
!
  use M_config, only : config !(unchanged) 
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  private
  public :: fix_group_center
  public :: save_w_group_centers
!
  contains
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Fix the group center for the 'gid'-th group
!
  subroutine fix_group_center(gid, mode)
!
    use M_group_id_setting, only : num_groups       !(unchanged)
    use M_group_id_setting, only : group_center_ini !(unchanged)
    use M_group_id_setting, only : set_weight_center_sub !(routine)
    implicit none
    integer,                intent(in)  :: gid
    character(len=*),       intent(in)  :: mode
    real(DOUBLE_PRECISION)              :: group_center_wrk(3)
    real(DOUBLE_PRECISION)              :: diff_center(3)
    real(DOUBLE_PRECISION)              :: shift_vector(3)
    integer                             :: lu
    integer                             :: step_count
!
    lu=config%calc%distributed%log_unit
    step_count=config%system%structure%mdstep
!
    if (lu > 0) write(lu,'(a)')'@@@ fix_group_center'
!
    if ( (gid < 1) .or. (gid > num_groups) ) then
      write(*,*)'ERROR(fix_group_center):gid=',gid
    endif
!
    call set_weight_center_sub(gid, group_center_wrk)
!
    if (lu > 0) write(lu,'(a,i10, a, i10, 3f15.10)') 'step_count=', step_count, & 
&                                                    '  group_center(old) = ',gid, group_center_wrk(1:3)
    if (lu > 0) write(lu,'(a,i10, a, i10, 3f15.10)') 'step_count=', step_count, & 
&                                                    '  group_center(ini) = ',gid, group_center_ini(1:3,gid)
!
    diff_center(1:3)=group_center_wrk(1:3)-group_center_ini(1:3,gid)
!
    shift_vector(:)=0.0d0
!
    select case(trim(mode))
      case('xyz')
        shift_vector(1:3)=-diff_center(1:3)
      case('x')
        shift_vector(1)=-diff_center(1)
      case('y')
        shift_vector(2)=-diff_center(2)
      case('z')
        shift_vector(3)=-diff_center(3)
      case default
        write(*,*)'ERROR(shift_group_center):mode=',trim(mode)
    end select
!
    if (lu > 0) write(lu,'(a,3f15.10)') '  diff_center      = ',  diff_center(1:3)
    if (lu > 0) write(lu,'(a,3f15.10)') ' shift_center      = ',  shift_vector(1:3)
!
    call shift_group_center(gid, shift_vector)
!
  end subroutine fix_group_center
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Shift the group center for the 'gid'-th group
!
  subroutine shift_group_center(gid, shift_vector)
!
    use M_group_id_setting, only : num_group_mem, group_mem    !(unchanged)
    use elses_mod_txp,      only : txp, typ, tzp               !(CHANGED!)
    implicit none
    integer,                intent(in)  :: gid
    real(DOUBLE_PRECISION), intent(in)  :: shift_vector(3)
    integer :: j, jsv, js
!
    do j = 1, num_group_mem(gid)
      jsv=group_mem(j, gid)
      js=jsv
      txp(js)=txp(js)+shift_vector(1)
      typ(js)=typ(js)+shift_vector(2)
      tzp(js)=tzp(js)+shift_vector(3)
    enddo
!
    call elses_gene_tx
!
  end subroutine shift_group_center
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Save structure with group centers
!
  subroutine save_w_group_centers
    use M_config,             only : config                  !(unchanged)
    use elses_mod_file_io,    only : vacant_unit             !(function)
    use elses_mod_phys_const, only : angst                   !(parameter)
    use M_group_id_setting,   only : num_groups              !(unchanged)
    use M_group_id_setting,   only : group_center            !(unchanged)
!
    implicit none
    integer :: gid, j, jsv
    integer :: n_atom
    integer :: lu
    integer :: ierr
    real(DOUBLE_PRECISION) :: mass_sum, mass
    integer ::fd
    character(len=*), parameter :: filename='test.xyz'
    character(len=*), parameter :: elem_name_of_center='X   ' 
    real(DOUBLE_PRECISION) :: cell_size_angst(3)
!
    fd=vacant_unit()
    lu=config%calc%distributed%log_unit
    n_atom = config%system%structure%natom
!
    cell_size_angst(1)=config%system%structure%unitcell%vectorA(1)*angst
    cell_size_angst(2)=config%system%structure%unitcell%vectorB(2)*angst
    cell_size_angst(3)=config%system%structure%unitcell%vectorC(3)*angst
!
    open(fd,file=filename)

    write(fd,'(i10,3f10.5)') n_atom+num_groups, cell_size_angst(1:3)
    write(fd,'(a)') 'sample with group centers'

    do jsv=1,n_atom
      write(fd,'(a4,3f10.5)') config%system%structure%vatom(jsv)%name, &
&       config%system%structure%vatom(jsv)%position(1:3)*angst
    enddo
!
    do gid = 1, num_groups
      write(fd,'(a4,3f10.5)') elem_name_of_center, &
&                             group_center(1:3,gid)*cell_size_angst(1:3)
    enddo
!
    close(fd)

  end subroutine save_w_group_centers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_gid_basic_routines

