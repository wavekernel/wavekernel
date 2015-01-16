module M_group_id_setting
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  integer              :: num_groups        ! number of groups
  integer, allocatable :: num_group_mem(:)  ! number of group members --> num_group_mem(num_groups)
  integer, allocatable :: group_mem(:,:)    ! group members or atom id's in each group
!
  real(DOUBLE_PRECISION), allocatable :: group_center(:,:)      ! group center 
  real(DOUBLE_PRECISION), allocatable :: group_center_ini(:,:)  ! initial group center
!
  private
  public :: num_groups
  public :: num_group_mem
  public :: group_mem
  public :: group_center
  public :: group_center_ini
!
  public :: ini_group_id
  public :: set_weight_center_sub
!
  contains
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Show group_id information
!
  subroutine ini_group_id
    use M_config, only: config  ! unchanged
!
    implicit none
    integer :: lu
!
    lu=config%calc%distributed%log_unit
!
    if ( .not. config%calc%use_group_id ) return
!
    if (lu > 0) write(lu,*) '@@ init_group_id: use group id =', config%calc%use_group_id
!
    call set_group_member
!      ----> set              : num_groups
!      ----> allocate and set : num_group_mem(num_groups)
!      ----> allocate and set : group_mem(max_group_mem,num_groups)
!
    call set_weight_center
!      ----> allocate and set : group_center(3,num_groups)
!
    call set_weight_center_ini
!      ----> allocate and set : group_center_ini(3,num_groups)
!
    call output_w_group_centers
!
!   call add_constraint_w_groups
!
  end subroutine ini_group_id
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Plot the output file with group centers
!
  subroutine output_w_group_centers
    use M_config,             only : config                    !(unchanged)
    use elses_mod_file_io,    only : vacant_unit             !(function)
    use elses_mod_phys_const, only : angst                   !(parameter)
!   use M_qm_domain,        only : atm_position, atm_element !(unchanged)
!   use elses_mod_mass,     only : awt                       !(unchanged)
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

  end subroutine output_w_group_centers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set the weight center for each group : COMPATIBLE ONLY TO THE NON-PERIODIC CASES
!
  subroutine set_weight_center
    use M_config,           only : config                    !(unchanged)
!
    implicit none
    integer :: gid, j, jsv
    integer :: n_atom
    integer :: lu
    integer :: ierr
    real(DOUBLE_PRECISION) :: mass_sum, mass
!
    lu=config%calc%distributed%log_unit
    n_atom = config%system%structure%natom
!
    if (lu > 0) write(lu,'(a)')'@@ set_weight_center'
!
    ierr=0
    if (config%system%boundary%periodic_x) ierr=1
    if (config%system%boundary%periodic_y) ierr=1
    if (config%system%boundary%periodic_z) ierr=1
    if (ierr == 1) then
      write(*,*)'ERROR(set_weight_center):Periodic cases are not supported now'
      stop
    endif
!
    allocate (group_center(3,num_groups),stat=ierr)
    if ( ierr /= 0 ) stop 'Alloc. Error:group_cener'
    group_center(:,:)=0.0d0
!
    do gid = 1, num_groups
      call set_weight_center_sub(gid,group_center(1:3,gid))
    enddo
!
    do gid = 1, num_groups
      if (lu > 0) write(lu,'(a,i10,3f10.5)')'gid, group center =', gid, group_center(1:3,gid)
    enddo
!
  end subroutine set_weight_center
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set group_center_ini
!
  subroutine set_weight_center_ini
    use M_config,           only : config                    !(unchanged)
!
    implicit none
    integer :: gid, j, jsv
    integer :: n_atom
    integer :: lu
    integer :: ierr
    real(DOUBLE_PRECISION) :: mass_sum, mass
!
    allocate (group_center_ini(3,num_groups),stat=ierr)
    if ( ierr /= 0 ) stop 'Alloc. Error:group_cener_ini'
    group_center_ini(:,:)=group_center(:,:)
!
  end subroutine set_weight_center_ini
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set the weight center for 'gid'-th group : COMPATIBLE ONLY TO THE NON-PERIODIC CASES
!
  subroutine set_weight_center_sub(gid, group_center_wrk)
    use M_config,           only : config                    !(unchanged)
    use M_qm_domain,        only : atm_position, atm_element !(unchanged)
    use elses_mod_mass,     only : awt                       !(unchanged)
!
    implicit none
    integer,                intent(in)  :: gid
    real(DOUBLE_PRECISION), intent(out) :: group_center_wrk(3)
    integer :: j, jsv
    integer :: n_atom
    integer :: ierr
    real(DOUBLE_PRECISION) :: mass, mass_sum
!
    n_atom = config%system%structure%natom
!
    mass_sum=0.0d0
    group_center_wrk(1:3)=0.0d0
    do j = 1, num_group_mem(gid)
      jsv=group_mem(j, gid)
      mass=awt(atm_element(jsv))
      mass_sum=mass_sum+mass
      if ( (jsv < 1) .or. (jsv > n_atom) ) then
        write(*,*)'ERROR(ini_group_id):jsv=',jsv
        stop
      endif
      group_center_wrk(1:3)=group_center_wrk(1:3)+mass*atm_position(1:3,jsv)
    enddo
    group_center_wrk(1:3)=group_center_wrk(1:3)/mass_sum
!
  end subroutine set_weight_center_sub
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set group member 
!
  subroutine set_group_member
    use M_config, only: config  ! unchanged
!
    implicit none
    integer :: i, j, gid
    integer :: n_atom
    integer :: max_group_mem
    integer :: lu
    integer :: ierr
!
    lu=config%calc%distributed%log_unit
!
    if ( .not. config%calc%use_group_id ) return
!
    n_atom = config%system%structure%natom
!
    if (lu > 0) write(lu,'(a,i10)')'@@ set_group_member: n_atom = ', n_atom
!
    num_groups=maxval(config%system%structure%vatom(:)%group_id,1)
    if (lu > 0) write(lu,'(a,i10)')' number of groups =', num_groups
!
    if ( (num_groups < 1) .or. (num_groups > n_atom ) ) then
      write(*,*)'ERROR:number of groups =', num_groups
      stop
    endif
!
    allocate (num_group_mem(num_groups),stat=ierr)
    if ( ierr /= 0 ) stop 'Alloc. Error:num_group_mem'
    num_group_mem(:)=0
!
    do i=1,n_atom
      gid = config%system%structure%vatom(i)%group_id
      if ( (gid < 1) .or. (gid > num_groups) ) then
        write(*,*)'ERROR:gid =', gid
        stop
      endif   
      num_group_mem(gid)=num_group_mem(gid)+1
    enddo
!
    do gid = 1, num_groups
      if (lu > 0) write(lu,'(a,2i10)') ' gid, num_group_mem = ', gid, num_group_mem(gid)
    enddo
!
    max_group_mem = maxval(num_group_mem, 1)
    if (lu > 0) write(lu,'(a,i10)') '  max_group_mem =', max_group_mem 
!
    allocate (group_mem(max_group_mem, num_groups), stat=ierr)
    if ( ierr /= 0 ) stop 'Alloc. Error:group_mem'
    group_mem(:,:)=-1
!
    num_group_mem(:)=0
    do i=1,n_atom
      gid = config%system%structure%vatom(i)%group_id
      num_group_mem(gid)=num_group_mem(gid)+1
      group_mem(num_group_mem(gid), gid)=i
    enddo
!
    do gid = 1, num_groups
      do j = 1, num_group_mem(gid)
        i=group_mem(j, gid)
        if ( (i < 1) .or. (i > n_atom) ) then
          write(*,*)'ERROR(ini_group_id):i=',i
          stop
        endif
        if (lu > 0) write(lu,'(a,3i10)') 'group mem : gid, mem_id, atom_id=', gid, j, i
      enddo
    enddo
!
  end subroutine set_group_member
!
end module M_group_id_setting
