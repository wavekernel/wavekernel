!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_output_matrices
!
  implicit none
  integer, allocatable :: called_first_h(:)
  integer, allocatable :: called_first_s(:)
  integer, allocatable :: called_first_e(:)
!
   private
!
! Public routines
!  public qm_output_matrices
   public output_levels_matrices
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Output for levels and matirces
!
  subroutine output_levels_matrices
    use M_lib_phys_const, only : ev4au  !(unchanged) 
    use M_config                        !(unchanged)
    use M_qm_output_levels, only : qm_output_levels !(routine)
    implicit none
!   character(len=1)  :: mode_for_matrix
    character(len=32) :: filename_wrk
    character(len=32) :: format_wrk
    real(8)           :: unit_conv
    logical           :: append_mode
    logical           :: active
    integer           :: ierr 
    integer           :: i_verbose, log_unit 
!
    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (config%calc%distributed%set) then
      if (i_verbose >=1) then
        if (log_unit > 0) then 
          write(log_unit,'(a)')'@@ output_levels_matrices...is skipped in the DST scheme'
        endif  
      endif   
      return
    endif   
!
    if (i_verbose >=1) then
      if (log_unit > 0) then 
        write(log_unit,'(a)')'@@ output_levels_matrices'
      endif  
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (config%output%matrix_hamiltonian%set) then
      filename_wrk=trim(config%output%matrix_hamiltonian%filename)
      format_wrk=trim(config%output%matrix_hamiltonian%format)
      if ( config%output%matrix_hamiltonian%unit == "eV" ) then
        unit_conv=ev4au
      else
        unit_conv=1.0d0
      endif   
      call set_active_status("H", active)
      write(*,*)'  matrix_hamiltonian:active=', active
      if (active) then
        call set_append_mode("H", append_mode) 
        call qm_output_matrices_mm("H", trim(filename_wrk), unit_conv, append_mode, format_wrk)
      endif   
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (config%output%matrix_overlap%set) then
      filename_wrk=trim(config%output%matrix_overlap%filename)
      format_wrk=trim(config%output%matrix_overlap%format)
      unit_conv=1.0d0
      call set_active_status("S", active)
      write(*,*)'  matrix_overlap:active=', active
      if (active) then
        call set_append_mode("S", append_mode) 
        call qm_output_matrices_mm("S", trim(filename_wrk), unit_conv, append_mode, format_wrk)
      endif   
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (config%output%eigen_level%set) then
      filename_wrk=trim(config%output%eigen_level%filename)
      format_wrk=trim(config%output%eigen_level%format)
      if ( config%output%eigen_level%unit == "eV" ) then
        unit_conv=ev4au
      else
        unit_conv=1.0d0
      endif   
      call set_active_status("E", active)
      write(*,*)'  eigen_level:active=', active
      if (active) then
        call set_append_mode("E", append_mode) 
        call qm_output_levels(trim(filename_wrk), unit_conv, append_mode, format_wrk)
      endif   
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  end subroutine output_levels_matrices
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set active status
!
  subroutine set_active_status(mode, active_mode)
    use M_config                        !(unchanged)
    use elses_mod_md_dat, only : itemd, itemdmx, final_iteration !(unchanged)
    implicit none
    character(len=*), intent(in) :: mode
    logical,         intent(out) :: active_mode
    integer                      :: ierr 
!
    active_mode = .false.
!
    select case(mode)
      case ("H")
        select case(config%output%matrix_hamiltonian%mode)
          case ("last") 
            if (itemd == itemdmx) active_mode = .true.
            if (final_iteration)  active_mode = .true.
          case ("periodic") 
            if (config%output%matrix_hamiltonian%interval /= 0) then
              if (mod(itemd-1,config%output%matrix_hamiltonian%interval) == 0 ) active_mode = .true.
            else
              write(*,*)'ERROR(set_active_status):config%matrix_hamiltonian%mode     =', config%output%matrix_hamiltonian%mode
              write(*,*)'ERROR(set_active_status):config%matrix_hamiltonian%interval =', config%output%matrix_hamiltonian%interval
              stop
            endif  
          case default
            write(*,*) 'ERROR(set_active_status):mode=',mode 
            stop
        end select   
      case ("S")
        select case(config%output%matrix_overlap%mode)
          case ("last") 
            if (itemd == itemdmx) active_mode = .true.
            if (final_iteration)  active_mode = .true.
          case ("periodic") 
            if (config%output%matrix_overlap%interval /= 0) then
              if (mod(itemd-1,config%output%matrix_overlap%interval) == 0 ) active_mode = .true.
            else
              write(*,*)'ERROR(set_active_status):config%matrix_overlap%mode     =', config%output%matrix_overlap%mode
              write(*,*)'ERROR(set_active_status):config%matrix_overlap%interval =', config%output%matrix_overlap%interval
              stop
            endif  
          case default
            write(*,*) 'ERROR(set_active_status):mode=',mode 
            stop
        end select   
      case ("E")
        select case(config%output%eigen_level%mode)
          case ("last") 
            if (itemd == itemdmx) active_mode = .true.
            if (final_iteration)  active_mode = .true.
          case ("periodic") 
            if (config%output%eigen_level%interval /= 0) then
              if (mod(itemd-1,config%output%eigen_level%interval) == 0 ) active_mode = .true.
            else
              write(*,*)'ERROR(set_active_status):config%output%eigen_level%mode     =', config%output%eigen_level%mode
              write(*,*)'ERROR(set_active_status):config%output%eigen_level%interval =', config%output%eigen_level%interval 
              stop
            endif  
          case default
            write(*,*) 'ERROR(set_active_status):mode=',mode 
            stop
        end select   
      case default
        write(*,*)'STOP(set_active_status)' 
        stop
    end select     
!
  end subroutine set_active_status
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set append mode
!
  subroutine set_append_mode(mode, append_mode)
    implicit none
    character(len=*), intent(in) :: mode
    logical,         intent(out) :: append_mode
    integer                      :: ierr 
!
    select case(mode)
      case ("H")
        if (.not. allocated(called_first_h)) then
          append_mode = .false.
          allocate(called_first_h(1), stat=ierr)
          if (ierr /=0 ) stop 'Alloc. error (qm_output_matrices):called_first_h'
        else
          append_mode = .true.
        endif   
      case ("S")
        if (.not. allocated(called_first_s)) then
          append_mode = .false.
          allocate(called_first_s(1), stat=ierr)
          if (ierr /=0 ) stop 'Alloc. error (qm_output_matrices):called_first_s'
        else
          append_mode = .true.
        endif   
      case ("E")
        if (.not. allocated(called_first_e)) then
          append_mode = .false.
          allocate(called_first_e(1), stat=ierr)
          if (ierr /=0 ) stop 'Alloc. error (qm_output_matrices):called_first_e'
        else
          append_mode = .true.
        endif   
      case default
        write(*,*)'STOP(set_append_mode' 
        stop
    end select   


  end subroutine set_append_mode
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Output for Hamiltonian or overlap matirces
!     into the Matrix-Market format
!
  subroutine qm_output_matrices_mm(mode, filename_wrk, unit_conv, append_mode, format_mode)
!
    use M_qm_domain ,  only : i_verbose    !(unchanged)
    use M_config,            only : config !(unchanged)
    use elses_mod_orb2,      only : js2j   !(unchanged)
    use elses_mod_js4jsv,    only : js4jsv !(unchanged)
    use elses_mod_md_dat,    only : itemd  !(unchanged)
!
    use elses_mod_file_io, only : vacant_unit !(function)
    use M_qm_domain ,  only : dhij, dsij, atm_element, nval, &
&                             jsv4jsd, njsd, noav, atm_element !(unchanged) 
!
    implicit none
    character(len=*), intent(in) :: mode
    character(len=*), intent(in) :: filename_wrk
    real(8),          intent(in) :: unit_conv
    logical,          intent(in) :: append_mode
    character(len=*), intent(in) :: format_mode
!
    integer      :: jj, prc_index, ierr, iunit, ict4h
    integer      :: js2, jsv2, nss2, nval2, ja2, jj2, jsd1
    integer      :: js1, jsv1, nss1, nval1, ja1, jj1
    integer      :: mat_size, num_of_non_zero
    real(kind=8) :: ddd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    ict4h=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Show the title (non-essential)
!
    if (i_verbose >= 1) then
      write(*,*)'@@ Plot a matrix into the matrix-market form : mode = ',trim(mode)
      write(*,*)'  output file = ',trim(filename_wrk)
    endif 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Trivial Checking 
!
    select case(trim(mode)) 
      case ('H') 
        if (.not. allocated(dhij)) then
          write(*,*)'INFO:(qm_output_matrices):dhij is not allocated-->skipped'
          return
        endif
      case ('S')
        if (.not. allocated(dsij)) then
          write(*,*)'INFO:(qm_output_matrices):dhij is not allocated-->skipped'
          return
        endif
      case default
      write(*,*)'ERROR(qm_output_matrices):mode=',trim(mode)
      stop
    end select
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open the output file
!
    iunit=vacant_unit()
!
    write(*,*)' iunit=',iunit
!
    if ( append_mode ) then
      if (i_verbose >= 1) then
        write(*,*)' The output file is initialized : ',trim(filename_wrk)
      endif
      open(iunit, file=filename_wrk, status='unknown', position='append')
    else
      open(iunit, file=filename_wrk, status='unknown')
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count up the matrix size
!
    jj=0
    do jsv2=1,noav
      jj=jj+nval(atm_element(jsv2))
    enddo
    mat_size=jj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   do prc_index=1,2
     jj=0
     do jsv2=1,noav
       nss2=atm_element(jsv2)
       nval2=nval(nss2)
       js2=js4jsv(jsv2)
       do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         js1=js4jsv(jsv1)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         js1=js4jsv(jsv1)
         do ja2=1,nval2
           jj2=js2j(ja2,js2)
           do ja1=1,nval1
             jj1=js2j(ja1,js1)
             if (trim(mode) == 'H') ddd=dhij(ja1,ja2,jsd1,jsv2)*unit_conv
             if (trim(mode) == 'S') ddd=dsij(ja1,ja2,jsd1,jsv2)*unit_conv
             if ( dabs(ddd) < 1.0d-16 ) cycle
             if ( trim(format_mode) == "MatrixMarket_sym" ) then
               if (jj1 < jj2) cycle
             endif   
             if (prc_index == 1) then
               jj=jj+1
             else
               write(iunit,'(2i10,f30.20)')jj1, jj2, ddd
             endif
          enddo  
        enddo  
       enddo  
     enddo  
     if (prc_index == 1) then 
       num_of_non_zero=jj
       if ( trim(format_mode) == "MatrixMarket_sym" ) then
         write(iunit,'(a)')'%%MatrixMarket matrix coordinate real symmetric'
       else   
         write(iunit,'(a)')'%%MatrixMarket matrix coordinate real general'
       endif   
       write(iunit,'(a,i10)')'% step_count=', config%system%structure%mdstep
       write(iunit,*)mat_size, mat_size, num_of_non_zero
       write(*,*)'  Matrix size, num of non-zero elements =',mat_size, num_of_non_zero
       write(*,*)'  Average num of non-zero elements per row =', dble(num_of_non_zero)/dble(mat_size)
     endif
   enddo
!
    close(iunit)

  end subroutine qm_output_matrices_mm
!
end module M_qm_output_matrices

