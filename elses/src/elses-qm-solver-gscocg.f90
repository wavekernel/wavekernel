!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_solver_gscocg
!
!
  implicit none
!  
  private
!
! Public routines
  public qm_solver_gscocg
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine qm_solver_gscocg
!
    use elses_mod_orb2,     only : js2j ! (unchanged)
    use M_qm_domain,        only : noav, atm_element, nval ! (unchanged)
    use elses_mod_noav,     only : noav ! (unchanged)
    use M_qm_projection,    only : proj_init_end, proj_get_mat_size, &
&                                  proj_get_list, get_interac_list_num_proj_atom, convert_dm_and_edm ! (routines)
    use M_la_gscocg_main,   only : gscocg_main ! (routines)
!
    implicit none
    integer :: atm_index, orb_index, mat_size, num_energy_points, j_src
    integer :: num_atom_proj, m_int
    integer :: nss2, nval2, m, n, ierr, i, k
    integer :: prc_index
    integer :: imode
    real(8) :: sigma
    real(8) :: elec_num,miu,miu1,miu2,bieps,Nelec,Eelec,nf,nf1,nf2,tau,dk,temp !mended by Teng

    integer, allocatable :: jsv4jsk(:)
    integer, allocatable :: jsk4jsv(:)
    integer, allocatable :: jjkset(:)
!
    complex(8), allocatable  :: b(:)
    complex(8)               :: z_seed
    complex(8), allocatable  :: z_shift(:)
    real(8), allocatable     :: nj(:)
    real(8), allocatable     :: nt(:)
    real(8), allocatable     :: dos(:,:)
    real(8), allocatable     :: denma(:,:), edenma(:,:)
!
    write(*,*)'@@ qm_solver_gscocg'
!
    num_energy_points=3001
    z_seed=dcmplx(0.0D0,0.005D0) 
    dk=0.001d0
!
    n=num_energy_points
    allocate (nt(n),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (z_shift(n),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    do i=1,n
       sigma= (-1.0d0-dreal(z_seed))+dble(i-1)*dk
       z_shift(i)=z_seed+sigma
    enddo
!
    imode=1
    call proj_init_end(imode)
!    
    miu=0.0d0 ! dummy value for chemical potential (non essential)
!
    do prc_index=1,2
!
      Nelec=0.0d0 
      nt(:)=0.0d0
!     do atm_index=1,2
      do atm_index=1,noav
!
        call get_elec_num_for_atom(atm_index,elec_num)
        Nelec=Nelec+elec_num ! summing up of electron numbers
!
        call proj_get_mat_size(atm_index,mat_size,num_atom_proj)
        write(*,*)'atm_index     =',atm_index
        write(*,*)'mat_size      =',mat_size
        write(*,*)'num_atom_proj =',num_atom_proj
        nss2=atm_element(atm_index)
        nval2=nval(nss2)
        write(*,*)'nss2          =',nss2
        write(*,*)'nval2         =',nval2
!
        m=mat_size
        n=num_energy_points
!
        allocate (jsv4jsk(num_atom_proj),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        allocate (jsk4jsv(noav),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        allocate (jjkset(num_atom_proj),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        call proj_get_list(atm_index,jsv4jsk,jsk4jsv,jjkset)
!
        allocate (b(m),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
!
        allocate (nj(n),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        call get_interac_list_num_proj_atom(atm_index,m_int)
!           ---> get number of atom in the interaction list : m_int
!
        allocate (dos(m_int,n),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        do orb_index=1,nval2
          write(*,*)"atom=",atm_index,"orbit=",orb_index 
          j_src=js2j(orb_index,atm_index)
          b(:)=0.0d0
          b(orb_index)=1.0d0
          call gscocg_main(j_src,jsv4jsk,jsk4jsv,jjkset,b,z_seed,z_shift,nj,dos)
          nt(:)=nt(:)+nj(:)
        enddo   
!
        do i=1,n
          write(20,*)dreal(z_shift(i)),nt(i)
        enddo
!
        if (prc_index == 2) then ! (calculation of density matrix and energy denxity matrix)
          allocate ( denma(m_int,nval2),stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in alloc'
          allocate (edenma(m_int,nval2),stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in alloc'
          do orb_index=1,nval2
            do i=1,m_int
              do k=1,n
               denma(i, orb_index) = denma(i,orb_index)+ dos(i,k)*Fermi(dreal(z_shift(k)),miu,tau)*dk
               edenma(i,orb_index)= edenma(i,orb_index)+dreal(z_shift(k))*dos(i,k)*Fermi(dreal(z_shift(k)),miu,tau)*dk
              enddo ! (loop end for k)
              write(90+orb_index,*)i,denma(i,orb_index),edenma(i,orb_index)
            enddo ! (loop end for i)
            call convert_dm_and_edm(atm_index,orb_index,denma(:,orb_index),edenma(:,orb_index),jsv4jsk,jsk4jsv,jjkset)
          enddo ! ( loop end for orb_index )  
          deallocate(denma,stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in dealloc'
          deallocate(edenma,stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in dealloc'
        endif  
!
!
        deallocate(jsv4jsk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
        deallocate(jsk4jsv,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
        deallocate(jjkset,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
        deallocate(b,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'

        deallocate(dos,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
        deallocate(nj,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      enddo ! loop end for atm_index 
!
!
    if (prc_index == 1) then ! (calculation of chemical potential)
      write(*,*)'Nelec=',Nelec
!!!!!!mended by Teng 2009 Nov.
!!!!! solve the chemical potential with the bisection method
 

      Nelec=Nelec/2.0d0
      miu1=-1.0D+05;miu2=+1.0D+05
      tau=0.01D0;bieps=1.0D-10
            
      nf1=0.0D0;nf2=0.0D0
      do i=1,n
         nf1=nf1+nt(i)*Fermi(dreal(z_shift(i)),miu1,tau)*dk
         nf2=nf2+nt(i)*Fermi(dreal(z_shift(i)),miu2,tau)*dk
      enddo
      nf1=nf1-Nelec
      nf2=nf2-Nelec
      write(*,*)nf1,nf2

      if(nf1*nf2.ge.0.0D0) stop 'Abort: reset initial'
      do while((miu2-miu1).gt.bieps)
         miu=0.5D0*(miu1+miu2)
     
         nf=0.0D0
         do i=1,n
            nf=nf+nt(i)*Fermi(dreal(z_shift(i)),miu,tau)*dk
         enddo
         nf=nf-Nelec

         if(nf.le.0.0D0)then
            miu1=miu
         else
            miu2=miu
         endif
      enddo
      write(*,*) "the chemical potential is ",miu


!!!the end of the chemical potential solver

    endif   
!
    enddo ! loop end for prc_index

    deallocate(z_shift,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
    deallocate(nt,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
!   write(*,*)'Stop by T. Hoshi'
!   stop
!
    imode=2
    call proj_init_end(imode)
!
  end subroutine qm_solver_gscocg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function Fermi(energy,miu,tau)
    implicit none
    real(8) :: energy,miu,fvalue,tau, Fermi
    Fermi=1.0d0/(1.0d0+exp((energy-miu)/tau))
    !Fermi=1.0d0/(1.0d0+exp((energy-miu)/tau))
  end function Fermi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ subroutine for obtaining the valence elec. num. for the given atom
!
  subroutine get_elec_num_for_atom(atm_index,elec_num)
!
    use elses_mod_val_elec,     only : val_elec_atm ! (unchanged)
    use M_qm_domain,            only : atm_element  ! (unchanged)
!
    implicit none
    integer,          intent(in)  :: atm_index
    real(8),          intent(out) :: elec_num
    integer :: nss, size1, size2
!
    size1=size(atm_element,1)
    size2=size(val_elec_atm,1)
!
    if ((atm_index <= 0) .or. (atm_index > size1)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'atm_index=',atm_index
      stop
    endif   
!
    nss=atm_element(atm_index)
!
    if ((nss <= 0) .or. (nss > size2)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'nss=',nss
      stop
    endif   
!
    elec_num=val_elec_atm(nss)
!
    if ((elec_num <= 1.0d-10) .or. (elec_num >= 1000.0d0)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'  elec_num=',elec_num
      stop
    endif   
!
  end subroutine get_elec_num_for_atom
!
end module M_qm_solver_gscocg

