!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

program elses_xml_generage_nrl_param
  implicit none

  type onsite_type
     real(8) a_s, b_s, c_s, d_s
     real(8) a_p, b_p, c_p, d_p
     real(8) a_t2g, b_t2g, c_t2g, d_t2g
     real(8) a_eg, b_eg, c_eg, d_eg
  end type onsite_type

  type hop_type
     real(8) e_sss, e_sps, e_pps, e_ppp, e_sds, e_pds, e_pdp, e_dds, e_ddp, e_ddd
     real(8) f_sss, f_sps, f_pps, f_ppp, f_sds, f_pds, f_pdp, f_dds, f_ddp, f_ddd
     real(8) fb_sss, fb_sps, fb_pps, fb_ppp, fb_sds, fb_pds, fb_pdp, fb_dds, fb_ddp, fb_ddd
     real(8) g_sss, g_sps, g_pps, g_ppp, g_sds, g_pds, g_pdp, g_dds, g_ddp, g_ddd
  end type hop_type

  type element_type
     integer num_orbital
     real(8) occupancy(3)
     real(8) mass
     real(8) lambda
     character(20) name
     type( onsite_type ) :: onsite
  end type element_type

  type pair_type
     real(8) cutoff
     real(8) rscreen
     type( hop_type ) ::  hamiltonian
     type( hop_type ) :: overlap
  end type pair_type

  type param_type
     integer num_species
     character(7) style
     type(element_type), pointer, dimension(:) :: element
     type(pair_type), pointer, dimension(:,:) :: pair
  end type param_type

  type(param_type) :: param

  call main

contains

  subroutine main

    implicit none
    integer :: i
    integer       :: iargc  
    character(len=50) flnm_in
    character(len=50) flnm_out
    character(len=50) species_name
    
    if( iargc() /= 3 ) then
       write(*,'(a)') "# Usage of this utility program"
       write(*,'(a)') "#   elses-xml-generate-nrl-param nrl-formatted-file xml-formatted-file species-name"
       stop
    end if
    
    call getarg(1,flnm_in)
    call getarg(2,flnm_out)
    call getarg(3,species_name)

    call nrlparam_load(flnm_in,species_name )

    do i = 1, param%num_species
       call nrlparam_save(flnm_out,i)
    end do

    deallocate( param%element )
    deallocate( param%pair )
    return
  end subroutine main

  subroutine nrlparam_load(filename,name)
    implicit none
    integer, parameter ::  fd = 1
    integer n
    character(len=*) filename
    character(len=*) name
    character(256) buffer

    open(fd,file=filename)

    ! read header
    read(fd,*) 
    read(fd,*) param%style

    read(fd,*) param%num_species
    n = param%num_species
    allocate( param%element(n) )
    allocate( param%pair(n,n) )
    param%element(1:n)%name = name

    read(fd,*) param%pair(1,1)%cutoff, param%pair(1,1)%rscreen
    read(fd,*) param%element(1)%num_orbital
    read(fd,*) param%element(1)%mass
    read(fd,*) param%element(1)%occupancy(1:3)
    read(fd,*) param%element(1)%lambda

    ! onsite energy
    read(fd,*) param%element(1)%onsite%a_s
    read(fd,*) param%element(1)%onsite%b_s
    read(fd,*) param%element(1)%onsite%c_s
    read(fd,*) param%element(1)%onsite%d_s

    read(fd,*) param%element(1)%onsite%a_p
    read(fd,*) param%element(1)%onsite%b_p
    read(fd,*) param%element(1)%onsite%c_p
    read(fd,*) param%element(1)%onsite%d_p

    read(fd,*) param%element(1)%onsite%a_t2g
    read(fd,*) param%element(1)%onsite%b_t2g
    read(fd,*) param%element(1)%onsite%c_t2g
    read(fd,*) param%element(1)%onsite%d_t2g

    read(fd,*) param%element(1)%onsite%a_eg
    read(fd,*) param%element(1)%onsite%b_eg
    read(fd,*) param%element(1)%onsite%c_eg
    read(fd,*) param%element(1)%onsite%d_eg

    ! hamiltonian
    read(fd,*) param%pair(1,1)%hamiltonian%e_sss
    read(fd,*) param%pair(1,1)%hamiltonian%f_sss
    read(fd,*) param%pair(1,1)%hamiltonian%fb_sss
    read(fd,*) param%pair(1,1)%hamiltonian%g_sss

    read(fd,*) param%pair(1,1)%hamiltonian%e_sps
    read(fd,*) param%pair(1,1)%hamiltonian%f_sps
    read(fd,*) param%pair(1,1)%hamiltonian%fb_sps
    read(fd,*) param%pair(1,1)%hamiltonian%g_sps

    read(fd,*) param%pair(1,1)%hamiltonian%e_pps
    read(fd,*) param%pair(1,1)%hamiltonian%f_pps
    read(fd,*) param%pair(1,1)%hamiltonian%fb_pps
    read(fd,*) param%pair(1,1)%hamiltonian%g_pps

    read(fd,*) param%pair(1,1)%hamiltonian%e_ppp
    read(fd,*) param%pair(1,1)%hamiltonian%f_ppp
    read(fd,*) param%pair(1,1)%hamiltonian%fb_ppp
    read(fd,*) param%pair(1,1)%hamiltonian%g_ppp

    read(fd,*) param%pair(1,1)%hamiltonian%e_sds
    read(fd,*) param%pair(1,1)%hamiltonian%f_sds
    read(fd,*) param%pair(1,1)%hamiltonian%fb_sds
    read(fd,*) param%pair(1,1)%hamiltonian%g_sds

    read(fd,*) param%pair(1,1)%hamiltonian%e_pds
    read(fd,*) param%pair(1,1)%hamiltonian%f_pds
    read(fd,*) param%pair(1,1)%hamiltonian%fb_pds
    read(fd,*) param%pair(1,1)%hamiltonian%g_pds

    read(fd,*) param%pair(1,1)%hamiltonian%e_pdp
    read(fd,*) param%pair(1,1)%hamiltonian%f_pdp
    read(fd,*) param%pair(1,1)%hamiltonian%fb_pdp
    read(fd,*) param%pair(1,1)%hamiltonian%g_pdp

    read(fd,*) param%pair(1,1)%hamiltonian%e_dds
    read(fd,*) param%pair(1,1)%hamiltonian%f_dds
    read(fd,*) param%pair(1,1)%hamiltonian%fb_dds
    read(fd,*) param%pair(1,1)%hamiltonian%g_dds

    read(fd,*) param%pair(1,1)%hamiltonian%e_ddp
    read(fd,*) param%pair(1,1)%hamiltonian%f_ddp
    read(fd,*) param%pair(1,1)%hamiltonian%fb_ddp
    read(fd,*) param%pair(1,1)%hamiltonian%g_ddp

    read(fd,*) param%pair(1,1)%hamiltonian%e_ddd
    read(fd,*) param%pair(1,1)%hamiltonian%f_ddd
    read(fd,*) param%pair(1,1)%hamiltonian%fb_ddd
    read(fd,*) param%pair(1,1)%hamiltonian%g_ddd

    ! overlap
    read(fd,*) param%pair(1,1)%overlap%e_sss
    read(fd,*) param%pair(1,1)%overlap%f_sss
    read(fd,*) param%pair(1,1)%overlap%fb_sss
    read(fd,*) param%pair(1,1)%overlap%g_sss

    read(fd,*) param%pair(1,1)%overlap%e_sps
    read(fd,*) param%pair(1,1)%overlap%f_sps
    read(fd,*) param%pair(1,1)%overlap%fb_sps
    read(fd,*) param%pair(1,1)%overlap%g_sps

    read(fd,*) param%pair(1,1)%overlap%e_pps
    read(fd,*) param%pair(1,1)%overlap%f_pps
    read(fd,*) param%pair(1,1)%overlap%fb_pps
    read(fd,*) param%pair(1,1)%overlap%g_pps

    read(fd,*) param%pair(1,1)%overlap%e_ppp
    read(fd,*) param%pair(1,1)%overlap%f_ppp
    read(fd,*) param%pair(1,1)%overlap%fb_ppp
    read(fd,*) param%pair(1,1)%overlap%g_ppp

    read(fd,*) param%pair(1,1)%overlap%e_sds
    read(fd,*) param%pair(1,1)%overlap%f_sds
    read(fd,*) param%pair(1,1)%overlap%fb_sds
    read(fd,*) param%pair(1,1)%overlap%g_sds

    read(fd,*) param%pair(1,1)%overlap%e_pds
    read(fd,*) param%pair(1,1)%overlap%f_pds
    read(fd,*) param%pair(1,1)%overlap%fb_pds
    read(fd,*) param%pair(1,1)%overlap%g_pds

    read(fd,*) param%pair(1,1)%overlap%e_pdp
    read(fd,*) param%pair(1,1)%overlap%f_pdp
    read(fd,*) param%pair(1,1)%overlap%fb_pdp
    read(fd,*) param%pair(1,1)%overlap%g_pdp

    read(fd,*) param%pair(1,1)%overlap%e_dds
    read(fd,*) param%pair(1,1)%overlap%f_dds
    read(fd,*) param%pair(1,1)%overlap%fb_dds
    read(fd,*) param%pair(1,1)%overlap%g_dds

    read(fd,*) param%pair(1,1)%overlap%e_ddp
    read(fd,*) param%pair(1,1)%overlap%f_ddp
    read(fd,*) param%pair(1,1)%overlap%fb_ddp
    read(fd,*) param%pair(1,1)%overlap%g_ddp

    read(fd,*) param%pair(1,1)%overlap%e_ddd
    read(fd,*) param%pair(1,1)%overlap%f_ddd
    read(fd,*) param%pair(1,1)%overlap%fb_ddd
    read(fd,*) param%pair(1,1)%overlap%g_ddd

    close(fd)
    return
  end subroutine nrlparam_load

  subroutine nrlparam_save(filename,is)
    implicit none
    integer is
    integer, parameter ::  fd = 1
    character(len=*) filename

    open(fd,file=filename)

    write(fd,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'

    write(fd,'(a,a,a)') '<element name="', trim(param%element(is)%name),'">'

    write(fd,*) '<classical_mechanics>'
    write(fd,'(a20,e23.15,a7)') '  <mass unit="a.u.">', param%element(is)%mass, '</mass>'
    write(fd,'(a16,e23.15,a15)') '  <ionic_charge>', 0.d0, '</ionic_charge>'
    write(fd,*) '</classical_mechanics>'

    write(fd,*) '<tight_binding_pair type="self">'
    write(fd,'(a20,e23.15,a7)') '  <rcut unit="a.u.">', param%pair(is,is)%cutoff, '</rcut>'
    write(fd,'(a11,e23.15,a10)') '  <rscreen>', param%pair(is,is)%rscreen, '</rscreen>'
    write(fd,'(a11,i5,a10)') '  <num_orb>', param%element(is)%num_orbital, '</num_orb>'
    write(fd,'(a13,3e23.15,a12)') '  <occupancy>', param%element(is)%occupancy(1:3), '</occupancy>'
    write(fd,'(a22,e23.15,a12)') '  <lambda unit="a.u.">', param%element(is)%lambda, '</lambda>'
    write(fd,'(a32,4e23.15,a9)') '  <onsite orbital="s" unit="Ry">', param%element(is)%onsite%a_s, &
         param%element(is)%onsite%b_s, param%element(is)%onsite%c_s, param%element(is)%onsite%d_s, '</onsite>'
    write(fd,'(a32,4e23.15,a9)') '  <onsite orbital="p" unit="Ry">', param%element(is)%onsite%a_p, &
         param%element(is)%onsite%b_p, param%element(is)%onsite%c_p, param%element(is)%onsite%d_p, '</onsite>'
    write(fd,'(a34,4e23.15,a9)') '  <onsite orbital="t2g" unit="Ry">', param%element(is)%onsite%a_t2g, &
         param%element(is)%onsite%b_t2g, param%element(is)%onsite%c_t2g, param%element(is)%onsite%d_t2g, '</onsite>'
    write(fd,'(a34,4e23.15,a9)') '  <onsite orbital="eg" unit="Ry">', param%element(is)%onsite%a_eg, &
         param%element(is)%onsite%b_eg, param%element(is)%onsite%c_eg, param%element(is)%onsite%d_eg, '</onsite>'

    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="sss" unit="Ry">', param%pair(is,is)%hamiltonian%e_sss, &
         param%pair(is,is)%hamiltonian%f_sss, param%pair(is,is)%hamiltonian%fb_sss, &
         param%pair(is,is)%hamiltonian%g_sss, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="sps" unit="Ry">', param%pair(is,is)%hamiltonian%e_sps, &
         param%pair(is,is)%hamiltonian%f_sps, param%pair(is,is)%hamiltonian%fb_sps, &
         param%pair(is,is)%hamiltonian%g_sps, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="pps" unit="Ry">', param%pair(is,is)%hamiltonian%e_pps, &
         param%pair(is,is)%hamiltonian%f_pps, param%pair(is,is)%hamiltonian%fb_pps, &
         param%pair(is,is)%hamiltonian%g_pps, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="ppp" unit="Ry">', param%pair(is,is)%hamiltonian%e_ppp, &
         param%pair(is,is)%hamiltonian%f_ppp, param%pair(is,is)%hamiltonian%fb_ppp, &
         param%pair(is,is)%hamiltonian%g_ppp, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="sds" unit="Ry">', param%pair(is,is)%hamiltonian%e_sds, &
         param%pair(is,is)%hamiltonian%f_sds, param%pair(is,is)%hamiltonian%fb_sds, &
         param%pair(is,is)%hamiltonian%g_sds, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="pds" unit="Ry">', param%pair(is,is)%hamiltonian%e_pds, &
         param%pair(is,is)%hamiltonian%f_pds, param%pair(is,is)%hamiltonian%fb_pds, &
         param%pair(is,is)%hamiltonian%g_pds, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="pdp" unit="Ry">', param%pair(is,is)%hamiltonian%e_pdp, &
         param%pair(is,is)%hamiltonian%f_pdp, param%pair(is,is)%hamiltonian%fb_pdp, &
         param%pair(is,is)%hamiltonian%g_pdp, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="dds" unit="Ry">', param%pair(is,is)%hamiltonian%e_dds, &
         param%pair(is,is)%hamiltonian%f_dds, param%pair(is,is)%hamiltonian%fb_dds, &
         param%pair(is,is)%hamiltonian%g_dds, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="ddp" unit="Ry">', param%pair(is,is)%hamiltonian%e_ddp, &
         param%pair(is,is)%hamiltonian%f_ddp, param%pair(is,is)%hamiltonian%fb_ddp, &
         param%pair(is,is)%hamiltonian%g_ddp, '</hamiltonian>'
    write(fd,'(a36,4e23.15,a14)') '  <hamiltonian bond="ddd" unit="Ry">', param%pair(is,is)%hamiltonian%e_ddd, &
         param%pair(is,is)%hamiltonian%f_ddd, param%pair(is,is)%hamiltonian%fb_ddd, &
         param%pair(is,is)%hamiltonian%g_ddd, '</hamiltonian>'

    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="sss" unit="Ry">', param%pair(is,is)%overlap%e_sss, &
         param%pair(is,is)%overlap%f_sss, param%pair(is,is)%overlap%fb_sss, &
         param%pair(is,is)%overlap%g_sss, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="sps" unit="Ry">', param%pair(is,is)%overlap%e_sps, &
         param%pair(is,is)%overlap%f_sps, param%pair(is,is)%overlap%fb_sps, &
         param%pair(is,is)%overlap%g_sps, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="pps" unit="Ry">', param%pair(is,is)%overlap%e_pps, &
         param%pair(is,is)%overlap%f_pps, param%pair(is,is)%overlap%fb_pps, &
         param%pair(is,is)%overlap%g_pps, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="ppp" unit="Ry">', param%pair(is,is)%overlap%e_ppp, &
         param%pair(is,is)%overlap%f_ppp, param%pair(is,is)%overlap%fb_ppp, &
         param%pair(is,is)%overlap%g_ppp, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="sds" unit="Ry">', param%pair(is,is)%overlap%e_sds, &
         param%pair(is,is)%overlap%f_sds, param%pair(is,is)%overlap%fb_sds, &
         param%pair(is,is)%overlap%g_sds, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="pds" unit="Ry">', param%pair(is,is)%overlap%e_pds, &
         param%pair(is,is)%overlap%f_pds, param%pair(is,is)%overlap%fb_pds, &
         param%pair(is,is)%overlap%g_pds, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="pdp" unit="Ry">', param%pair(is,is)%overlap%e_pdp, &
         param%pair(is,is)%overlap%f_pdp, param%pair(is,is)%overlap%fb_pdp, &
         param%pair(is,is)%overlap%g_pdp, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="dds" unit="Ry">', param%pair(is,is)%overlap%e_dds, &
         param%pair(is,is)%overlap%f_dds, param%pair(is,is)%overlap%fb_dds, &
         param%pair(is,is)%overlap%g_dds, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="ddp" unit="Ry">', param%pair(is,is)%overlap%e_ddp, &
         param%pair(is,is)%overlap%f_ddp, param%pair(is,is)%overlap%fb_ddp, &
         param%pair(is,is)%overlap%g_ddp, '</overlap>'
    write(fd,'(a32,4e23.15,a10)') '  <overlap bond="ddd" unit="Ry">', param%pair(is,is)%overlap%e_ddd, &
         param%pair(is,is)%overlap%f_ddd, param%pair(is,is)%overlap%fb_ddd, &
         param%pair(is,is)%overlap%g_ddd, '</overlap>'

    write(fd,*) '</tight_binding_pair>'

    write(fd,'(a)') '</element>'
    close(fd)

    return
  end subroutine nrlparam_save
end program elses_xml_generage_nrl_param
