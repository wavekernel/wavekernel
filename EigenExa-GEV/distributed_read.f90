module distributed_read_m

  implicit none

contains

  subroutine distributed_read_header(filename, dim)
    character(*), intent(in) :: filename
    integer, intent(out) :: dim
    integer :: iunit, rows, cols, nnz
    character :: rep*10, field*7, symm*19

    open(iunit, file = trim(filename), status='old')
    call mminfo(iunit, rep, field, symm, rows, cols, nnz)
    dim = rows
    close(iunit)
  end subroutine distributed_read_header


  subroutine distributed_read(filename, desc, mat)
    character(*), intent(in) :: filename
    integer, intent(in) :: desc(9)
    double precision, intent(out) :: mat(:)
    integer, parameter :: iunit = 8
    integer :: rows, cols, nnz, info, i
    integer :: m_global, n_global, m_local, n_local
    integer :: nprow, npcol, myprow, mypcol, prow, pcol
    double precision :: val, work(1)
    character :: rep*10, field*7, symm*19
    character(len=256) :: linebuf
    integer, parameter :: dtype_ = 1, ctxt_ = 2, m_ = 3, n_ = 4, &
         mb_ = 5, nb_ = 6, rsrc_ = 7, csrc_ = 8, lld_ = 9
    ! Functions
    integer indxg2l, indxg2p

    call blacs_gridinfo(desc(ctxt_), nprow, npcol, myprow, mypcol)
    !print *, 'procs: ', nprow, npcol, myprow, mypcol

    open(iunit, file = trim(filename), status='old')
    call mminfo(iunit, rep, field, symm, rows, cols, nnz)
    rewind(iunit)

    mat(:) = 0.0d0
    do
      read (iunit, '(a)') linebuf
      if (linebuf(1:1) /= '%') then
        exit
      end if
    end do
    do i = 1, nnz
      read (iunit, *) m_global, n_global, val
      prow = indxg2p(m_global, desc(mb_), 0, desc(rsrc_), nprow)
      pcol = indxg2p(n_global, desc(nb_), 0, desc(csrc_), npcol)
      if (myprow == prow .and. mypcol == pcol) then
        !print *, prow, pcol, ": ", m_global, n_global, ": ", val
        m_local = indxg2l(m_global, desc(mb_), 0, desc(rsrc_), nprow)
        n_local = indxg2l(n_global, desc(nb_), 0, desc(csrc_), npcol)
        mat(desc(lld_) * (n_local - 1) + m_local) = val
      end if
      !  Transposed
      prow = indxg2p(n_global, desc(mb_), 0, desc(rsrc_), nprow)
      pcol = indxg2p(m_global, desc(nb_), 0, desc(csrc_), npcol)
      if (myprow == prow .and. mypcol == pcol) then
        m_local = indxg2l(n_global, desc(mb_), 0, desc(rsrc_), nprow)
        n_local = indxg2l(m_global, desc(nb_), 0, desc(csrc_), npcol)
        mat(desc(lld_) * (n_local - 1) + m_local) = val
      end if
    end do
    close(iunit)

    !call pdlaprnt(desc(m_), desc(n_), mat, 1, 1, desc, 0, 0, 'B', 6, work)
  end subroutine distributed_read
end module distributed_read_m
