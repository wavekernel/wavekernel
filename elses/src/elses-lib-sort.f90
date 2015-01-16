module M_lib_sort
contains
  !***********************************************************************
  !     Heap Sort Program
  !         written by S. Yamamoto, Dec 12, 1998
  !         last revised Jun 23, 2003
  ! When there are some elements which have the same value,
  ! the original order isn't conserved.
  ! This program doesn't need an array for the temporary work area.
  !     Inputs
  !     N:   number of list elements
  !     value: given list of values
  !     Outputs (On exit, values are not altered.)
  !     index: (value(index(i)),i=1,N) is a non-decreasing series.
  !***********************************************************************
  !**   Heap sort is not suitable for a recursive algorithm, becase    ***
  !** a hierarchy of heap spreads over whole range of storage.         ***
  !***********************************************************************
  subroutine dhsort(N,value,index)
    implicit none
    integer, parameter :: DOUBLE_PRECISION=kind(0d0)
    integer N,i
    integer index(N),p,q,r,s,t
    intent(in)  :: N,value
    intent(out) :: index
    real(DOUBLE_PRECISION)::  value(N),u
    !     Making Heap ( Ordered Binary Tree )
    !     (i)Every node has greater value than those of its descendents.
    !     (ii)When index(p) is parent node, its child nodes are 
    !     index(2*p) and index(2*p+1).
    index(1)=1
    do i=2, N
       !     A heap is stored in index(1:i-1).
       !     appending the i-th node to the heap
       p=i
       !     If the value of the parent node is less than that of the current
       !     node, exchange their positions with each other.
       do while(p .gt. 1)
          q=p/2
          r=index(q)
          !     If all the ancestors of node p is greater than node p,
          !     quit inner loop.
          if(value(r).ge.value(i)) exit
          !     If not, exchange the positions with each other.
          index(p)=r
          p=q
       end do
       index(p)=i
    end do
    do i=N, 2, -1
       !     Moving the last node to the root.
       r=index(i)
       u=value(r)
       !     Deleting the root node.
       !     (Now, the size of heap decreases by one.)
       index(i)=index(1)
       p=1
       q=2
       !     To keep the rule of heap, the root node (previously,
       !     the last node) will be stored in index(p).
       !     index(q) and index(q+1) are the children of index(p).
       do while(q+1 .lt. i)
          s=index(q)
          t=index(q+1)
          if(value(t).ge.value(s)) then
             !     Node q+1 is greater than node q.
             !     Node q should have the smallest value among two children,
             !     because we exchange the node q with parent node p later.
             q=q+1
             s=t
          end if
          !     If all the descendents of the node p is smaller than the node p,
          !     quit inner loop.
          if(u.ge.value(s)) exit
          !     If not, exchange the positions.
          index(p)=s
          p=q
          q=2*q
       end do
       !     for the case that the last node have only one child
       if(q+1.eq.i) then
          s=index(q)
          if(u.lt.value(s)) then
             index(p)=s
             p=q
          end if
       end if
       index(p)=r
    end do
  end subroutine dhsort
end module M_lib_sort
