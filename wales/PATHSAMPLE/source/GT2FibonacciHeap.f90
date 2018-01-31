!
!  Copyright (C) 2005-2006 Semen Trygubenko
!  This file is part of GT.
!
!  GT is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  GT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with GT; see the file COPYING. If not, write to the
!  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
!  MA 02111-1307 USA
!
module FibonacciHeapModule
     use DataModule,only: Node,NodeArray
     implicit none
     integer,parameter :: MaxDegree=100
     integer :: n
     type(Node),pointer :: FHmin
     contains
     subroutine FHinit
          implicit none
          n=0
          nullify(FHmin)
     end subroutine FHinit
     subroutine FHinsert(x)
          implicit none
          type(Node),pointer :: x
          n=n+1
          x%degree=0
          nullify(x%parent,x%child)
          x%left=>x
          x%right=>x
          x%mark=.False.
          x%FH=.True.
          if (associated(FHmin)) then
               call FHconcatenate(FHmin,x)
               if (x%d<FHmin%d) FHmin=>x
          else
               FHmin=>x
          endif
     end subroutine FHinsert
     subroutine FHconcatenate(x,a)
          implicit none
          type(Node),pointer :: x,a,t
          t=>x%left
          x%left=>a%left
          a%left=>t
          nullify(t)
          a%left%right=>a
          x%left%right=>x
     end subroutine FHconcatenate
     subroutine FHremove(x,out)
          implicit none
          type(Node),pointer :: x,out
          if (x%index==x%right%index) then
               nullify(out)
          else
               x%left%right=>x%right
               x%right%left=>x%left
               out=>x%right
          endif
     end subroutine FHremove
     subroutine FHextractMin(z)
          implicit none
          type(Node),pointer :: z,child
          if (associated(FHmin)) then
               z=>FHmin
               if (associated(z%child)) then
                    child=>z%child
                    nullify(child%parent)
                    do
                         child=>child%left
                         if (associated(child%parent)) then
                              nullify(child%parent)
                         else
                              exit
                         endif
                    enddo
                    call FHconcatenate(z,z%child)
               endif
               call FHremove(z,FHmin)
               if (associated(FHmin)) call FHconsolidate()
               n=n-1
               z%FH=.False.
          else
               nullify(z)
          endif
     end subroutine FHextractMin
     subroutine FHlink(y,x)
          implicit none
          type(Node),pointer :: x,y,z
          call FHremove(y,z)
          nullify(z)
          y%parent=>x
          y%left=>y
          y%right=>y
          if (associated(x%child)) then
               call FHconcatenate(x%child,y)
          else
               x%child=>y
          endif
          x%degree=x%degree+1
          y%mark=.False.
     end subroutine FHlink
     subroutine FHconsolidate
          implicit none
          integer :: i,j,k,d
          type(NodeArray),allocatable :: A(:),B(:)
          type(Node),pointer :: x,y
          allocate(A(0:MaxDegree))
          do i=0,MaxDegree
               nullify(A(i)%p)
          enddo
          x=>FHmin
          j=x%index
          i=0
          do
               if (i==0.or..not.x%index==j) then
                    i=i+1
                    x=>x%right
               else
                    nullify(x)
                    exit
               endif
          enddo
          k=i
          allocate(B(k))
          x=>FHmin
          j=x%index
          i=0
          do
               if (i==0.or..not.x%index==j) then
                    i=i+1
                    B(i)%p=>x
                    x=>x%right
               else
                    nullify(x)
                    exit
               endif
          enddo
          do i=1,k
               x=>B(i)%p
               d=x%degree
               do
                    if (associated(A(d)%p)) then
                         y=>A(d)%p
                         if (x%d>y%d) then
                              y=>x
                              x=>A(d)%p
                         endif
                         call FHlink(y,x)
                         nullify(A(d)%p)
                         d=d+1
                         if (d>MaxDegree) then
                              print *, 'FHconsolidate> ERROR: MaxDegree parameter must be increased'; stop
                         endif
                    else
                         exit
                    endif
               enddo
               A(d)%p=>x
          enddo
          nullify(x)
          do i=1,k
               nullify(B(i)%p)
          enddo
          deallocate(B)
          nullify(FHmin)
          do i=0,MaxDegree
               if (associated(A(i)%p)) then
                    A(i)%p%left=>A(i)%p
                    A(i)%p%right=>A(i)%p
                    if (associated(FHmin)) then
                         call FHconcatenate(FHmin,A(i)%p)
                         if (A(i)%p%d<FHmin%d) FHmin=>A(i)%p
                    else
                         FHmin=>A(i)%p
                    endif
               endif
          enddo
          do i=0,MaxDegree
               nullify(A(i)%p)
          enddo
          deallocate(A)
     end subroutine FHconsolidate
     subroutine FHdecreaseKey(x,d)
          implicit none
          type(Node),pointer :: x,y
          integer,intent(in) :: d
          if (x%FH) then
               if (d>x%d) then
                    print *, 'FHdecreaseKey> ERROR: new key is greater than current key'; stop
               else
                    x%d=d
                    y=>x%parent
                    if (associated(y)) then
                         if (x%d<y%d) then
                              call FHcut(x,y)
                              call FHcascadingCut(y)
                         endif
                    endif
                    if (x%d<FHmin%d) FHmin=>x
               endif
          else
               print *, 'FHdecreaseKey> ERROR: node is not on the heap'; stop
          endif
     end subroutine FHdecreaseKey
     subroutine FHcut(x,y)
          implicit none
          type(Node),pointer :: x,y
          call FHremove(x,y%child)
          y%degree=y%degree-1
          x%left=>x
          x%right=>x
          call FHconcatenate(FHmin,x)
          nullify(x%parent)
          x%mark=.False.
     end subroutine FHcut
     recursive subroutine FHcascadingCut(y)
          implicit none
          type(Node),pointer :: y,z
          z=>y%parent
          if (associated(z)) then
               if (y%mark) then
                    call FHcut(y,z)
                    call FHcascadingCut(z)
               else
                    y%mark=.True.
               endif
          endif
     end subroutine FHcascadingCut
     subroutine FHincreaseKey(x,d)
          implicit none
          type(Node),pointer :: x,z
          integer,intent(in) :: d
          if (x%FH) then
               if (d<x%d) then
                    print *, 'FHdecreaseKey> ERROR: new key is smaller than current key'; stop
               else
                    call FHdecreaseKey(x,-huge(d))
                    call FHextractMin(z)
                    nullify(z)
                    x%d=d
                    call FHinsert(x)
               endif
          else
               print *, 'FHdecreaseKey> ERROR: node is not on the heap'; stop
          endif
     end subroutine FHincreaseKey
end module FibonacciHeapModule
