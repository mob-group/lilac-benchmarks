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
module DLLModule
     use DataModule,only:IntegerList,RealList,NodeList,Node
     implicit none
     contains
     subroutine DLLinitR(l,value)
          implicit none
          type(RealList),pointer :: l
          double precision,intent(in) :: value
          if (associated(l)) then
               print *, 'DLLinitI> ERROR: list pointer is already associated!'; stop
          else
               allocate(l)
               l%f=>l
               l%b=>l
               l%value=value
          endif
     end subroutine DLLinitR
     subroutine DLLinitI(l,value)
          implicit none
          type(IntegerList),pointer :: l
          integer,intent(in) :: value
          if (associated(l)) then
               print *, 'DLLinitI> ERROR: list pointer is already associated!'; stop
          else
               allocate(l)
               l%f=>l
               l%b=>l
               l%value=value
          endif
     end subroutine DLLinitI
     subroutine DLLinitN(l,value)
          implicit none
          type(NodeList),pointer :: l
          type(Node),pointer :: value
          if (associated(l)) then
               print *, 'DLLinitN> ERROR: list pointer is already associated!'; stop
          else
               allocate(l)
               l%f=>l
               l%b=>l
               l%value=>value
          endif
     end subroutine DLLinitN
     subroutine DLLaddR(l,value)
          implicit none
          type(RealList),pointer :: l
          double precision,intent(in) :: value
          type(RealList),pointer :: new
          if (associated(l)) then
               if (.not.associated(l%f)) then
                    print *, 'DLLaddI> ERROR: forward pointer is not associated!'; stop
               elseif (.not.associated(l%b)) then
                    print *, 'DLLaddI> ERROR: backward pointer is not associated!'; stop
               else
                    allocate(new)
                    new%value=value
                    new%f=>l
                    new%b=>l%b
                    l%b%f=>new
                    l%b=>new
                    l=>new
                    nullify(new)
               endif
          else
               call DLLinitR(l,value)
          endif
     end subroutine DLLaddR
     subroutine DLLaddI(l,value)
          implicit none
          type(IntegerList),pointer :: l
          integer,intent(in) :: value
          type(IntegerList),pointer :: new
          if (associated(l)) then
               if (.not.associated(l%f)) then
                    print *, 'DLLaddI> ERROR: forward pointer is not associated!'; stop
               elseif (.not.associated(l%b)) then
                    print *, 'DLLaddI> ERROR: backward pointer is not associated!'; stop
               else
                    allocate(new)
                    new%value=value
                    new%f=>l
                    new%b=>l%b
                    l%b%f=>new
                    l%b=>new
                    l=>new
                    nullify(new)
               endif
          else
               call DLLinitI(l,value)
          endif
     end subroutine DLLaddI
     subroutine DLLaddN(l,value)
          implicit none
          type(NodeList),pointer :: l
          type(Node),pointer :: value
          type(NodeList),pointer :: new
          if (associated(l)) then
               if (.not.associated(l%f)) then
                    print *, 'DLLaddN> ERROR: forward pointer is not associated!'; stop
               elseif (.not.associated(l%b)) then
                    print *, 'DLLaddI> ERROR: backward pointer is not associated!'; stop
               else
                    allocate(new)
                    new%value=>value
                    new%f=>l
                    new%b=>l%b
                    l%b%f=>new
                    l%b=>new
                    l=>new
                    nullify(new)
               endif
          else
               call DLLinitN(l,value)
          endif
     end subroutine DLLaddN
     subroutine DLLprintR(l)
          implicit none
          type(RealList),pointer :: l
          integer :: i
          type(RealList),pointer :: s
          if (associated(l)) then
               print *, 'DLLprintR> List contains:'
               s=>l%f
               nullify(s%b)
               i=-1
               do
                    if (associated(s%b).or.i==-1) then
                         i=i+1
                         print *, "DLLprintR> ",i,s%value
                         s=>s%f
                    else
                         print *, "DLLprintR> List length =",i+1
                         l%f%b=>l
                         nullify(s)
                         exit
                    endif
               enddo
          else
               print *, 'DLLprintR> The list is empty'
          endif
     end subroutine DLLprintR
     subroutine DLLprintI(l)
          implicit none
          type(IntegerList),pointer :: l
          integer :: i
          type(IntegerList),pointer :: s
          if (associated(l)) then
               print *, 'DLLprintI> List contains:'
               s=>l%f
               nullify(s%b)
               i=-1
               do
                    if (associated(s%b).or.i==-1) then
                         i=i+1
                         print *, "DLLprintI> ",i,s%value
                         s=>s%f
                    else
                         print *, "DLLprintI> List length =",i+1
                         l%f%b=>l
                         nullify(s)
                         exit
                    endif
               enddo
          else
               print *, 'DLLprintI> The list is empty'
          endif
     end subroutine DLLprintI
     subroutine DLLprintN(l)
          implicit none
          type(NodeList),pointer :: l
          integer :: i
          type(NodeList),pointer :: s
          if (associated(l)) then
               print *, 'DLLprintN> List contains:'
               s=>l%f
               nullify(s%b)
               i=-1
               do
                    if (associated(s%b).or.i==-1) then
                         i=i+1
                         print *, "DLLprintN> ",i,s%value%index
                         s=>s%f
                    else
                         print *, "DLLprintN> List length =",i+1
                         l%f%b=>l
                         nullify(s)
                         exit
                    endif
               enddo
          else
               print *, 'DLLprintN> The list is empty'
          endif
     end subroutine DLLprintN
     subroutine DLLdelR(l)
          implicit none
          type(RealList),pointer :: l
          type(RealList),pointer :: s
          if (associated(l)) then
               s=>l%f
               nullify(l%f)
               if (associated(s%f)) then
                    s%b=>l%b
                    s%b%f=>s
                    nullify(l%b)
                    deallocate(l)
                    l=>s
                    nullify(s)
               else
                    nullify(l%b,s)
                    deallocate(l)
               endif
          endif
     end subroutine DLLdelR
     subroutine DLLdelI(l)
          implicit none
          type(IntegerList),pointer :: l
          type(IntegerList),pointer :: s
          if (associated(l)) then
               s=>l%f
               nullify(l%f)
               if (associated(s%f)) then
                    s%b=>l%b
                    s%b%f=>s
                    nullify(l%b)
                    deallocate(l)
                    l=>s
                    nullify(s)
               else
                    nullify(l%b,s)
                    deallocate(l)
               endif
          endif
     end subroutine DLLdelI
     subroutine DLLdelN(l)
          implicit none
          type(NodeList),pointer :: l
          type(NodeList),pointer :: s
          if (associated(l)) then
               s=>l%f
               nullify(l%f)
               if (associated(s%f)) then
                    s%b=>l%b
                    s%b%f=>s
                    nullify(l%b,l%value)
                    deallocate(l)
                    l=>s
                    nullify(s)
               else
                    nullify(l%b,s,l%value)
                    deallocate(l)
               endif
          endif
     end subroutine DLLdelN
     subroutine DLLsumR(l,r)
          implicit none
          type(RealList),pointer :: l
          integer :: i
          double precision,intent(out) :: r
          type(RealList),pointer :: s
          r=0.0d0
          if (associated(l)) then
               s=>l%f
               nullify(s%b)
               i=-1
               do
                    if (associated(s%b).or.i==-1) then
                         i=i+1
                         r=r+s%value
                         s=>s%f
                    else
                         l%f%b=>l
                         nullify(s)
                         exit
                    endif
               enddo
          endif
          !print *, 'DLLsumR> Sum =',r
     end subroutine DLLsumR
     subroutine DLLdivR(l,r)
          implicit none
          type(RealList),pointer :: l
          double precision,intent(in) :: r
          integer :: i
          type(RealList),pointer :: s
          if (associated(l)) then
               s=>l%f
               nullify(s%b)
               i=-1
               do
                    if (associated(s%b).or.i==-1) then
                         i=i+1
                         s%value=s%value/r
                         s=>s%f
                    else
                         l%f%b=>l
                         nullify(s)
                         exit
                    endif
               enddo
          endif
     end subroutine DLLdivR
     subroutine DLLmulR(l,r)
          implicit none
          type(RealList),pointer :: l
          double precision,intent(in) :: r
          integer :: i
          type(RealList),pointer :: s
          if (associated(l)) then
               s=>l%f
               nullify(s%b)
               i=-1
               do
                    if (associated(s%b).or.i==-1) then
                         i=i+1
                         s%value=s%value*r
                         s=>s%f
                    else
                         l%f%b=>l
                         nullify(s)
                         exit
                    endif
               enddo
          endif
     end subroutine DLLmulR
end module DLLModule
