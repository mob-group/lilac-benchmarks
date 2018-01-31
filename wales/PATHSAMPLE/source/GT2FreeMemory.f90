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
module FreeMemoryModule
     implicit none
     contains
     subroutine FreeMemory
          use DataModule,only:G,Nnodes
          use DLLModule,only:NodeList,RealList
          implicit none
          integer :: i
          type(NodeList),pointer :: A
          type(RealList),pointer :: P
          do i=1,Nnodes
               if (associated(G(i)%p%A)) then
                    A=>G(i)%p%A
                    P=>G(i)%p%P
                    nullify(A%b%f,P%b%f,A%value)
                    nullify(A%b,P%b)
                    do
                         if (associated(A%f)) then
                              G(i)%p%A=>A%f
                              G(i)%p%P=>P%f
                              nullify(A%f%b,P%f%b)
                              nullify(A%f,P%f)
                              deallocate(A,P)
                              A=>G(i)%p%A
                              P=>G(i)%p%P
                         else
                              deallocate(A,P)
                              exit
                         endif
                    enddo
               endif
               deallocate(G(i)%p)
          enddo
          deallocate(G)
     end subroutine FreeMemory
end module FreeMemoryModule
