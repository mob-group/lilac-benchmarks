!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
subroutine eigensort_val_asc(eval, evec, dof1, dof2)
      implicit none

      integer,intent(in)         :: dof1, dof2
      ! e.g. EVAL is an array of eigenvalues, EVEC is the corresponding eigenvector array.
      DOUBLE PRECISION,INTENT(INOUT)  :: EVAL(DOF1), EVEC(DOF1,DOF2)

      integer :: i,j,k
      DOUBLE PRECISION :: TMP 

      ! Step through the array, each time considering the subarray with indices ranging from i to the end of the array
      do i=1,dof1-1
         k=i
         tmp=eval(i)
         ! Loop through all eigenvalues further down the array, and find the largest.
         do j=i+1, dof1
            if (eval(j)>=tmp) then
               k=j
               tmp=eval(j)
            endif
         enddo
         ! If the largest one isn't already at the front of the sub-array that we're still searching:
         if (.not.k==i) then
            ! Swap the largest e-value with the one at the front of the sub-array
            eval(k)=eval(i)
            eval(i)=tmp
            ! Transpose the corresponding eigenvectors as well
            do j=1,dof1
               tmp = evec(j,i)
               evec(j,i) = evec(j,k)
               evec(j,k) = tmp
            enddo
         endif
      enddo
 end subroutine eigensort_val_asc
