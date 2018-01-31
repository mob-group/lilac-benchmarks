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

! SAT&TVB  Mon Apr 18 21:08:26 BST 2005
MODULE BinaryIO
     implicit none
     contains

     subroutine ReadInpFile(Q) 
          use key,only: filth2
          use commons,only: natoms
          use porfuncs
          implicit none
          double precision,dimension(:),intent(inout)::Q

          integer :: i,j
          character(len=256) :: filename

          ! generate filename
          if (filth2==0) then
               filename='points1.inp'
          else
               WRITE(filename,*) filth2
               filename='points1.inp.'//trim(adjustl(filename))
          endif

          ! read in the coordinates
          inquire(iolength=i) (Q(j),j=1,3*Natoms)
          open(113,file=trim(filename),access='direct',form='unformatted',status='old',recl=i)
          read(113,rec=1) (Q(j),j=1,3*Natoms)
          close(113)

          ! check the coordinates look sensible
          if (MaxVal(Q)==0.0d0) then
              print *, 'Zero coordinates - stop'
              STOP
          endif
     end subroutine ReadInpFile

     subroutine WriteOutFile(Q,filename) 
          use commons,only: natoms
          implicit none
          double precision,dimension(:),intent(inout)::Q
          character(len=*),intent(in) :: filename

          integer :: i,j

          ! write out the coordinates
          inquire(iolength=i) (Q(j),j=1,3*Natoms)
          open(113,file=filename,access='direct',form='unformatted',status='unknown',recl=i)
          write(113,rec=1) (Q(j),j=1,3*Natoms)
          close(113)
     end subroutine WriteOutFile
END MODULE BinaryIO
