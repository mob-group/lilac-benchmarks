!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
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
MODULE KEYTAU
     IMPLICIT NONE
     SAVE
     INTEGER :: TANTYPE=1
     CONTAINS

     SUBROUTINE KEYTAUPRINT
          IMPLICIT NONE

          SELECT CASE (TANTYPE)
          CASE (1)
               WRITE(*,'(1x,a)') "KeyTau> Using Henkelman and Jonsson's improved tangent"
          CASE (2)
               WRITE(*,'(1x,a)') 'KeyTau> Using normalized line segment'
          CASE (4)
               WRITE(*,'(1x,a)') 'KeyTau> UNRES tangent (in internals?) <- please correct this print statement'
          END SELECT
     END SUBROUTINE KEYTAUPRINT
END MODULE KEYTAU
