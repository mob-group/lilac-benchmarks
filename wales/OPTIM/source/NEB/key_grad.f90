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
MODULE KEYGRAD
     IMPLICIT NONE
     SAVE
     LOGICAL :: ORT               = .FALSE.
     CHARACTER(LEN=5) :: GRADTYPE = "dneb"
     CONTAINS

     SUBROUTINE KEYGRADPRINT
          USE CHARUTILS
          IMPLICIT NONE
          IF (ORT) THEN
               WRITE(*,'(1x,a)') 'KeyGrad> Overall rotation and translation will be removed'
          ELSE
               WRITE(*,'(1x,a)') 'KeyGrad> Overall rotation and translation will NOT be removed'
          ENDIF
          
          SELECT CASE (GRADTYPE)
               CASE ("jnew")
                    WRITE(*,'(1x,a)') 'KeyGrad> Gradient will be constructed according to new recommendations of Jonsson'
               CASE ("jold")
                    WRITE(*,'(1x,a)') 'KeyGrad> Gradient will be constructed according to old recommendations of Jonsson'
               CASE ("dneb")
                    WRITE(*,'(1x,a)') 'KeyGrad> Using doubly nudged elastic band gradient'
          END SELECT
     END SUBROUTINE KEYGRADPRINT
END MODULE KEYGRAD
