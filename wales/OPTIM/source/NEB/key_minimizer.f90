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
MODULE KEYMINIMIZER
     IMPLICIT NONE
     SAVE
     INTEGER :: NITERMIN         = 0
     INTEGER :: NITERMAX         = 1
     DOUBLE PRECISION :: RMSTOL           = 0.01D0
     CHARACTER(LEN=5) :: MINTYPE = "lbfgs"
     CONTAINS
     
     SUBROUTINE KEYMINIMIZERPRINT(VARIABLEIN)
          USE CHARUTILS
          USE KEYLBFGS,ONLY:KEYLBFGSPRINT
          USE KEYSQVV,ONLY:KEYSQVVPRINT
          IMPLICIT NONE
          LOGICAL,INTENT(IN),OPTIONAL :: VARIABLEIN

          LOGICAL :: VARIABLE

          VARIABLE=.FALSE.
          IF (PRESENT(VARIABLEIN)) VARIABLE=VARIABLEIN

          IF (VARIABLE) THEN
               WRITE(*,'(1x,a)') 'KeyMin> Maximal number of iterations will vary,'//&
               &' depending on the number of images in the band'
          ELSE
               INTSTR=WI(NITERMIN)
               INTSTR2=WI(NITERMAX)
               WRITE(*,'(1x,a)') 'KeyMin> Number of iterations: min, max = '//trim(IntStr)//', '//trim(IntStr2)
          ENDIF

          REALSTR=RM0S(WR(RMSTOL,9))
          WRITE(*,'(1x,a)') 'KeyMin> RMS convergence criterion is set to '//trim(RealStr)

          SELECT CASE (MINTYPE)
          CASE ("lbfgs")
               WRITE(*,'(1x,a)') 'KeyMin> L-BFGS minimization'
               CALL KEYLBFGSPRINT
          CASE ("sqvv")
               WRITE(*,'(1x,a)') 'KeyMin> SQVV minimisation'
               CALL KEYSQVVPRINT
          END SELECT
     END SUBROUTINE KEYMINIMIZERPRINT
END MODULE KEYMINIMIZER
