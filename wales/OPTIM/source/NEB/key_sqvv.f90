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
MODULE KEYSQVV
     IMPLICIT NONE
     SAVE

     DOUBLE PRECISION :: DT                = 0.01D0
     
     LOGICAL :: LIMITSQVVSTEP     = .FALSE.
     DOUBLE PRECISION :: STEPDOFMAX        = 0.01D0

     LOGICAL :: SQVVGUESS         = .FALSE.
     INTEGER :: NITERSQVVGUESSMAX = 300
     DOUBLE PRECISION :: SQVVGUESSRMSTOL   = 2.0D0

     LOGICAL :: READVEL           = .FALSE.
     CONTAINS

     SUBROUTINE KEYSQVVPRINT
          USE CHARUTILS
          IMPLICIT NONE
          REALSTR=RM0S(WR(DT,9))          
          WRITE(*,'(1x,a)') 'KeySQVV> Using time integration step of '//trim(RealStr)

          IF (LIMITSQVVSTEP) THEN
               REALSTR=RM0S(WR(STEPDOFMAX,9))
               WRITE(*,'(1x,a)') 'KeySQVV> Maximal step size for each degree of freedom = '//trim(RealStr)
          ELSE
               WRITE(*,'(1x,a)') 'KeySQVV> Unlimited step'
          ENDIF

          IF (SQVVGUESS) THEN
               INTSTR=WI(NITERSQVVGUESSMAX)
               WRITE(*,'(1x,a)') 'KeySQVV> Up to '//trim(IntStr)//' steps of SQVV initial relaxation are allowed'
               REALSTR=RM0S(WR(SQVVGUESSRMSTOL,9))
               WRITE(*,'(1x,a)') 'KeySQVV> RMS convergence tolerance is '//trim(RealStr)
          ELSE
               WRITE(*,'(1x,a)') 'KeySQVV> No preoptimization was requested'
          ENDIF

          IF (READVEL) THEN
               WRITE(*,'(1x,a)') 'KeySQVV> Reading in velocity file'
          ENDIF
     END SUBROUTINE KEYSQVVPRINT
END MODULE KEYSQVV
