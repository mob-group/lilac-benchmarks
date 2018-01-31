
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
! jmc Note nothing has been done here to fix the unres pathlength coordinate resetting problem...
!
SUBROUTINE GOLDEN(VECS,QINIT,ETS,PINIT,POPT,PENERGY)
USE COMMONS, ONLY : NOPT, DEBUG
USE KEY, ONLY : PUSHOPTMAX, PUSHOPTCONV
IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: GRATIO=1.618033989D0
DOUBLE PRECISION ASTEP, BSTEP, CSTEP, DSTEP, RMS, VECS(NOPT), QINIT(NOPT), ENERGYC, ENERGYD
DOUBLE PRECISION QC(NOPT), QD(NOPT), VNEW(NOPT), ETS, POPT, PENERGY, PINIT
INTEGER J1
!
! Initially astep=0 the ts at QINIT(J1), bstep=PUSHOFF, corresponding to current Q(J1)
!
ASTEP=0.0D0
BSTEP=PINIT
CSTEP = BSTEP - (BSTEP - ASTEP) / GRATIO
DSTEP = ASTEP + (BSTEP - ASTEP) / GRATIO
IF (DEBUG) WRITE(6,'(A,4G20.10,A)') 'path> golden a, b, c, d=',ASTEP,BSTEP,CSTEP,DSTEP,' initially'
DO J1=1,PUSHOPTMAX
   QC(1:NOPT)=QINIT(1:NOPT)+CSTEP*VECS(1:NOPT)
   QD(1:NOPT)=QINIT(1:NOPT)+DSTEP*VECS(1:NOPT)
   CALL POTENTIAL(QC,ENERGYC,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   CALL POTENTIAL(QD,ENERGYD,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   IF (DEBUG) WRITE(6,'(A,4G25.15)') 'path> energy at ts, C and D and diff=', &
     &                    ETS,ENERGYC,ENERGYD,MIN(ENERGYC,ENERGYD)-ETS
   IF (ENERGYC.LT.ENERGYD) THEN
      BSTEP=DSTEP
   ELSE
      ASTEP=CSTEP
   ENDIF
   CSTEP = BSTEP - (BSTEP - ASTEP) / GRATIO
   DSTEP = ASTEP + (BSTEP - ASTEP) / GRATIO
!  IF (DEBUG) WRITE(6,'(A,4G20.10)') 'path> golden a, b, c, d=',ASTEP,BSTEP,CSTEP,DSTEP
   IF (ABS(ASTEP-BSTEP).LT.PUSHOPTCONV) THEN
      ASTEP=(ASTEP+BSTEP)/2.0D0
      EXIT
   ENDIF
ENDDO

PENERGY=MIN(ENERGYC,ENERGYD)
POPT=ASTEP
 
END SUBROUTINE GOLDEN


