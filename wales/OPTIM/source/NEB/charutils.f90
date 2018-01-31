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
MODULE CHARUTILS
     IMPLICIT NONE
     SAVE
     INTEGER,PARAMETER :: INTSTRLENGTH=10,REALSTRLENGTH=20
     CHARACTER(LEN=INTSTRLENGTH) :: INTSTR,INTSTR2
     CHARACTER(LEN=REALSTRLENGTH) :: REALSTR,REALSTR2
     CONTAINS

     FUNCTION WI(I)
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I
          CHARACTER(LEN=INTSTRLENGTH) :: WI
          
          CHARACTER(LEN=INTSTRLENGTH) :: STR

          WRITE(STR,'(i10)') i
          WI = ADJUSTL(STR)
     END FUNCTION WI

     FUNCTION WR(R,D)
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: D !NUMBER OF DIGITS AFTER POINT TO WRITE; RANGE 1-9.
          DOUBLE PRECISION,INTENT(IN) :: R

          INTEGER :: I
          CHARACTER(LEN=REALSTRLENGTH) :: WR

          CHARACTER(LEN=REALSTRLENGTH) :: STR,FORMT
          WRITE(FORMT,'(a5,i1,a1)') '(f20.',d,')'
          WRITE(STR,FORMT) R
          WR = ADJUSTL(STR)
          IF (D==0) THEN
               DO I=20,1,-1
                    IF (WR(I:I)==".") then
                         WR(I:I)=" "; exit
                    ENDIF
               ENDDO
          ENDIF
     END FUNCTION WR

     FUNCTION RM0S(STR)
          IMPLICIT NONE

          CHARACTER(LEN=REALSTRLENGTH) :: STR,RM0S

          INTEGER :: I

          RM0S=ADJUSTR(STR)

          DO I=REALSTRLENGTH,2,-1
               IF (.NOT.RM0S(I:I)=="0") then
                    EXIT
               ENDIF
          ENDDO
          
          RM0S=ADJUSTL(RM0S(1:I))
     END FUNCTION RM0S
END MODULE CHARUTILS
