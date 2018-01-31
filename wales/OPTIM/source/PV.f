C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C  Uses two calls to potential per variable box length to make numerical first
C  and second derivatives for box length optimisation by eigenvector-following.
C  Zero second derivative results if the energy is linear in the box length if
C  fractional coordinates are not used.
C  This can occur if FIXIMAGE is turned off for a line minimisation and LBFGS
C  takes a daft step, giving silly ANV, which does not get fixed.
C
C  Analytic derivatives for normal and fractional coordinates added for 'LS' 
C  atom type 26/6/01.
C
      SUBROUTINE PVOPT(COORDS,EZERO,VNEW)
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      LOGICAL PVAGAIN, PVFLAG
      DOUBLE PRECISION EZERO, COORDS(3*NATOMS), DIFF, GRAD, SECOND, SMAX, VNEW(3*NATOMS), 
     1                 EPLUS, EMINUS, RMS, EXP, EXM, EYP, EYM, EZP, EZM, GRADX, GRADY, GRADZ, SECONDX, SECONDY, SECONDZ,
     2                 PPARAM1, PPARAM2, PPARAM3, STEPX, STEPY, STEPZ, EPREV, EEST
      INTEGER NCOUNT, J1
      LOGICAL CUBIC, NRLOCAL, ANALYTIC
      COMMON /PVF/ PVFLAG
      COMMON /CUB/ CUBIC
      COMMON /BOXDERIVS/ GRADX, GRADY, GRADZ, SECONDX, SECONDY, SECONDZ
      
      ANALYTIC=.FALSE.
      IF (ZSYM(NATOMS).EQ.'LS') ANALYTIC=.TRUE.
      DIFF=0.001D0  
      SMAX=0.05D0
      TRAD=1.0D0
      PVFLAG=.FALSE.
C     FIXIMAGE=.TRUE.
C     IF (NORESET) THEN
C        NRLOCAL=.TRUE.
C        NORESET=.FALSE.
C     ENDIF
      NCOUNT=0
10    PVAGAIN=.FALSE.

      IF ((DABS(PARAM1-PARAM2).LT.1.0D-10).AND.(DABS(PARAM3-PARAM2).LT.1.0D-10).AND.CUBIC) THEN
         IF (ANALYTIC) THEN
            CALL POTENTIAL(COORDS,EZERO,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.TRUE.)
            GRAD=GRADX
            SECOND=SECONDX
         ELSE
            PARAM1=PARAM1+DIFF
            PARAM2=PARAM2+DIFF
            PARAM3=PARAM3+DIFF
            CALL POTENTIAL(COORDS,EPLUS,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM1=PARAM1-2*DIFF
            PARAM2=PARAM2-2*DIFF
            PARAM3=PARAM3-2*DIFF
            CALL POTENTIAL(COORDS,EMINUS,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM1=PARAM1+DIFF
            PARAM2=PARAM2+DIFF
            PARAM3=PARAM3+DIFF

            GRAD=(EPLUS-EMINUS)/(2.0D0*DIFF) 
      
            SECOND=(EPLUS+EMINUS-2.0D0*EZERO)/(DIFF*DIFF) 
         ENDIF

         WRITE(*,'(A)') ' '
         WRITE(*,'(A,T50,F20.10)') 'Enthalpy in lattice optimisation is ',EZERO
         WRITE(*,'(A,T50,F20.10)') 'Gradient wrt box length=',GRAD
         WRITE(*,'(A,T50,F20.10)') 'Second derivative wrt box length=',SECOND
         WRITE(*,'(A,T50,F20.10)') 'Full step size=',-GRAD/SECOND
         IF (SECOND.LT.0.0D0) SECOND=-SECOND
         PPARAM1=PARAM1
         PPARAM2=PARAM2
         PPARAM3=PARAM3
         IF (DABS(GRAD/SECOND).GT.SMAX) THEN
            PARAM1=PARAM1-SMAX*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
            WRITE(*,'(A,T50,F20.10,A,F15.10)') 'Scaled step=',-SMAX*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD)),' new box length=',PARAM1
         ELSE
            PARAM1=PARAM1-GRAD/SECOND
            WRITE(*,'(A,T50,F20.10,A,F15.10)') 'Box length step=',-GRAD/SECOND,' new box length=',PARAM1
         ENDIF
         PARAM2=PARAM1
         PARAM3=PARAM1
         IF (ABS(GRAD).LT.PVCONV) PVFLAG=.TRUE.
         IF (ABS(GRAD).GT.PVTOL) PVAGAIN=.TRUE.
      ELSE
         IF (NCOUNT.GT.0) THEN
            EEST=EZERO+GRADX*STEPX+GRADY*STEPY+GRADZ*STEPZ
     1                   +(SECONDX*STEPX**2+SECONDY*STEPY**2+SECONDZ*STEPZ**2)/2.0D0
            EPREV=EZERO
         ENDIF
         IF (ANALYTIC) THEN
            CALL POTENTIAL(COORDS,EZERO,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.TRUE.)
         ELSE
            PARAM1=PARAM1+DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+1)=COORDS(3*(J1-1)+1)*PARAM1/(PARAM1-DIFF)
               ENDDO
            ENDIF
            CALL POTENTIAL(COORDS,EXP,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM1=PARAM1-2*DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+1)=COORDS(3*(J1-1)+1)*PARAM1/(PARAM1+2*DIFF)
               ENDDO
            ENDIF
            CALL POTENTIAL(COORDS,EXM,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM1=PARAM1+DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+1)=COORDS(3*(J1-1)+1)*PARAM1/(PARAM1-DIFF)
               ENDDO
            ENDIF
   
            PARAM2=PARAM2+DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+2)=COORDS(3*(J1-1)+2)*PARAM2/(PARAM2-DIFF)
               ENDDO
            ENDIF
            CALL POTENTIAL(COORDS,EYP,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM2=PARAM2-2*DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+2)=COORDS(3*(J1-1)+2)*PARAM2/(PARAM2+2*DIFF)
               ENDDO
            ENDIF
            CALL POTENTIAL(COORDS,EYM,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM2=PARAM2+DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+2)=COORDS(3*(J1-1)+2)*PARAM2/(PARAM2-DIFF)
               ENDDO
            ENDIF
   
            PARAM3=PARAM3+DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+3)=COORDS(3*(J1-1)+3)*PARAM3/(PARAM3-DIFF)
               ENDDO
            ENDIF
            CALL POTENTIAL(COORDS,EZP,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM3=PARAM3-2*DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+3)=COORDS(3*(J1-1)+3)*PARAM3/(PARAM3+2*DIFF)
               ENDDO
            ENDIF
            CALL POTENTIAL(COORDS,EZM,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PARAM3=PARAM3+DIFF
            IF (FRACTIONAL) THEN
               DO J1=1,NATOMS
                  COORDS(3*(J1-1)+3)=COORDS(3*(J1-1)+3)*PARAM3/(PARAM3-DIFF)
               ENDDO
            ENDIF
   
            GRADX=(EXP-EXM)/(2.0D0*DIFF)
            SECONDX=(EXP+EXM-2.0D0*EZERO)/(DIFF*DIFF) 
            GRADY=(EYP-EYM)/(2.0D0*DIFF)
            SECONDY=(EYP+EYM-2.0D0*EZERO)/(DIFF*DIFF) 
            GRADZ=(EZP-EZM)/(2.0D0*DIFF)
            SECONDZ=(EZP+EZM-2.0D0*EZERO)/(DIFF*DIFF) 
         ENDIF

C        IF (NCOUNT.GT.0) THEN
C           IF (ABS((EZERO-EEST)/(EZERO-EPREV)).LT.TRAD) THEN
C              SMAX=SMAX*1.05D0
C           ELSE
C              SMAX=SMAX/1.1D0
C           ENDIF
C           WRITE(*,'(A,4F20.10)') 'delta E actual and predicted, ratio, SMAX: ',EZERO-EPREV,EEST-EPREV,
C    1                             ABS((EZERO-EEST)/(EZERO-EPREV)),SMAX
C        ENDIF

         WRITE(*,'(A)') ' '
         WRITE(*,'(A,T50,F20.10)') 'Enthalpy in lattice optimisation is ',EZERO
         WRITE(*,'(A,T50,3F20.10)') 'Gradient wrt x,y,z box lengths=',GRADX,GRADY,GRADZ
         WRITE(*,'(A,T50,3F20.10)') 'Second derivatives wrt x,y,z box lengths=',SECONDX,SECONDY,SECONDZ
         IF ((.NOT.PVTS).AND.((MIN(MIN(SECONDX,SECONDY),SECONDZ).LT.MAX(MAX(SECONDX,SECONDY),SECONDZ)/10.0D0).OR.
     1        (MAX(MAX(SECONDX,SECONDY),SECONDZ).LT.1.0D0))) THEN
            PRINT*,'Small box length curvature detected - skipping box length optimisation'
            RETURN
         ENDIF

         IF (SECONDX.EQ.0.0D0) SECONDX=1.0D0
         IF (SECONDY.EQ.0.0D0) SECONDY=1.0D0
         IF (SECONDZ.EQ.0.0D0) SECONDZ=1.0D0
C
C  Fractional coordinates - steps seem to need damping
C    
         IF (FRACTIONAL) THEN
            SECONDX=SECONDX*1.2D0
            SECONDY=SECONDY*1.2D0
            SECONDZ=SECONDZ*1.2D0
         ENDIF

         STEPX=-2.0D0*GRADX/(DABS(SECONDX)*(1.0D0+DSQRT(1.0D0 + 4.0D0*(GRADX/SECONDX)**2)))
         STEPY=-2.0D0*GRADY/(DABS(SECONDY)*(1.0D0+DSQRT(1.0D0 + 4.0D0*(GRADY/SECONDY)**2)))
         STEPZ=-2.0D0*GRADZ/(DABS(SECONDZ)*(1.0D0+DSQRT(1.0D0 + 4.0D0*(GRADZ/SECONDZ)**2)))
         WRITE(*,'(A,T50,3F20.10)') 'Full   step size for x,y,z =',STEPX,STEPY,STEPZ

         IF (PVTS) THEN
            IF (NBOXTS.EQ.1) THEN
               STEPX=-STEPX
               IF (SECONDX.GT.0.0D0) STEPX=PUSHOFF*ABS(STEPX)/STEPX
            ELSE IF (NBOXTS.EQ.2) THEN
               STEPY=-STEPY
               IF (SECONDY.GT.0.0D0) STEPY=PUSHOFF*ABS(STEPY)/STEPY
            ELSE 
               STEPZ=-STEPZ
               IF (SECONDZ.GT.0.0D0) STEPZ=PUSHOFF*ABS(STEPZ)/STEPZ
            ENDIF
         ENDIF

         PPARAM1=PARAM1
         IF (DABS(STEPX).GT.SMAX) THEN
            PARAM1=PARAM1+SMAX*DABS(STEPX)/STEPX
         ELSE
            PARAM1=PARAM1+STEPX
         ENDIF

         PPARAM2=PARAM2
         IF (DABS(STEPY).GT.SMAX) THEN
            PARAM2=PARAM2+SMAX*DABS(STEPY)/STEPY
         ELSE
            PARAM2=PARAM2+STEPY
         ENDIF

         PPARAM3=PARAM3
         IF (DABS(STEPZ).GT.SMAX) THEN
            PARAM3=PARAM3+SMAX*DABS(STEPZ)/STEPZ
         ELSE
            PARAM3=PARAM3+STEPZ
         ENDIF
         STEPX=PARAM1-PPARAM1
         STEPY=PARAM2-PPARAM2
         STEPZ=PARAM3-PPARAM3
         WRITE(*,'(A,T50,3F20.10)') 'Actual step size for x,y,z =',PARAM1-PPARAM1,PARAM2-PPARAM2,PARAM3-PPARAM3
         WRITE(*,'(A,T50,3F20.10)') 'New x,y,z box lengths:',PARAM1,PARAM2,PARAM3
         IF ((ABS(GRADX).LT.PVCONV).AND.(ABS(GRADY).LT.PVCONV).AND.(ABS(GRADZ).LT.PVCONV)) PVFLAG=.TRUE.
         IF ((ABS(GRADX).GT.PVTOL).OR.(ABS(GRADY).GT.PVTOL).OR.(ABS(GRADZ).GT.PVTOL)) PVAGAIN=.TRUE.
      ENDIF

      IF (FRACTIONAL) THEN
         DO J1=1,NATOMS
            COORDS(3*(J1-1)+1)=COORDS(3*(J1-1)+1)*PARAM1/PPARAM1
            COORDS(3*(J1-1)+2)=COORDS(3*(J1-1)+2)*PARAM2/PPARAM2
            COORDS(3*(J1-1)+3)=COORDS(3*(J1-1)+3)*PARAM3/PPARAM3
         ENDDO
      ENDIF

      IF (PVAGAIN) THEN
C        FIXIMAGE=.FALSE.
         PRINT*,'WARNING - NRLOCAL was unset - not tested!'
C        IF (NRLOCAL) NORESET=.TRUE.
         IF (.NOT.ANALYTIC) CALL POTENTIAL(COORDS,EZERO,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
C        FIXIMAGE=.TRUE.
C        IF (NRLOCAL) NORESET=.FALSE.
         NCOUNT=NCOUNT+1
         IF (NCOUNT.LT.PVSTEPS) GOTO 10
         WRITE(*,'(A)') 'Warning, maximum number of box length steps exceeded in PVopt'
      ENDIF
C     FIXIMAGE=.FALSE.
C     IF (NRLOCAL) NORESET=.TRUE.

      RETURN
      END
