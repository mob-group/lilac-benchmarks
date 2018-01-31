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
      SUBROUTINE SHIFTH(Q,UNDO,NOPT,NATOMS,ATMASS)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NOPT, NATOMS, NIADD
      DOUBLE PRECISION EV(3*NATOMS,6), RNORM(6), DUMMY, CMX, CMY, CMZ, TMASS,
     1                 AMASS(NATOMS), RMASS(NATOMS), Q(3*NATOMS), ATMASS(NATOMS),
     &                 IT(3,3), IV(3,3), OLDV(3*NATOMS), NEWV(3*NATOMS), QTEMP(3*NATOMS)

      LOGICAL UNDO
      INTEGER IPVT(6)
      DOUBLE PRECISION DET(2), ZWORK(6), RCOND

C     PRINT*,'in shifth shiftl=',shiftl
C
C  If MASST then the coordinates, Q, have already been multiplied by the square root of the mass
C
      SHIFTED=.TRUE.

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      TMASS=0.0D0
      NIADD=0
      IF (RIGIDBODY) THEN
         NATOMS=NATOMS/2
         NIADD=3*NATOMS
      ENDIF
C     PRINT*,'NATOMS,NIADD=',NATOMS,NIADD
      DO J2=1,NATOMS
         AMASS(J2)=1.0D0
         RMASS(J2)=1.0D0
         IF (MASST) AMASS(J2)=ATMASS(J2)
         IF (MASST) RMASS(J2)=SQRT(ATMASS(J2))
         TMASS=TMASS+AMASS(J2)
      ENDDO
      DO J2=1,NATOMS
         CMX=CMX+Q(3*(J2-1)+1)*RMASS(J2)
         CMY=CMY+Q(3*(J2-1)+2)*RMASS(J2)
         CMZ=CMZ+Q(3*(J2-1)+3)*RMASS(J2)
      ENDDO
      CMX=CMX/TMASS
      CMY=CMY/TMASS
      CMZ=CMZ/TMASS
C     PRINT*,'CMX,CMY,CMZ, NATOMS=',CMX,CMY,CMZ,NATOMS

      DO J1=1,6
         RNORM(J1)=0.0D0
         DO J2=1,3*NATOMS+NIADD
            EV(J2,J1)=0.0D0
         ENDDO
      ENDDO

      IF (MASST) THEN

         DO J1=1,NATOMS
            J3=3*(J1-1)
            DUMMY=SQRT(ATMASS(J1))

            EV(J3+1,1)=DUMMY
            RNORM(1)=RNORM(1)+DUMMY

            EV(J3+2,2)=DUMMY
            RNORM(2)=RNORM(2)+DUMMY

            EV(J3+3,3)=DUMMY
            RNORM(3)=RNORM(3)+DUMMY
   
            EV(J3+2,4)= Q(J3+3)*DUMMY
            EV(J3+3,4)=-Q(J3+2)*DUMMY
            RNORM(4)=RNORM(4)+(Q(J3+3)**2+Q(J3+2)**2)*DUMMY**2
     
            EV(J3+1,5)=-Q(J3+3)*DUMMY
            EV(J3+3,5)= Q(J3+1)*DUMMY
            RNORM(5)=RNORM(5)+(Q(J3+3)**2+Q(J3+1)**2)*DUMMY**2
    
            EV(J3+1,6)= Q(J3+2)*DUMMY
            EV(J3+2,6)=-Q(J3+1)*DUMMY
            RNORM(6)=RNORM(6)+(Q(J3+2)**2+Q(J3+1)**2)*DUMMY**2
         ENDDO

      ELSE

         DO J1=1,NATOMS
            J3=3*(J1-1)

            EV(J3+1,1)=1.0D0
            RNORM(1)=RNORM(1)+1.0D0

            EV(J3+2,2)=1.0D0
            RNORM(2)=RNORM(2)+1.0D0

            EV(J3+3,3)=1.0D0
            RNORM(3)=RNORM(3)+1.0D0
  
            EV(J3+2,4)= Q(J3+3)-CMZ
            EV(J3+3,4)=-Q(J3+2)+CMY
            RNORM(4)=RNORM(4)+(Q(J3+3)-CMZ)**2+(Q(J3+2)-CMY)**2

            EV(J3+1,5)=-Q(J3+3)+CMZ
            EV(J3+3,5)= Q(J3+1)-CMX
            RNORM(5)=RNORM(5)+(Q(J3+3)-CMZ)**2+(Q(J3+1)-CMX)**2

            EV(J3+1,6)= Q(J3+2)-CMY
            EV(J3+2,6)=-Q(J3+1)+CMX
            RNORM(6)=RNORM(6)+(Q(J3+2)-CMY)**2+(Q(J3+1)-CMX)**2
         ENDDO

      ENDIF

      DO J1=1,6
         RNORM(J1)=1.0D0/SQRT(RNORM(J1))
      ENDDO

      DO J1=1,NATOMS
         J3=3*(J1-1)

         EV(J3+1,1)=EV(J3+1,1)*RNORM(1)
         EV(J3+2,2)=EV(J3+2,2)*RNORM(2)
         EV(J3+3,3)=EV(J3+3,3)*RNORM(3)
 
         EV(J3+2,4)=EV(J3+2,4)*RNORM(4)
         EV(J3+3,4)=EV(J3+3,4)*RNORM(4)

         EV(J3+1,5)=EV(J3+1,5)*RNORM(5)
         EV(J3+3,5)=EV(J3+3,5)*RNORM(5)

         EV(J3+1,6)=EV(J3+1,6)*RNORM(6)
         EV(J3+2,6)=EV(J3+2,6)*RNORM(6)
      ENDDO

      IF (.TRUE.) THEN
C
C Build inertia tensor
C
         DO J2=1,NATOMS
            QTEMP(3*(J2-1)+1)=Q(3*(J2-1)+1)-CMX
            QTEMP(3*(J2-1)+2)=Q(3*(J2-1)+2)-CMY
            QTEMP(3*(J2-1)+3)=Q(3*(J2-1)+3)-CMZ
         ENDDO
         CALL INERTIA(IT,QTEMP)
C
C Diagonalize inertia tensor. The 0 flag reorders the e/values and
C e/vectors from smallest to largest.
C Second index of IV labels the eigenvectors.
C
         CALL EIG(IT,IV,3,3,0)

         DO J1=1,NATOMS
            DUMMY=1.0D0
            IF (MASST) DUMMY=SQRT(ATMASS(J1))
            EV(3*(J1-1)+1,1)=IV(1,1)*DUMMY
            EV(3*(J1-1)+2,1)=IV(2,1)*DUMMY
            EV(3*(J1-1)+3,1)=IV(3,1)*DUMMY
            EV(3*(J1-1)+1,2)=IV(1,2)*DUMMY
            EV(3*(J1-1)+2,2)=IV(2,2)*DUMMY
            EV(3*(J1-1)+3,2)=IV(3,2)*DUMMY
            EV(3*(J1-1)+1,3)=IV(1,3)*DUMMY
            EV(3*(J1-1)+2,3)=IV(2,3)*DUMMY
            EV(3*(J1-1)+3,3)=IV(3,3)*DUMMY
            EV(3*(J1-1)+1,4)=(IV(3,1)*QTEMP(3*(J1-1)+2)-IV(2,1)*QTEMP(3*(J1-1)+3))*DUMMY
            EV(3*(J1-1)+2,4)=(IV(1,1)*QTEMP(3*(J1-1)+3)-IV(3,1)*QTEMP(3*(J1-1)+1))*DUMMY
            EV(3*(J1-1)+3,4)=(IV(2,1)*QTEMP(3*(J1-1)+1)-IV(1,1)*QTEMP(3*(J1-1)+2))*DUMMY
            EV(3*(J1-1)+1,5)=(IV(3,2)*QTEMP(3*(J1-1)+2)-IV(2,2)*QTEMP(3*(J1-1)+3))*DUMMY
            EV(3*(J1-1)+2,5)=(IV(1,2)*QTEMP(3*(J1-1)+3)-IV(3,2)*QTEMP(3*(J1-1)+1))*DUMMY
            EV(3*(J1-1)+3,5)=(IV(2,2)*QTEMP(3*(J1-1)+1)-IV(1,2)*QTEMP(3*(J1-1)+2))*DUMMY
            EV(3*(J1-1)+1,6)=(IV(3,3)*QTEMP(3*(J1-1)+2)-IV(2,3)*QTEMP(3*(J1-1)+3))*DUMMY
            EV(3*(J1-1)+2,6)=(IV(1,3)*QTEMP(3*(J1-1)+3)-IV(3,3)*QTEMP(3*(J1-1)+1))*DUMMY
            EV(3*(J1-1)+3,6)=(IV(2,3)*QTEMP(3*(J1-1)+1)-IV(1,3)*QTEMP(3*(J1-1)+2))*DUMMY
         ENDDO
         DO J1=1,6
            RNORM(J1)=0.0D0
            DO J2=1,3*NATOMS
               RNORM(J1)=RNORM(J1)+EV(J2,J1)**2
            ENDDO
         ENDDO
         DO J1=1,6
            RNORM(J1)=1.0D0/SQRT(RNORM(J1))
            DO J2=1,3*NATOMS
               EV(J2,J1)=EV(J2,J1)*RNORM(J1)
            ENDDO
         ENDDO
      ENDIF

C     PRINT*,'unshifted hess:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)

      IF (FREEZE) THEN
         DO J1=1,NATOMS+NIADD/3
            IF (.NOT.FROZEN(J1))  CYCLE
            J3=3*(J1-1)
            HESS(J3+1,J3+1)=HESS(J3+1,J3+1)+SHIFTL(1)
            HESS(J3+2,J3+2)=HESS(J3+2,J3+2)+SHIFTL(1)
            HESS(J3+3,J3+3)=HESS(J3+3,J3+3)+SHIFTL(1)
         ENDDO
      ELSEIF (.NOT.UNDO) THEN
!        DO J1=1,NATOMS+NIADD/3
!           J3=3*(J1-1)
!           DO J2=1,NATOMS+NIADD/3
!              J4=3*(J2-1)
!  
!              HESS(J4+1,J3+1)=HESS(J4+1,J3+1)-SHIFTL(1)*EV(J4+1,1)*EV(J3+1,1)
!    1                                        -SHIFTL(5)*EV(J4+1,5)*EV(J3+1,5)-SHIFTL(6)*EV(J4+1,6)*EV(J3+1,6)
!              HESS(J4+1,J3+2)=HESS(J4+1,J3+2)-SHIFTL(6)*EV(J4+1,6)*EV(J3+2,6)
!              HESS(J4+1,J3+3)=HESS(J4+1,J3+3)-SHIFTL(5)*EV(J4+1,5)*EV(J3+3,5)
!              HESS(J4+2,J3+1)=HESS(J4+2,J3+1)-SHIFTL(6)*EV(J4+2,6)*EV(J3+1,6)
!              HESS(J4+2,J3+2)=HESS(J4+2,J3+2)-SHIFTL(2)*EV(J4+2,2)*EV(J3+2,2)
!    1                                        -SHIFTL(4)*EV(J4+2,4)*EV(J3+2,4)-SHIFTL(6)*EV(J4+2,6)*EV(J3+2,6)
!              HESS(J4+2,J3+3)=HESS(J4+2,J3+3)-SHIFTL(4)*EV(J4+2,4)*EV(J3+3,4)
!              HESS(J4+3,J3+1)=HESS(J4+3,J3+1)-SHIFTL(5)*EV(J4+3,5)*EV(J3+1,5)
!              HESS(J4+3,J3+2)=HESS(J4+3,J3+2)-SHIFTL(4)*EV(J4+3,4)*EV(J3+2,4)
!              HESS(J4+3,J3+3)=HESS(J4+3,J3+3)-SHIFTL(3)*EV(J4+3,3)*EV(J3+3,3)
!    1                                        -SHIFTL(4)*EV(J4+3,4)*EV(J3+3,4)-SHIFTL(5)*EV(J4+3,5)*EV(J3+3,5)
!
!           ENDDO
!        ENDDO
         DO J1=1,3*NATOMS+NIADD
            DO J2=1,3*NATOMS+NIADD
               DO J3=1,6
                  HESS(J1,J2)=HESS(J1,J2)-SHIFTL(J3)*EV(J1,J3)*EV(J2,J3)
               ENDDO
            ENDDO
         ENDDO
      ELSE
!        DO J1=1,NATOMS+NIADD/3
!           J3=3*(J1-1)
!           DO J2=1,NATOMS+NIADD/3
!              J4=3*(J2-1)
!              HESS(J4+1,J3+1)=HESS(J4+1,J3+1)+SHIFTL(1)*EV(J4+1,1)*EV(J3+1,1)
!    1                                        +SHIFTL(5)*EV(J4+1,5)*EV(J3+1,5)+SHIFTL(6)*EV(J4+1,6)*EV(J3+1,6)
!              HESS(J4+1,J3+2)=HESS(J4+1,J3+2)+SHIFTL(6)*EV(J4+1,6)*EV(J3+2,6)
!              HESS(J4+1,J3+3)=HESS(J4+1,J3+3)+SHIFTL(5)*EV(J4+1,5)*EV(J3+3,5)
!              HESS(J4+2,J3+1)=HESS(J4+2,J3+1)+SHIFTL(6)*EV(J4+2,6)*EV(J3+1,6)
!              HESS(J4+2,J3+2)=HESS(J4+2,J3+2)+SHIFTL(2)*EV(J4+2,2)*EV(J3+2,2)
!    1                                        +SHIFTL(4)*EV(J4+2,4)*EV(J3+2,4)+SHIFTL(6)*EV(J4+2,6)*EV(J3+2,6)
!              HESS(J4+2,J3+3)=HESS(J4+2,J3+3)+SHIFTL(4)*EV(J4+2,4)*EV(J3+3,4)
!              HESS(J4+3,J3+1)=HESS(J4+3,J3+1)+SHIFTL(5)*EV(J4+3,5)*EV(J3+1,5)
!              HESS(J4+3,J3+2)=HESS(J4+3,J3+2)+SHIFTL(4)*EV(J4+3,4)*EV(J3+2,4)
!              HESS(J4+3,J3+3)=HESS(J4+3,J3+3)+SHIFTL(3)*EV(J4+3,3)*EV(J3+3,3)
!    1                                        +SHIFTL(4)*EV(J4+3,4)*EV(J3+3,4)+SHIFTL(5)*EV(J4+3,5)*EV(J3+3,5)
!
!           ENDDO
!        ENDDO
         DO J1=1,3*NATOMS+NIADD
            DO J2=1,3*NATOMS+NIADD
               DO J3=1,6
                  HESS(J1,J2)=HESS(J1,J2)+SHIFTL(J3)*EV(J1,J3)*EV(J2,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (TWOD) THEN
         DO J2=3,3*NATOMS,3
            HESS(J2,J2)=SHIFTV
         ENDDO
      ENDIF

C     PRINT*,'shifted hess:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)
C     PRINT*,'coordinates:'
C     WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

      IF (RIGIDBODY) NATOMS=NATOMS*2

      RETURN
      END
