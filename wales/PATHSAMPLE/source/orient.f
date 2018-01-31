!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Put permutational isomers into a standard orientation.
C
      SUBROUTINE ORIENT(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION Q1(3*NATOMS), DIST(NATOMS), DMAX, RVEC(3), T1(3*NATOMS), CMX, CMY, CMZ,
     1                 COST, SINT, RDOTN, DMAX2, PROJ
      DOUBLE PRECISION, PARAMETER :: ORBDISTTOL=1.0D-2 ! used to be 0.001 - needs to be smaller for C60 to work
      INTEGER J1, I, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
C
C  Move centre of mass to the origin.
C
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+Q1(3*(I-1)+1)
         CMY=CMY+Q1(3*(I-1)+2)
         CMZ=CMZ+Q1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
C     PRINT*,'CMX,CMY,CMZ=',CMX,CMY,CMZ
      DO I=1,NATOMS
         Q1(3*(I-1)+1)=Q1(3*(I-1)+1)-CMX
         Q1(3*(I-1)+2)=Q1(3*(I-1)+2)-CMY
         Q1(3*(I-1)+3)=Q1(3*(I-1)+3)-CMZ
      ENDDO
      
      DMAX=-1.0D0
      NORBIT1=1
      DO J1=1,NATOMS
         DIST(J1)=SQRT(Q1(3*(J1-1)+1)**2+Q1(3*(J1-1)+2)**2+Q1(3*(J1-1)+3)**2)
         IF (ABS(DIST(J1)-DMAX).LT.ORBDISTTOL) THEN
            NORBIT1=NORBIT1+1
            IF (NORBIT1.EQ.NCHOOSE1) THEN
               JMAX1=J1
            ENDIF
         ELSE IF (DIST(J1).GT.DMAX) THEN
            DMAX=DIST(J1)
            NORBIT1=1
            JMAX1=J1
         ENDIF
      ENDDO
C
C  For tagged atoms the choice of the first atom matters if it belongs to an orbit of size > 1.
C
!     IF (DEBUG) PRINT*,'atom ',JMAX1,' will be moved onto the z axis, distance=',DMAX
      IF ((ABS(Q1(3*(JMAX1-1)+1)).LT.1.0D-8).AND.(ABS(Q1(3*(JMAX1-1)+2)).LT.1.0D-8)) THEN
         IF (Q1(3*(JMAX1-1)+3).GT.0.0D0) THEN
            T1(1:3*NATOMS)=Q1(1:3*NATOMS)
         ELSE  ! rotate about the x axis DO NOT INVERT!!
            DO J1=1,NATOMS
               T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)
               T1(3*(J1-1)+2)=-Q1(3*(J1-1)+2)
               T1(3*(J1-1)+3)=-Q1(3*(J1-1)+3)
            ENDDO
         ENDIF
         GOTO 10
      ENDIF
      COST=Q1(3*(JMAX1-1)+3)/DMAX
      SINT=SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)/DMAX
C
C  Rotate atom JMAX1 onto the z axis.
C  Rotate all the atoms through ANGLE about RVEC. Use rotation formula
C  from Goldstein p. 165.
C
      RVEC(1)= Q1(3*(JMAX1-1)+2)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
      RVEC(2)=-Q1(3*(JMAX1-1)+1)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
      RVEC(3)=0.0D0
      DO J1=1,NATOMS
C        IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=Q1(3*(J1-1)+1)*RVEC(1)+Q1(3*(J1-1)+2)*RVEC(2)+Q1(3*(J1-1)+3)*RVEC(3)
            T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)-(Q1(3*(J1-1)+2)*RVEC(3)-Q1(3*(J1-1)+3)*RVEC(2))*SINT
            T1(3*(J1-1)+2)=Q1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)-(Q1(3*(J1-1)+3)*RVEC(1)-Q1(3*(J1-1)+1)*RVEC(3))*SINT
            T1(3*(J1-1)+3)=Q1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)-(Q1(3*(J1-1)+1)*RVEC(2)-Q1(3*(J1-1)+2)*RVEC(1))*SINT
C        ENDIF
      ENDDO

C     DO J1=1,3*NATOMS
C        Q1(J1)=T1(J1)
C     ENDDO

10    CONTINUE
C     IF (DEBUG) THEN
C        WRITE(4,'(I5)') NATOMS
C        WRITE(4,'(A)') 'after first orient rotation:'
C        WRITE(4,'(A5,3F20.10)') 'LA   ',T1(1),T1(2),T1(3)
C        DO J1=2,NATOMS
C           WRITE(4,'(A5,3F20.10)') 'LB   ',T1(3*(J1-1)+1),T1(3*(J1-1)+2),T1(3*(J1-1)+3)
C        ENDDO
C     ENDIF
C
C  Find the atom with the largest distance from the z axis.
C
      DMAX=-1.0D0
      DO J1=1,NATOMS
         DIST(J1)=SQRT(T1(3*(J1-1)+1)**2+T1(3*(J1-1)+2)**2)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
      ENDDO
      DMAX2=-1.0D100
      DO J1=1,NATOMS
         IF (ABS(DIST(J1)-DMAX).LT.ORBDISTTOL) THEN
            CALL ROTXZ(NATOMS,J1,T1,PROJ,DIST)
            IF (ABS(PROJ-DMAX2).LT.ORBDISTTOL) THEN
               NORBIT2=NORBIT2+1
               IF (NORBIT2.EQ.NCHOOSE2) THEN
                  JMAX2=J1
               ENDIF
            ELSE IF (PROJ.GT.DMAX2) THEN
               NORBIT2=1
               DMAX2=PROJ
               JMAX2=J1
            ENDIF
         ENDIF
      ENDDO
!     IF (DEBUG) PRINT*,'NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX1,JMAX2=',NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX1,JMAX2
C
C  and rotate it into the xz plane.
C
      CALL ROTXZ(NATOMS,JMAX2,T1,DMAX2,DIST)

C     IF (DEBUG) PRINT*,'after'
C     IF (DEBUG) WRITE(*,'(3F20.10)') (T1(J1),J1=1,3*NATOMS)

      RETURN
      END

      SUBROUTINE ROTXZ(NATOMS,JDO,T1,PROJ,DIST)
      IMPLICIT NONE
      INTEGER NATOMS, JDO, J1
      DOUBLE PRECISION T1(3*NATOMS), PROJ, DIST(NATOMS), RVEC(3), COST, SINT, RDOTN, TX, TY, TZ

C     PRINT '(I6)',NATOMS
C     PRINT '(A,I6,G20.10)','before rotation'
C     WRITE(*,'(A2,2X,3G20.10)') ('LA',T1(3*(J3-1)+1),T1(3*(J3-1)+2),T1(3*(J3-1)+3),J3=1,NATOMS)

      IF (ABS(T1(3*(JDO-1)+2)).LT.1.0D-8) THEN
         IF (T1(3*(JDO-1)+1).LT.0.0D0) THEN ! rotate about the z axis DO NOT INVERT!!
            DO J1=1,NATOMS
               T1(3*(J1-1)+1)=-T1(3*(J1-1)+1)
               T1(3*(J1-1)+2)=-T1(3*(J1-1)+2)
            ENDDO
         ENDIF
         GOTO 20
      ENDIF

      COST=T1(3*(JDO-1)+1)/DIST(JDO)
      SINT=T1(3*(JDO-1)+2)/DIST(JDO)
C     PRINT '(A,4G20.10)','T1(3*(JDO-1)+2),COST,SINT,1=',T1(3*(JDO-1)+2),COST,SINT,COST**2+SINT**2
      IF (ABS(COST**2+SINT**2-1.0D0).GT.2.0D-3) THEN
         PRINT '(A,G20.10)','ERROR - in ROTXZ cos**2+sin**2=',COST**2+SINT**2
         STOP
      ENDIF
      RVEC(1)=0.0D0
      RVEC(2)=0.0D0
      RVEC(3)=1.0D0
      DO J1=1,NATOMS
         IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=T1(3*(J1-1)+1)*RVEC(1)+T1(3*(J1-1)+2)*RVEC(2)+T1(3*(J1-1)+3)*RVEC(3)
            TX=T1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)+(T1(3*(J1-1)+2)*RVEC(3)-T1(3*(J1-1)+3)*RVEC(2))*SINT
            TY=T1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)+(T1(3*(J1-1)+3)*RVEC(1)-T1(3*(J1-1)+1)*RVEC(3))*SINT
            TZ=T1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)+(T1(3*(J1-1)+1)*RVEC(2)-T1(3*(J1-1)+2)*RVEC(1))*SINT
            T1(3*(J1-1)+1)=TX
            T1(3*(J1-1)+2)=TY
            T1(3*(J1-1)+3)=TZ
         ENDIF
      ENDDO

20    CONTINUE

      PROJ=0.0D0
      DO J1=1,NATOMS
         IF (T1(3*(J1-1)+3).GT.1.0D-2) PROJ=PROJ+T1(3*(J1-1)+1)
      ENDDO
C     PRINT '(I6)',NATOMS
C     PRINT '(A,I6,G20.10)','after rotation atom,PROJ=',JDO,PROJ
C     WRITE(*,'(A2,2X,3G20.10)') ('LA',T1(3*(J3-1)+1),T1(3*(J3-1)+2),T1(3*(J3-1)+3),J3=1,NATOMS)

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Put permutational isomers into a standard orientation.
C
      SUBROUTINE ORIENT2D(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION Q1(3*NATOMS), DIST(NATOMS), DMAX, RVEC(3), T1(3*NATOMS), CMX, CMY, CMZ,
     1                 COST, SINT, RDOTN, TX, TY, TZ
      DOUBLE PRECISION, PARAMETER :: ORBDISTTOL=1.0D-2
      INTEGER J1, I, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
C
C  Move centre of mass to the origin.
C
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+Q1(3*(I-1)+1)
         CMY=CMY+Q1(3*(I-1)+2)
         CMZ=CMZ+Q1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
C     PRINT*,'CMX,CMY,CMZ=',CMX,CMY,CMZ
      DO I=1,NATOMS
         Q1(3*(I-1)+1)=Q1(3*(I-1)+1)-CMX
         Q1(3*(I-1)+2)=Q1(3*(I-1)+2)-CMY
         Q1(3*(I-1)+3)=Q1(3*(I-1)+3)-CMZ
      ENDDO
      
      NORBIT1=1

      DO J1=1,3*NATOMS
         T1(J1)=Q1(J1)
      ENDDO

C     IF (DEBUG) THEN
C        PRINT*,'before'
C        WRITE(*,'(3F20.10)') (Q1(J1),J1=1,3*NATOMS)
C        PRINT*,'after'
C        WRITE(*,'(3F20.10)') (T1(J1),J1=1,3*NATOMS)
C     ENDIF
C
C  Find the atom with the largest distance from the z axis.
C
      DMAX=-1.0D0
      NORBIT2=1
      DO J1=1,NATOMS
         DIST(J1)=SQRT(T1(3*(J1-1)+1)**2+T1(3*(J1-1)+2)**2)
C        WRITE(*,'(A,I5,3F20.10)') 'J1,DIST,DMAX,ABS(diff)=',J1,DIST(J1),DMAX,ABS(DIST(J1)-DMAX)
         IF (ABS(DIST(J1)-DMAX).LT.ORBDISTTOL) THEN
            NORBIT2=NORBIT2+1
            IF (NORBIT2.EQ.NCHOOSE2) THEN
               JMAX2=J1
            ENDIF
         ELSE IF (DIST(J1).GT.DMAX) THEN
            DMAX=DIST(J1)
            NORBIT2=1
            JMAX2=J1
         ENDIF
      ENDDO
C     IF (DEBUG) PRINT*,'NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX2=',NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX2
C
C  and rotate it into the xz plane.
C
C     IF (DEBUG) PRINT*,'atom ',JMAX2,' is now the furthest from the  z axis, distance=',DMAX,' rotate into xz plane'
      IF (T1(3*(JMAX2-1)+2).EQ.0.0D0) GOTO 20
      COST=T1(3*(JMAX2-1)+1)/DMAX
      SINT=T1(3*(JMAX2-1)+2)/DMAX
C     PRINT*,'COST,SINT=',COST,SINT
      RVEC(1)=0.0D0
      RVEC(2)=0.0D0
      RVEC(3)=1.0D0
      DO J1=1,NATOMS
         IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=T1(3*(J1-1)+1)*RVEC(1)+T1(3*(J1-1)+2)*RVEC(2)+T1(3*(J1-1)+3)*RVEC(3)
            TX=T1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)+(T1(3*(J1-1)+2)*RVEC(3)-T1(3*(J1-1)+3)*RVEC(2))*SINT
            TY=T1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)+(T1(3*(J1-1)+3)*RVEC(1)-T1(3*(J1-1)+1)*RVEC(3))*SINT
            TZ=T1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)+(T1(3*(J1-1)+1)*RVEC(2)-T1(3*(J1-1)+2)*RVEC(1))*SINT
            T1(3*(J1-1)+1)=TX
            T1(3*(J1-1)+2)=TY
            T1(3*(J1-1)+3)=TZ
         ENDIF
      ENDDO

20    CONTINUE

C     IF (DEBUG) PRINT*,'after'
C     IF (DEBUG) WRITE(*,'(3F20.10)') (T1(J1),J1=1,3*NATOMS)

      RETURN
      END
