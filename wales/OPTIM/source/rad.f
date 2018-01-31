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
C  Add energy and gradient correction terms for the container separately.
C
C     SUBROUTINE RAD(X,ENERGY,GRAD,GTEST)
      SUBROUTINE RAD(X)
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      INTEGER J1, J3, NEIGH, J2
      DOUBLE PRECISION X(3*NATOMS), DIST, CMX, CMY, CMZ, NEARD, DUMMY, ALPHA
      LOGICAL RADMOVED, OVERLP
      COMMON /DISCON/ RADMOVED, OVERLP

      RADMOVED=.FALSE.
      IF (BULKT) RETURN

      ALPHA=1.0D0

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO J1=1,NATOMS
         CMX=CMX+X(3*(J1-1)+1)
         CMY=CMY+X(3*(J1-1)+2)
         CMZ=CMZ+X(3*(J1-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
C     CMX=0.0D0
C     CMY=0.0D0
C     CMZ=0.0D0
      DO J1=1,NATOMS
         J3=3*J1
         DIST=(X(J3-2)-CMX)**2+(X(J3-1)-CMY)**2+(X(J3)-CMZ)**2
         IF (DIST.GT.RADIUS) THEN
            PRINT*,'In RAD DIST,RADIUS=',SQRT(DIST),SQRT(RADIUS)

C            ENERGY=ENERGY+ALPHA*(SQRT(DIST)-SQRT(RADIUS))**2/2.0D0

C            IF (GTEST) THEN
C               GRAD(J3-2)=GRAD(J3-2)-ALPHA*(SQRT(DIST)-SQRT(RADIUS))*(X(J3-2)-CMX)/SQRT(DIST)
C               GRAD(J3-1)=GRAD(J3-1)-ALPHA*(SQRT(DIST)-SQRT(RADIUS))*(X(J3-1)-CMY)/SQRT(DIST)
C               GRAD(J3)=  GRAD(J3)  -ALPHA*(SQRT(DIST)-SQRT(RADIUS))*(X(J3)-CMZ)/SQRT(DIST)
C            ENDIF

            NEARD=1.0D90
            DO J2=1,NATOMS
               IF (J2.NE.J1) THEN
                  DUMMY=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2+(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2+(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
                  IF (DUMMY.LT.NEARD) THEN
                     NEARD=DUMMY
                     NEIGH=J2
                  ENDIF
               ENDIF
            ENDDO
            NEARD=SQRT(NEARD)
C           IF (NEARD.GT.1.2D0) THEN
C              PRINT*,'coords before:'
C              WRITE(*,'(3F20.10)') (X(J2),J2=1,3*NATOMS)
C              PRINT*,'moving atom ',J1,' nearer to atom ',NEIGH,' NEARD=',NEARD
C              PRINT*,'atom ',J1,' initially at ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3)
C              PRINT*,'atom ',NEIGH,'           at ',X(3*(NEIGH-1)+1),X(3*(NEIGH-1)+2),X(3*(NEIGH-1)+3)
               X(3*(J1-1)+1)=X(3*(NEIGH-1)+1)+(X(3*(J1-1)+1)-X(3*(NEIGH-1)+1))*1.2D0/NEARD
               X(3*(J1-1)+2)=X(3*(NEIGH-1)+2)+(X(3*(J1-1)+2)-X(3*(NEIGH-1)+2))*1.2D0/NEARD
               X(3*(J1-1)+3)=X(3*(NEIGH-1)+3)+(X(3*(J1-1)+3)-X(3*(NEIGH-1)+3))*1.2D0/NEARD
               X(3*(J1-1)+1)=X(3*(NEIGH-1)+1)+(X(3*(J1-1)+1)-X(3*(NEIGH-1)+1))*0.95D0
               X(3*(J1-1)+2)=X(3*(NEIGH-1)+2)+(X(3*(J1-1)+2)-X(3*(NEIGH-1)+2))*0.95D0
               X(3*(J1-1)+3)=X(3*(NEIGH-1)+3)+(X(3*(J1-1)+3)-X(3*(NEIGH-1)+3))*0.95D0
               RADMOVED=.TRUE.
C              PRINT*,'atom ',J1,'       now at ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3)
C              PRINT*,'DIST is now=',SQRT((X(J3-2)-CMX)**2+(X(J3-1)-CMY)**2+(X(J3)-CMZ)**2)
C              PRINT*,'coords after:'
C              WRITE(*,'(3F20.10)') (X(J2),J2=1,3*NATOMS)
C           ENDIF
C           DIST=SQRT(RADIUS*0.7D0/DIST)
C           X(J3-2)=CMX+(X(J3-2)-CMX)*DIST
C           X(J3-1)=CMY+(X(J3-1)-CMY)*DIST
C           X(J3)=CMZ+(X(J3)-CMZ)*DIST
         ENDIF
      ENDDO

      RETURN
      END
