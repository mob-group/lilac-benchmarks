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
C jmc hacked to work in internal coordinates (with unres)
C NOPT will be set to the value for Cartesians, NINTS is for internals
C DUMMY3, COORDS, VEC will be in internals
C GRAD1,GRAD2 will be returned from potential in internals
C Therefore only thing that should really be affected is the call to potential 

      SUBROUTINE INTSECDIAG(VEC,COORDS,ENERGY,GL,DIAG,GTEST,XRMS)
      USE COMMONS
      USE KEY
      IMPLICIT NONE

      INTEGER J1
      LOGICAL GTEST, FPLUS, FMINUS
      DOUBLE PRECISION DUMMY3(NINTS), DIFF, VEC(3*NATOMS), COORDS(3*NATOMS), DIAG2, DIAG3,
     1                 EPLUS,EMINUS,ENERGY,GRAD1(3*NATOMS),LOCALV(NINTS),PROJ,
     2                 RMS,DIAG,XRMS,GL(3*NATOMS),GRAD2(3*NATOMS), VECL, ZETA, DUMCART(3*NATOMS)

      DIFF=1.0D-3
      IF (UNRST) DIFF=5.0D-4

      LOCALV(1:NINTS)=VEC(1:NINTS)
      CALL VECNORM(LOCALV(1:NINTS),NINTS)

C jmc      IF (NFREEZE.LT.3) CALL ORTHOGOPT(VEC,COORDS,.TRUE.)
C to orthogonalise the displacement to translations and rotations

      DUMCART=1.0D0

C     VECL=0.0D0
C     DO J1=1,NINTS
C        VECL=VECL+VEC(J1)**2
C     ENDDO
C     ZETA=DIFF/SQRT(VECL)

      VECL=1.0D0
      ZETA=DIFF

      DUMMY3(1:NINTS)=COORDS(1:NINTS)+ZETA*LOCALV(1:NINTS)

C     WRITE(*,'(6F15.5)') (DUMMY3(J1),J1=1,NINTS)

      IF (UNRST) THEN
         CALL var_to_geom(NINTS,DUMMY3)
         CALL chainbuild
      END IF
      CALL POTENTIAL(DUMCART,EPLUS,GRAD1,GTEST,.FALSE.,RMS,.FALSE.,.FALSE.)

      DUMMY3(1:NINTS)=COORDS(1:NINTS)-ZETA*LOCALV(1:NINTS)

C     WRITE(*,'(6F15.5)') (DUMMY3(J1),J1=1,NINTS)

      IF (UNRST) THEN
         CALL var_to_geom(NINTS,DUMMY3)
         CALL chainbuild
      END IF
C jmc remember for unres, passing of coords array (first arg) is irrelevant...
      CALL POTENTIAL(DUMCART,EMINUS,GRAD2,GTEST,.FALSE.,RMS,.FALSE.,.FALSE.)

      DIAG=(EPLUS+EMINUS-2.0D0*ENERGY)/((ZETA**2)*VECL) ! jmc DIAG is lambda*vecl (but vecl=1.0d0 anyway...)

      DIAG2=0.0D0
      DO J1=1,NINTS
         DIAG2=DIAG2+(GRAD1(J1)-GRAD2(J1))*LOCALV(J1)
C        WRITE(*,'(A,I4,4F20.10)') 'J1,GRAD1,GRAD2,LOCALV,DIAG2=',J1,GRAD1(J1),GRAD2(J1),LOCALV(J1),DIAG2
      ENDDO
      DIAG2=DIAG2/(2.0D0*ZETA)
      DIAG3=2*(DIAG-DIAG2/2)
C     WRITE(*,'(A,6F20.10)') 'D,D2,D3,E+,E-,E=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY
C     IF (.NOT.GTEST) WRITE(*,'(A,6F20.10)') 'D,D2,D3,E+,E-,E=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY
C
C  Although DIAG3 is a more accurate estimate of the diagonal second derivative, it
C  cannot be differentiated analytically.
C
      IF (GTEST) THEN
C        DO J1=1,NINTS
C this is from eqn 6.20 in Chapter 6 of 'The Book'...
C           GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG*LOCALV(J1)/VECL**2
C           WRITE(*,'(A,I4,4G16.7)') 'J1,GRAD1,GRAD2,VEC,GL=',J1,GRAD1(J1),GRAD2(J1),LOCALV(J1),GL(J1)
C        ENDDO
         GL(1:NINTS)=(GRAD1(1:NINTS)-GRAD2(1:NINTS))/(ZETA*VECL**2)-2.0D0*DIAG*LOCALV(1:NINTS)/VECL**2
C jmc so GL is vecl*dlambda/dx...
C to orthogonalise the **gradient** to translations and rotations (obviously not necessary...)
C        CALL ORTHOGOPT(GL,COORDS,.FALSE.)

C
C  Project out any component of the gradient along LOCALV (which is a unit vector).
C
         PROJ=0.0D0
         DO J1=1,NINTS
            PROJ=PROJ+GL(J1)*LOCALV(J1)
         ENDDO
         DO J1=1,NINTS 
            GL(J1)=GL(J1)-PROJ*LOCALV(J1)
         ENDDO

         XRMS=0.0D0
         DO J1=1,NINTS
            XRMS=XRMS+GL(J1)**2
         ENDDO
         XRMS=DSQRT(XRMS/NINTS)
         IF (DEBUG) WRITE(*,'(A,3G15.5,3G20.12,G10.3)') 'D,D2,D3,E+,E-,E,RMS=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY,XRMS
         IF (DEBUG) WRITE(*,'(A,G20.10)') 'predicted gradient component=',(EPLUS-EMINUS)/(2*ZETA)
C        WRITE(*,'(A,5G20.10)') 'D,E+,E-,E,RMS=',DIAG,EPLUS,EMINUS,ENERGY,XRMS
      ENDIF
CC
CC  Project out any component of the gradient along LOCALV (which is a unit vector).
CC
C      PROJ=0.0D0
C      DO J1=1,NINTS
C         PROJ=PROJ+GL(J1)*LOCALV(J1)
C      ENDDO
C      DO J1=1,NINTS
C         GL(J1)=GL(J1)-PROJ*LOCALV(J1)
C      ENDDO

      RETURN
      END
