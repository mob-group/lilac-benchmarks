C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales GNU General Public License
C   This file is part of OPTIM. not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C     SUBROUTINE THOMSON(X,V,ETHOMSON,GTEST)
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C     DOUBLE PRECISION X(*), DIST, V(*), ETHOMSON, DUMMY, CT1, ST1, CT2, ST2, CPDIFF, SPDIFF
C   You should have received a copy of the GNU General Public License(NATOMS)
C   along with this program; if not, write to the Free Software0
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C  Energy, gradient and second derivatives for the Thomson problem in theta, phi coordinates.
C
      SUBROUTINE THOMSON(NATOMS,X,V,ETHOMSON,GTEST,STEST)
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST, STEST
      INTEGER J1, J2, J3, J4, NATOMS, LNATOMS
      DOUBLE PRECISION X(*), DIST, V(*), ETHOMSON, DUMMY, CT1, ST1, CT2, ST2, CPDIFF, SPDIFF, CP1, SP1, CP2, SP2
      DOUBLE PRECISION COST((NATOMS/2)*3), SINT((NATOMS/2)*3), COSP((NATOMS/2)*3), SINP((NATOMS/2)*3)
      DOUBLE PRECISION, PARAMETER :: SR2=1.4142135623730950488D0

      LNATOMS=(NATOMS/2)*3
      DO J1=1,LNATOMS
         J3=2*J1
         COST(J1)=COS(X(J3-1))
         SINT(J1)=SIN(X(J3-1))
         COSP(J1)=COS(X(J3))
         SINP(J1)=SIN(X(J3))
      ENDDO

      ETHOMSON=0.0D0
      V(1:2*LNATOMS)=0.0D0
      DO J1=1,LNATOMS
         J3=2*J1
         CT1=COST(J1)
         ST1=SINT(J1)
         DO J2=J1+1,LNATOMS
            J4=2*J2
            CT2=COST(J2)
            ST2=SINT(J2)
            CPDIFF=COSP(J1)*COSP(J2)+SINP(J1)*SINP(J2)
            SPDIFF=COSP(J2)*SINP(J1)-SINP(J2)*COSP(J1)
            DIST=1.0D0/SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
            ETHOMSON=ETHOMSON+DIST
            DIST=(DIST/SR2)**3
            DUMMY=SPDIFF*ST1*ST2*DIST
            V(J3-1)=V(J3-1)+(CPDIFF*CT1*ST2-CT2*ST1)*DIST
            V(J3)=V(J3)    -DUMMY
            V(J4-1)=V(J4-1)+(CPDIFF*CT2*ST1-CT1*ST2)*DIST
            V(J4)=V(J4)    +DUMMY
         ENDDO
      ENDDO
      ETHOMSON=ETHOMSON/SR2

      IF (.NOT.STEST) RETURN

      HESS(1:2*LNATOMS,1:2*LNATOMS)=0.0D0
      DO J1=1,LNATOMS
         J3=2*J1
         CT1=COST(J1)
         ST1=SINT(J1)
         CP1=COSP(J1)
         SP1=SINP(J1)
         DO J2=1,J1-1
            J4=2*J2
            CT2=COST(J2)
            ST2=SINT(J2)
            CP2=COSP(J2)
            SP2=SINP(J2)
            CPDIFF=COSP(J1)*COSP(J2)+SINP(J1)*SINP(J2)
            SPDIFF=COSP(J2)*SINP(J1)-SINP(J2)*COSP(J1)
            DIST=SR2*SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
            DIST=DIST**5
            HESS(J3-1,J3-1)=HESS(J3-1,J3-1)+
     &           (-(ct2*(4*ct1 + ct1**2*ct2 - ct2*(5 + st1**2))) - 4*cpdiff*(1 + ct1*ct2)*st1*st2 + 
     &            cpdiff**2*(5 + ct1**2 - st1**2)*st2**2)/(2.*dist)
            HESS(J3,J3)=HESS(J3,J3)+
     &           (st1*st2*(4*cpdiff*(-1 + ct1*ct2) + (5 - (sp1*(-cp2 + sp2) + cp1*(cp2 + sp2))*(cp1*(cp2 - sp2) + 
     &           sp1*(cp2 + sp2)))*st1*st2))/(2.*dist)
            HESS(J3-1,J3)=HESS(J3-1,J3)+
     &             (spdiff*st2*(2*ct1**2*ct2 + 3*ct2*st1**2 - ct1*(2 + cpdiff*st1*st2)))/dist
         ENDDO
         DO J2=J1+1,LNATOMS
            J4=2*J2
            CT2=COST(J2)
            ST2=SINT(J2)
            CP2=COSP(J2)
            SP2=SINP(J2)
            CPDIFF=COSP(J1)*COSP(J2)+SINP(J1)*SINP(J2)
            SPDIFF=COSP(J2)*SINP(J1)-SINP(J2)*COSP(J1)
            DIST=SR2*SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
            DIST=DIST**2
            HESS(J3-1,J3-1)=HESS(J3-1,J3-1)+
     &           (-(ct2*(4*ct1 + ct1**2*ct2 - ct2*(5 + st1**2))) - 4*cpdiff*(1 + ct1*ct2)*st1*st2 + 
     &            cpdiff**2*(5 + ct1**2 - st1**2)*st2**2)/(2.*dist)
            HESS(J3,J3)=HESS(J3,J3)+
     &           (st1*st2*(4*cpdiff*(-1 + ct1*ct2) + (5 - (sp1*(-cp2 + sp2) + cp1*(cp2 + sp2))*(cp1*(cp2 - sp2) + 
     &           sp1*(cp2 + sp2)))*st1*st2))/(2.*dist)
            HESS(J3-1,J3)=HESS(J3-1,J3)+
     &             (spdiff*st2*(2*ct1**2*ct2 + 3*ct2*st1**2 - ct1*(2 + cpdiff*st1*st2)))/dist
            HESS(J3-1,J4-1)=(cpdiff*(-5 + ct2*(ct1*(4 + ct1*ct2) - ct2*st1**2)) + 2*(2 + (1 + cpdiff**2)*ct1*ct2)*st1*st2 - 
     &              cpdiff*(ct1 - st1)*(ct1 + st1)*st2**2)/(2.*dist)
            HESS(J3-1,J4)=(spdiff*st2*(-2*ct1**2*ct2 - 3*ct2*st1**2 + ct1*(2 + cpdiff*st1*st2)))/dist
            HESS(J3,J4-1)=(spdiff*st1*(2*ct1*ct2**2 + 3*ct1*st2**2 - ct2*(2 + cpdiff*st1*st2)))/dist
            HESS(J3,J4)=(st1*st2*(cpdiff*(4 - 4*ct1*ct2) + (-5 + (sp1*(-cp2 + sp2) + cp1*(cp2 + sp2))*(cp1*(cp2 - sp2) + 
     &                  sp1*(cp2 + sp2)))*st1*st2))/(2.*dist)
         ENDDO
      ENDDO
      DO J1=1,2*LNATOMS
         DO J2=J1+1,2*LNATOMS
            HESS(J2,J1)=HESS(J1,J2)
         ENDDO
      ENDDO

      RETURN
      END

C
C  Orthogonalise VEC1 to overall rotations about the x, y, and z axes.
C
      SUBROUTINE ORTHOGTH(VEC1,Q,OTEST)
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2, J3
      DOUBLE PRECISION Q(*), VEC1(*), DUMMY1, VECX(3*NATOMS), VECY(3*NATOMS), VECZ(3*NATOMS)
      LOGICAL OTEST

      VECX(1:3*NATOMS)=0.0D0; VECY(1:3*NATOMS)=0.0D0; VECZ(1:3*NATOMS)=0.0D0
      DO J1=1,3*NATOMS,2
         VECX(J1)=SIN(Q(J1+1))
         VECY(J1)=COS(Q(J1+1))
      ENDDO
      DO J1=2,3*NATOMS,2
         IF (SIN(Q(J1)).NE.0.0D0) THEN
            VECX(J1)= COS(Q(J1-1))*COS(Q(J1))/SIN(Q(J1-1))
            VECY(J1)=-COS(Q(J1-1))*SIN(Q(J1))/SIN(Q(J1-1))
         ENDIF
         VECZ(J1)=1.0D0
      ENDDO
      CALL VECNORM(VECX,3*NATOMS)
      CALL VECNORM(VECY,3*NATOMS)
      CALL VECNORM(VECZ,3*NATOMS)

      DUMMY1=0.0D0
      DO J2=1,3*NATOMS
         DUMMY1=DUMMY1+VEC1(J2)*VECX(J2)
      ENDDO
      DO J2=1,3*NATOMS
         VEC1(J2)=VEC1(J2)-DUMMY1*VECX(J2)
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,3*NATOMS)

      DUMMY1=0.0D0
      DO J2=1,3*NATOMS
         DUMMY1=DUMMY1+VEC1(J2)*VECY(J2)
      ENDDO
      DO J2=1,3*NATOMS
         VEC1(J2)=VEC1(J2)-DUMMY1*VECY(J2)
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,3*NATOMS)

      DUMMY1=0.0D0
      DO J2=1,3*NATOMS
         DUMMY1=DUMMY1+VEC1(J2)*VECZ(J2)
      ENDDO
      DO J2=1,3*NATOMS
         VEC1(J2)=VEC1(J2)-DUMMY1*VECZ(J2)
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,3*NATOMS)

      RETURN
      END
C
C  Eigenvalue shifting for the Thomson problem
C
      SUBROUTINE SHIFTHTH(Q,NATOMS)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, NATOMS
      DOUBLE PRECISION DUMMY, Q(3*NATOMS), VECX(3*NATOMS), VECY(3*NATOMS), VECZ(3*NATOMS)

      SHIFTED=.TRUE.

      VECX(1:3*NATOMS)=0.0D0; VECY(1:3*NATOMS)=0.0D0; VECZ(1:3*NATOMS)=0.0D0
      DO J1=1,3*NATOMS,2
         VECX(J1)=SIN(Q(J1+1))
         VECY(J1)=COS(Q(J1+1))
      ENDDO
      DO J1=2,3*NATOMS,2
         IF (SIN(Q(J1)).NE.0.0D0) THEN
            VECX(J1)= COS(Q(J1-1))*COS(Q(J1))/SIN(Q(J1-1))
            VECY(J1)=-COS(Q(J1-1))*SIN(Q(J1))/SIN(Q(J1-1))
         ENDIF
         VECZ(J1)=1.0D0
      ENDDO
      CALL VECNORM(VECX,3*NATOMS)
      CALL VECNORM(VECY,3*NATOMS)
      CALL VECNORM(VECZ,3*NATOMS)

      DO J1=1,3*NATOMS
         DO J2=1,3*NATOMS
            HESS(J2,J1)=HESS(J2,J1)+SHIFTL(1)*VECX(J2)*VECX(J1)+SHIFTL(2)*VECY(J2)*VECY(J1)+SHIFTL(3)*VECZ(J2)*VECZ(J1)
         ENDDO
      ENDDO

      RETURN
      END
