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
C*************************************************************************
C
C  Subroutine LJADD calculates the cartesian gradient and second
C  derivative matrix analytically for LJ with addressable epsilon values. Reduced units.
C
C*************************************************************************
C
      SUBROUTINE LJADD(N, X, V, ENERGY, GTEST, STEST)
      USE KEY, ONLY : LJADDEPS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, R6,
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N), XG(N,N),
     3                 R14(N,N), F(N,N), DUMMY, DUMMYX, DUMMYY, DUMMYZ, DIST, XMUL2
C 
C  Store distance matrices.
C
      ENERGY=0.0D0
      IF (GTEST.AND.(.NOT.STEST)) THEN
         DO J1=1,N
            J3=3*J1
            XG(J1,J1)=0.0D0
            DO J2=J1+1,N
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=LJADDEPS(J2,J1)*R6*(R6-1.0D0)
               ENERGY=ENERGY+DUMMY
               DIST=DIST*R6
               XG(J2,J1)=-LJADDEPS(J2,J1)*24.0D0*(2.0D0*R6-1.0D0)*DIST
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO
      ELSEIF (GTEST) THEN
         DO J1=1,N
            XG(J1,J1)=0.0D0
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               ENERGY=ENERGY+LJADDEPS(J2,J1)*R6*(R6-1.0D0)
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               XG(J2,J1)=-LJADDEPS(J2,J1)*24.0D0*(2.0D0*R6-1.0D0)*R2(J1,J2)*R6
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO 
      ELSE
         DO J1=1,N
            J3=3*(J1-1)
            DO J2=J1+1,N
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6=R2T**3
               ENERGY=ENERGY+LJADDEPS(J2,J1)*R6*(R6-1.0D0)
            ENDDO
         ENDDO

      ENDIF
      ENERGY=4.0D0*ENERGY

      IF (.NOT.GTEST) RETURN
      DO J1=1,N
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,N
            J2=3*J4
            XMUL2=XG(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO

      IF (.NOT.STEST) RETURN
      CALL LJADDS(G,F,R2,R14,R8,X,N)

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJADDS(G,F,R2,R14,R8,X,N)
      USE MODHESS
      USE KEY, ONLY : LJADDEPS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), 
     2                 X(3*N),DUMMY

C
C  Calculate the g tensor.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            G(J2,J1)=-LJADDEPS(J2,J1)*24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=LJADDEPS(J2,J1)*672.0D0*R14(J2,J1)-LJADDEPS(J2,J1)*192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END
C
C*************************************************************************
C
C  Subroutine LJADD2 calculates the cartesian gradient and second
C  derivative matrix analytically for LJ with addressable epsilon values. Reduced units.
C  This routine treats multiple copies of a target of cluster size NADDTARGET.
C  The epsilon values are replicated via the MOD function.
C
C*************************************************************************
C
      SUBROUTINE LJADD2(N, X, V, ENERGY, GTEST, STEST)
      USE KEY, ONLY : LJADDEPS, NADDTARGET
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, MJ1, MJ2
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, R6,
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N), XG(N,N),
     3                 R14(N,N), F(N,N), DUMMY, DUMMYX, DUMMYY, DUMMYZ, DIST, XMUL2
C
C  Store distance matrices.
C
!     WRITE(MYUNIT,'(A)') 'coords in LJADD2:'
!     WRITE(MYUNIT,'(3G20.10)') X(1:3*N)
      ENERGY=0.0D0
      IF (GTEST.AND.(.NOT.STEST)) THEN
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            J3=3*J1
            XG(J1,J1)=0.0D0
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               IF (MJ1.EQ.MJ2) THEN ! repulsive part only
                  DUMMY=LJADDEPS(MJ2,MJ1)*R6*R6
               ELSE
                  DUMMY=LJADDEPS(MJ2,MJ1)*R6*(R6-1.0D0)
               ENDIF
               ENERGY=ENERGY+DUMMY
!              WRITE(MYUNIT,'(A,2I8,2G20.10)') 'J1,J2,DUMMY,ENERGY=',J1,J2,DUMMY,ENERGY
               DIST=DIST*R6
               IF (MJ1.EQ.MJ2) THEN ! repulsive part only
                  XG(J2,J1)=-LJADDEPS(MJ2,MJ1)*24.0D0*2.0D0*R6*DIST
               ELSE
                  XG(J2,J1)=-LJADDEPS(MJ2,MJ1)*24.0D0*(2.0D0*R6-1.0D0)*DIST
               ENDIF
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO
      ELSEIF (GTEST) THEN
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            XG(J1,J1)=0.0D0
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               IF (MJ1.EQ.MJ2) THEN ! repulsive part only
                  DUMMY=LJADDEPS(MJ2,MJ1)*R6*R6
               ELSE
                  DUMMY=LJADDEPS(MJ2,MJ1)*R6*(R6-1.0D0)
               ENDIF
               ENERGY=ENERGY+DUMMY
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               IF (MJ1.EQ.MJ2) THEN ! repulsive part only
                  XG(J2,J1)=-LJADDEPS(MJ2,MJ1)*24.0D0*2.0D0*R6*R2(J1,J2)*R6
               ELSE
                  XG(J2,J1)=-LJADDEPS(MJ2,MJ1)*24.0D0*(2.0D0*R6-1.0D0)*R2(J1,J2)*R6
               ENDIF
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            J3=3*(J1-1)
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6=R2T**3
               IF (MJ1.EQ.MJ2) THEN ! repulsive part only
                  ENERGY=ENERGY+LJADDEPS(MJ2,MJ1)*R6*R6
               ELSE
                  ENERGY=ENERGY+LJADDEPS(MJ2,MJ1)*R6*(R6-1.0D0)
               ENDIF
            ENDDO
         ENDDO

      ENDIF
      ENERGY=4.0D0*ENERGY

      IF (.NOT.GTEST) RETURN
      DO J1=1,N
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,N
            J2=3*J4
            XMUL2=XG(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO

      IF (.NOT.STEST) RETURN
      CALL LJADDS2(G,F,R2,R14,R8,X,N)

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJADDS2(G,F,R2,R14,R8,X,N)
      USE MODHESS
      USE KEY, ONLY : LJADDEPS, NADDTARGET
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, MJ1, MJ2
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N),
     2                 X(3*N),DUMMY

C
C  Calculate the g tensor.
C
      DO J1=1,N
         MJ1=MOD(J1-1,NADDTARGET)+1
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            MJ2=MOD(J2-1,NADDTARGET)+1
            IF (MJ1.EQ.MJ2) THEN ! repulsive part only
               G(J2,J1)=-LJADDEPS(MJ2,MJ1)*24.0D0*2.0D0*R14(J2,J1)
            ELSE
               G(J2,J1)=-LJADDEPS(MJ2,MJ1)*24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
            ENDIF
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,N
         MJ1=MOD(J1-1,NADDTARGET)+1
         F(J1,J1)=0.0D0
         DO J2=J1+1,N
            MJ2=MOD(J2-1,NADDTARGET)+1
            IF (MJ1.EQ.MJ2) THEN ! repulsive part only
               F(J2,J1)=LJADDEPS(MJ2,MJ1)*672.0D0*R14(J2,J1)
            ELSE
               F(J2,J1)=LJADDEPS(MJ2,MJ1)*672.0D0*R14(J2,J1)-LJADDEPS(MJ2,MJ1)*192.0D0*R8(J2,J1)
            ENDIF
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)*
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1)
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END

C
C*************************************************************************
C
C  Subroutine LJADD3 calculates the cartesian gradient and second
C  derivative matrix analytically for LJ with addressable epsilon values. Reduced units.
C  This routine treats multiple copies of a target of cluster size NADDTARGET.
C  The epsilon values are replicated via the MOD function.
C  In LJADD3 we have separate scaling for repulsion and attraction.
C
C*************************************************************************
C
      SUBROUTINE LJADD3(N, X, V, ENERGY, GTEST, STEST)
      USE KEY, ONLY : LJADDREP, LJADDATT, NADDTARGET, MYUNIT
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, MJ1, MJ2
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, R6,
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N), XG(N,N),
     3                 R14(N,N), F(N,N), DUMMY, DUMMYX, DUMMYY, DUMMYZ, DIST, XMUL2
C 
C  Store distance matrices.
C
!     WRITE(MYUNIT,'(A)') 'coords in LJADD3:'
!     WRITE(MYUNIT,'(3G20.10)') X(1:3*N)
      ENERGY=0.0D0
      IF (GTEST.AND.(.NOT.STEST)) THEN
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            J3=3*J1
            XG(J1,J1)=0.0D0
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATT(MJ2,MJ1))
               ENERGY=ENERGY+DUMMY
               DIST=DIST*R6
               XG(J2,J1)=-24.0D0*(2.0D0*R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATT(MJ2,MJ1))*DIST
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO
      ELSEIF (GTEST) THEN
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            XG(J1,J1)=0.0D0
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               DUMMY=R6*(R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATT(MJ2,MJ1))
               ENERGY=ENERGY+DUMMY
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               XG(J2,J1)=-24.0D0*(2.0D0*R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATT(MJ2,MJ1))*R2(J1,J2)*R6
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO 
      ELSE
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            J3=3*(J1-1)
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6=R2T**3
               ENERGY=ENERGY+R6*(R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATT(MJ2,MJ1))
            ENDDO
         ENDDO

      ENDIF
      ENERGY=4.0D0*ENERGY

      IF (.NOT.GTEST) RETURN
      DO J1=1,N
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,N
            J2=3*J4
            XMUL2=XG(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO

      IF (.NOT.STEST) RETURN
      CALL LJADDS3(G,F,R2,R14,R8,X,N)

      RETURN
      END

C
C*************************************************************************
C
C  Subroutine LJADD4 calculates the cartesian gradient and second
C  derivative matrix analytically for LJ with addressable epsilon values. Reduced units.
C  This routine treats multiple copies of a target of cluster size NADDTARGET.
C  The epsilon values are replicated via the MOD function.
C  In LJADD4 we have separate scaling for repulsion and attraction.
C  The nearest-neighbour attractive term is scaled in LJADD4 by a cooperative term.
C
C*************************************************************************
C
      SUBROUTINE LJADD4(N, X, V, ENERGY, GTEST, STEST)
      USE KEY, ONLY : LJADDREP, LJADDATT, NADDTARGET, MYUNIT, NUMNN, TANHFAC, LJADDCUTOFF, LJADDNN, LJADDREFNORM
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, MJ1, MJ2, MJ3, NCOPIES, K1, K2
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, R6,
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N), XG(N,N),
     3                 R14(N,N), F(N,N), DUMMY, DUMMYX, DUMMYY, DUMMYZ, DIST, XMUL2,
     &                 LJADDATTLOCAL(NADDTARGET,NADDTARGET), DUMMY2, DUMMYEPS
      DOUBLE PRECISION GEXTRA(3*N)

      NCOPIES=N/NADDTARGET
      DUMMYEPS=0.0D0
      GEXTRA(1:3*N)=0.0D0
      DO J1=1,NUMNN
         DO J2=1,NCOPIES
            MJ2=LJADDNN(J1,1)+(J2-1)*NADDTARGET
            DO J3=1,NCOPIES
               MJ3=LJADDNN(J1,2)+(J3-1)*NADDTARGET
               DIST=SQRT((X(3*(MJ2-1)+1)-X(3*(MJ3-1)+1))**2 + 
     &                   (X(3*(MJ2-1)+2)-X(3*(MJ3-1)+2))**2 +
     &                   (X(3*(MJ2-1)+3)-X(3*(MJ3-1)+3))**2)
               DUMMYEPS=DUMMYEPS+TANH(TANHFAC*(LJADDCUTOFF-DIST))+1.0D0
               DUMMY2=(1.0D0/COSH(TANHFAC*(LJADDCUTOFF-DIST))**2)/DIST
               GEXTRA(3*(MJ2-1)+1)=GEXTRA(3*(MJ2-1)+1)-DUMMY2*(X(3*(MJ2-1)+1)-X(3*(MJ3-1)+1))
               GEXTRA(3*(MJ2-1)+2)=GEXTRA(3*(MJ2-1)+2)-DUMMY2*(X(3*(MJ2-1)+2)-X(3*(MJ3-1)+2))
               GEXTRA(3*(MJ2-1)+3)=GEXTRA(3*(MJ2-1)+3)-DUMMY2*(X(3*(MJ2-1)+3)-X(3*(MJ3-1)+3))
               GEXTRA(3*(MJ3-1)+1)=GEXTRA(3*(MJ3-1)+1)-DUMMY2*(X(3*(MJ3-1)+1)-X(3*(MJ2-1)+1))
               GEXTRA(3*(MJ3-1)+2)=GEXTRA(3*(MJ3-1)+2)-DUMMY2*(X(3*(MJ3-1)+2)-X(3*(MJ2-1)+2))
               GEXTRA(3*(MJ3-1)+3)=GEXTRA(3*(MJ3-1)+3)-DUMMY2*(X(3*(MJ3-1)+3)-X(3*(MJ2-1)+3))
            ENDDO
         ENDDO
      ENDDO
       
      DUMMYEPS=DUMMYEPS/(2.0D0*NCOPIES*LJADDREFNORM)
      GEXTRA(1:3*N)=GEXTRA(1:3*N)*TANHFAC/(2.0D0*NCOPIES*LJADDREFNORM)
!     PRINT '(A,G20.10)','eps att NN scale factor=',DUMMYEPS
      LJADDATTLOCAL(1:NADDTARGET,1:NADDTARGET)=LJADDATT(1:NADDTARGET,1:NADDTARGET)
!!!! debug
!     DUMMYEPS=1.0D0
!     GEXTRA(1:3*N)=0.0D0
!!!! debug
      DUMMY2=0.0D0
      DO J1=1,NUMNN
         K1=LJADDNN(J1,1)
         K2=LJADDNN(J1,2)
         DO J2=1,NCOPIES
            MJ2=K1+(J2-1)*NADDTARGET
            DO J3=1,NCOPIES
               MJ3=K2+(J3-1)*NADDTARGET
               LJADDATTLOCAL(K1,K2)=LJADDATT(K1,K2)*DUMMYEPS
               LJADDATTLOCAL(K2,K1)=LJADDATT(K2,K1)*DUMMYEPS
               DIST=(X(3*(MJ2-1)+1)-X(3*(MJ3-1)+1))**2 + 
     &              (X(3*(MJ2-1)+2)-X(3*(MJ3-1)+2))**2 +
     &              (X(3*(MJ2-1)+3)-X(3*(MJ3-1)+3))**2
               DIST=DIST**3
!
! We are assuming a constant NN att factor here, indepdendent of particle ids.
!
               DUMMY2=DUMMY2-4.0D0*LJADDATT(K1,K2)/DIST
            ENDDO
         ENDDO
      ENDDO
      GEXTRA(1:3*N)=GEXTRA(1:3*N)*DUMMY2
      ENERGY=0.0D0
      IF (GTEST.AND.(.NOT.STEST)) THEN
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            J3=3*J1
            XG(J1,J1)=0.0D0
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATTLOCAL(MJ2,MJ1))
               ENERGY=ENERGY+DUMMY
               DIST=DIST*R6
               XG(J2,J1)=-24.0D0*(2.0D0*R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATTLOCAL(MJ2,MJ1))*DIST
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO
      ELSEIF (GTEST) THEN
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            XG(J1,J1)=0.0D0
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               DUMMY=R6*(R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATTLOCAL(MJ2,MJ1))
               ENERGY=ENERGY+DUMMY
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               XG(J2,J1)=-24.0D0*(2.0D0*R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATTLOCAL(MJ2,MJ1))*R2(J1,J2)*R6
               XG(J1,J2)=XG(J2,J1)
            ENDDO
         ENDDO 
      ELSE
         DO J1=1,N
            MJ1=MOD(J1-1,NADDTARGET)+1
            J3=3*(J1-1)
            DO J2=J1+1,N
               MJ2=MOD(J2-1,NADDTARGET)+1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6=R2T**3
               ENERGY=ENERGY+R6*(R6*LJADDREP(MJ2,MJ1)-1.0D0*LJADDATTLOCAL(MJ2,MJ1))
            ENDDO
         ENDDO

      ENDIF
      ENERGY=4.0D0*ENERGY

      IF (.NOT.GTEST) RETURN
      DO J1=1,N
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,N
            J2=3*J4
            XMUL2=XG(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX+GEXTRA(J3-2)
         V(J3-1)=DUMMYY+GEXTRA(J3-1)
         V(J3)=DUMMYZ+GEXTRA(J3)
      ENDDO

      IF (.NOT.STEST) RETURN
      PRINT *,'ERROR *** second derivatives not yet coded'
      STOP
      CALL LJADDS4(G,F,R2,R14,R8,X,N,LJADDATTLOCAL)
!
! Now we need to add the extra terms.
!

      RETURN
      END
C*****************************************************************************

      SUBROUTINE LJADDS4(G,F,R2,R14,R8,X,N,LJADDATTLOCAL)
      USE MODHESS
      USE KEY, ONLY : LJADDREP, LJADDATT, NADDTARGET
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, MJ1, MJ2
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), LJADDATTLOCAL(N,N),
     2                 X(3*N),DUMMY

C
C  Calculate the g tensor.
C
      DO J1=1,N
         MJ1=MOD(J1-1,NADDTARGET)+1
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            MJ2=MOD(J2-1,NADDTARGET)+1
               G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)*LJADDREP(J2,J1)-R8(J2,J1)*LJADDATTLOCAL(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,N
         MJ1=MOD(J1-1,NADDTARGET)+1
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            MJ2=MOD(J2-1,NADDTARGET)+1
            F(J2,J1)=LJADDREP(MJ2,MJ1)*672.0D0*R14(J2,J1)-LJADDATTLOCAL(MJ2,MJ1)*192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJADDS3(G,F,R2,R14,R8,X,N)
      USE MODHESS
      USE KEY, ONLY : LJADDREP, LJADDATT, NADDTARGET
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, MJ1, MJ2
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), 
     2                 X(3*N),DUMMY

C
C  Calculate the g tensor.
C
      DO J1=1,N
         MJ1=MOD(J1-1,NADDTARGET)+1
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            MJ2=MOD(J2-1,NADDTARGET)+1
               G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)*LJADDREP(J2,J1)-R8(J2,J1)*LJADDATT(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,N
         MJ1=MOD(J1-1,NADDTARGET)+1
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            MJ2=MOD(J2-1,NADDTARGET)+1
            F(J2,J1)=LJADDREP(MJ2,MJ1)*672.0D0*R14(J2,J1)-LJADDATT(MJ2,MJ1)*192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END


