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
C  Subroutine LJPDIFF calculates the energy, cartesian gradient and second
C  derivative matrix analytically for pure Lennard-Jones in reduced units
C  (epsilon=sigma=1). If a negative cutoff is specified, each atom is
C  allowed to interact with the closest periodic image of the remaining
C  N-1 atoms regardless of separation (minimum image convention).
C
C  MM 9.ix.96
C
C*************************************************************************
C
      SUBROUTINE LJPDIFF(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST,PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3,
     1                 V(3*N), R2(N,N),  R6, R2T,
     2                 R8(N,N), G(N,N),
     3                 R14(N,N), F(N,N), 
     4                 POTEL, BOXLX, BOXLY, BOXLZ,
     5                 VEC(N,N,3), RCUT, IRCUT2, CUTOFF
      LOGICAL GTEST, STEST, PTEST
C
C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
      IF (CUTOFF.GT.0.D0) THEN
         IF ((BOXLX.LT.BOXLY).AND.(BOXLX.LT.BOXLZ)) THEN
           RCUT=BOXLX*CUTOFF
         ELSE IF (BOXLY.LT.BOXLZ) THEN
           RCUT=BOXLY*CUTOFF
         ELSE
           RCUT=BOXLZ*CUTOFF
         ENDIF
         IF (PTEST) PRINT*,'Cutoff used = ',RCUT
         IRCUT2 = 1.D0/RCUT**2
      ELSE
         IF (PTEST) PRINT*,'Minimum image convention with no cutoff.'
         IRCUT2=-1.0D0
         RCUT=-1.0D0
      ENDIF
C
C  Deal with any atoms that have left the box.
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
C
C  Calculate interatomic vectors using the minimum image convention.
C  VEC(i,j,alpha) is the alpha (x, y or z) component of the vector between
C  atoms i and j.
C
      POTEL=0.0D0
      IF (GTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
C 
C  Store distance matrices (unit of distance is sigma).
C
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               IF (R2(J2,J1).GT.IRCUT2) POTEL=POTEL + R6*(R6-1.0D0)
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1, N
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2T=VEC1**2+VEC2**2+VEC3**2
               R6=0.0D0
               IF (1.0D0/R2T.GT.IRCUT2) THEN
                  R2T=1.0D0/R2T
                  R6=R2T**3
                  POTEL=POTEL + R6*(R6-1.0D0)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF (.NOT.GTEST) RETURN
      CALL LJPG(G,R2,R14,R8,V,VEC,N,IRCUT2)
      
      IF (.NOT.STEST) RETURN
      CALL LJPS(G,F,R2,R14,R8,VEC,N,IRCUT2)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPG(G,R2,R14,R8,V,VEC,N,IRCUT2)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     2                 IRCUT2
C
C  Calculate the g tensor.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N 
            G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (R2(J4,J1).GT.IRCUT2) THEN
                  V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPS(G,F,R2,R14,R8,VEC,N,IRCUT2)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 VEC(N,N,3), R2(N,N),
     2                 IRCUT2, F(N,N)

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (R2(J4,J1).GT.IRCUT2) THEN
                  HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*
     1                 VEC(J4,J1,J2)**2 + G(J4,J1)   
               ENDIF
            ENDDO
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
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (R2(J5,J1).GT.IRCUT2) THEN
                     HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3) + 
     1               F(J5,J1)*R2(J5,J1)* 
     2               VEC(J5,J1,J2)*VEC(J5,J1,J4)
                  ENDIF
               ENDDO
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
               IF (R2(J4,J1).GT.IRCUT2) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           VEC(J4,J1,J2)**2-G(J4,J1)
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
               ENDIF
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
                  IF (R2(J4,J1).GT.IRCUT2) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *VEC(J4,J1,J2)*VEC(J4,J1,J5)
                  ELSE
                     HESS(J6,J3)=0.0D0
                  ENDIF
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
