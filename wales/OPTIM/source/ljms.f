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
c
C*************************************************************************
C
C  Subroutine LJMS (LJ-MS potentials) calculates the cartesian 
C  gradient and second
C  derivative matrix analytically. Reduced units, epsilon and rm.
C
C*************************************************************************
C
      SUBROUTINE LJMS (N, EPSILON, RM, GAMMA, X, V, ENERGY, GTEST, STEST)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, 
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N), ENERGY1,
     3                 ENERGY2, R14(N,N), F(N,N), 
     4                 GAMMA, R6(N,N), R4T, R6T,
     5                 epsilon,rm,rm4,rm6,rm12
C 
C  Store distance matrices.
C
      ENERGY1=0.0D0
      ENERGY2=0.0D0
C     epsilon=12.d0
C     rm=0.78d0
      rm12=rm**12
      rm6=rm**6
      rm4=rm**4

      IF (GTEST) THEN
         DO J1=1,N-1
            R2(J1,J1)=0.0D0
            R6(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N-1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6(J2,J1)=R2(J2,J1)**3
               ENERGY1=ENERGY1+R6(J2,J1)*(R6(J2,J1)-2.0D0)
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
            ENDDO
         ENDDO
         J1=N
            R2(J1,J1)=0.0D0
            R6(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=1,N-1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R2(J2,J1)=R2(J2,J1)
               R4T=R2(J2,J1)**2
               R6(J2,J1)=R2(J2,J1)**3
               ENERGY2=ENERGY2+(1.d0+GAMMA)*rm12*R6(J2,J1)*R6(J2,J1)
     1                -4.d0*GAMMA*rm6*R6(J2,J1)
     2                -3.d0*(1.d0-GAMMA)*rm4*R4T
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
            ENDDO
      ELSE
         DO J1=1,N-1
            J3=3*(J1-1)
            DO J2=J1+1,N-1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6T=R2T**3
               ENERGY1=ENERGY1+R6T*(R6T-2.0D0)
            ENDDO
         ENDDO
         J1=N
            J3=3*(J1-1)
            DO J2=1,N-1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R2T=R2T
               R4T=R2T**2
               R6T=R2T**3
               ENERGY2=ENERGY2+(1.d0+GAMMA)*rm12*R6T*R6T
     1                -4.d0*GAMMA*rm6*R6T
     2                -3.d0*(1.d0-GAMMA)*rm4*R4T
            ENDDO

      ENDIF
      ENERGY=ENERGY1+epsilon*0.5d0*ENERGY2

      IF (.NOT.GTEST) RETURN
      CALL LJMSG(G,R14,R8,R6,V,X,N,GAMMA,rm12,rm6,rm4,epsilon)
      
      IF (.NOT.STEST) RETURN
      CALL LJMSS(G,F,R2,R14,R8,R6,X,N,GAMMA,rm12,rm6,rm4,epsilon)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJMSG(G,R14,R8,R6,V,X,N,GAMMA,rm12,rm6,rm4,epsilon)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 V(3*N), X(3*N), DUMMY, R6(N,N), 
     2                 epsilon,GAMMA,rm12,rm6,rm4
C
C  Calculate the g tensor.
C
      DO J1=1,N-1
         G(J1,J1)=0.0D0
         DO J2=J1+1,N-1
            G(J2,J1)=-12.0D0*(R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      J1=N
         G(J1,J1)=0.0D0
         DO J2=1,N-1
            G(J2,J1)=-6.d0*((1.d0+GAMMA)*rm12*R14(J2,J1)
     1               -2.d0*GAMMA*rm6*R8(J2,J1)
     2               -(1.d0-GAMMA)*rm4*R6(J2,J1))       
            G(J2,J1)=epsilon*G(J2,J1)
            G(J1,J2)=G(J2,J1)
         ENDDO
C
C  From here on down the code is system-independent!
C  First calculate the gradient analytically.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJMSS(G,F,R2,R14,R8,R6,X,N,GAMMA,rm12,rm6,rm4,epsilon)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), X(3*N),
     3                 DUMMY,GAMMA,rm12,rm6,rm4,epsilon,R6(N,N)

      DO J1=1,N-1
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=168.0D0*R14(J2,J1)-96.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO

      J1=N
         F(J1,J1)=0.0D0
         DO J2=1,N-1
           F(J2,J1)=84.0D0*(1.d0+GAMMA)*rm12*R14(J2,J1)-16.0D0*GAMMA*rm6*
     1              R8(J2,J1)-36.d0*(1.d0-GAMMA)*rm4*R6(J2,J1)
           F(J2,J1)=epsilon*F(J2,J1)
           F(J1,J2)=F(J2,J1)
         END DO

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
