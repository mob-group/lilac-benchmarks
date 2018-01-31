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
C  Subroutine DOUBLE adds a double-well potential between the
C  first two atoms.
C
      SUBROUTINE DOUBLE(N, X, V, EDOUBLE, GTEST, STEST, CUTOFF, H, W)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), VD(6), R2, EDOUBLE,
     1                 V(3*N), R, G(2,2), H, W, CUTOFF, AD(6,6)

      R2=(X(1)-X(4))**2+(X(2)-X(5))**2+(X(3)-X(6))**2
      R=SQRT(R2)
      R2=1.0D0/R2
      EDOUBLE=H*(1.0D0-((R-CUTOFF-W)/W)**2)**2

      IF (.NOT.GTEST) RETURN
      CALL DG(G,R,VD,X,H,W,CUTOFF)
      DO J1=1,6
         V(J1)=V(J1)+VD(J1)
      ENDDO
      
      IF (.NOT.STEST) RETURN
      CALL DS(G,R,R2,AD,X,H,W,CUTOFF)
      DO J1=1,6
         DO J2=1,6
            HESS(J2,J1)=HESS(J2,J1)+AD(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE DG(G,R,V,X,H,W,RWCA)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION G(2,2), V(6), X(*), DUMMY, R, H, W, RWCA
C
C  Calculate the g tensor.
C
      G(1,1)=0.0D0
      G(2,2)=0.0D0
      G(1,2)=4*h*(r - rwca)*(r - rwca - 2*w)*(r - rwca - w)/(r*w**4)
      G(2,1)=G(1,2)
C     WRITE(*,'(A,F20.10)') 'G(1,2)=',G(1,2)

      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,2
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE DS(G,R,R2,A,X,H,W,RWCA)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(2,2), R, A(6,6), X(*), DUMMY, H, W, RWCA, R2, F(2,2)

      F(1,1)=0.0D0
      F(2,2)=0.0D0
      F(2,1)=4*h*(2*r**3 - 3*r**2*(rwca + w) + rwca*(rwca + w)*(rwca + 2*w))/(r*w**4)
      F(1,2)=F(2,1)
C     WRITE(*,'(A,F20.10)') 'F(1,2)=',F(1,2)
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,2
               DUMMY=DUMMY+F(J4,J1)*R2*(X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            A(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,2
                  DUMMY=DUMMY + F(J5,J1)*R2*(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               A(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,2
               A(3*(J4-1)+J2,J3)=-F(J4,J1)*R2*(X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,2
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  A(J6,J3)=-F(J4,J1)*R2*(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
                  A(3*(J4-1)+J2,3*(J1-1)+J5)=A(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,6
         DO J2=J1+1,6
            A(J1,J2)=A(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END
