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
C  Subroutine P4DIFF calculates the cartesian gradient and second
C  derivative matrix analytically. Reduced units.
C
C*************************************************************************
C
      SUBROUTINE P4DIFF(N, X, V, ENERGY, K, GTEST, STEST)
      IMPLICIT NONE
      INTEGER N, L, J
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(N), ENERGY, K, DUM, V(N)

C     K=0.0 ! now assigned via PARAM1 allowing us to set it in odata
      L=SQRT(N*1.0D0)
C 
C  Calculate energy.
C
      ENERGY=0.0D0
      DO J=1,N
         DUM=0.5*(X(J))**2
         IF (J+L.GT.N) THEN
            IF (J+1.GT.N) THEN
               ENERGY=ENERGY+DUM*(DUM-1.0)-K*X(J)*(X(J+1-N)+X(J+L-N))
            ELSE 
               ENERGY=ENERGY+DUM*(DUM-1.0)-K*X(J)*(X(J+1)+X(J+L-N))
            ENDIF
         ELSE
            ENERGY=ENERGY+DUM*(DUM-1.0)-K*X(J)*(X(J+1)+X(J+L))
         ENDIF
      ENDDO

      IF (.NOT.GTEST) RETURN
      CALL P4G(V,X,N,L,K)
      
      IF (.NOT.STEST) RETURN
      CALL P4S(X,N,L,K)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE P4G(V,X,N,L,K)
      IMPLICIT NONE
      INTEGER N, L, J
      DOUBLE PRECISION V(N), X(N), K

      DO J=1,N
         V(J)=X(J)*(X(J)*X(J)-1.0)
         IF (J+L .GT. N) THEN
            IF (J+1 .GT. N) THEN
               V(J)=V(J)-K*(X(J+1-N)+X(J+L-N)+X(J-1)+X(J-L))
            ELSE
               V(J)=V(J)-K*(X(J+1)+X(J+L-N)+X(J-1)+X(J-L))
            ENDIF
         ELSE IF (J-L .LT. 1) THEN
            IF (J-1 .LT. 1) THEN
               V(J)=V(J)-K*(X(J+1)+X(J+L)+X(J-1+N)+X(J-L+N))
            ELSE
               V(J)=V(J)-K*(X(J+1)+X(J+L)+X(J-1)+X(J-L+N))
            ENDIF
         ELSE
            V(J)=V(J)-K*(X(J+1)+X(J+L)+X(J-1)+X(J-L))
         ENDIF
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE P4S(X,N,L,K)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, L, J1, J2
      DOUBLE PRECISION X(N), K
C
C  Now do the hessian.
C
      HESS(1:N,1:N)=0.0D0
      DO J1=1,N
         HESS(J1,J1)=3.0*X(J1)*X(J1)-1.0
      ENDDO

      HESS(N,1)=-K
      HESS(1,N)=-K
      DO J1=1,N-1
         HESS(J1,J1+1)=-K
         HESS(J1+1,J1)=-K
      ENDDO

      DO J1=1,N-L
         HESS(J1,J1+L)=-K
         HESS(J1+L,J1)=-K
      ENDDO

      DO J1=N-L+1,N
         HESS(J1,J1-N+L)=-K
         HESS(J1-N+L,J1)=-K
      ENDDO

      RETURN
      END
