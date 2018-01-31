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
C  Here we calculate the two-body Morse potential 
C                                        
C*************************************************************************
C
      SUBROUTINE MORSE (N,X,P2,VNEW,RHO,GTEST,STEST)
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST,STEST
      DOUBLE PRECISION X(3*N), P2, RHO, TEMP, DIST, 
     1                 VNEW(3*N), R, RR(N,N), ERMR(N,N), ERMRM(N,N), ERMRT(N,N)

      P2=0.0D0
      DO J1=1,N
         RR(J1,J1)=0.0D0
         ERMR(J1,J1)=0.0D0
         ERMRM(J1,J1)=0.0D0
         ERMRT(J1,J1)=0.0D0
         J3=3*(J1-1)
         DO J2=J1+1,N
            J4=3*(J2-1)
            DIST=(X(J4+1)-X(J3+1))**2+(X(J4+2)-X(J3+2))**2+(X(J4+3)-X(J3+3))**2
            DIST=DSQRT(DIST)
            TEMP=DEXP(-RHO*(DIST-1.0D0))
            P2=P2+TEMP*(TEMP-2.0)
            IF (GTEST.OR.STEST) THEN
               R=RHO*DIST
               RR(J2,J1)=1.0D0/R
               ERMR(J2,J1)=TEMP
               ERMRM(J2,J1)=TEMP-1.0D0
               ERMRT(J2,J1)=2.0D0*TEMP-1.0D0
               RR(J1,J2)=RR(J2,J1)
               ERMR(J1,J2)=ERMR(J2,J1)
               ERMRM(J1,J2)=ERMRM(J2,J1)
               ERMRT(J1,J2)=ERMRT(J2,J1)
            ENDIF
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
      CALL MG(N, X, VNEW, RHO, RR, ERMR, ERMRM) 

      IF (.NOT.STEST) RETURN
      CALL MS(N, X, RHO, RR, ERMR, ERMRM, ERMRT) 

      RETURN
      END
C
C*************************************************************************
C
C  Subroutine MG calculates the cartesian gradient analytically for the Morse potential.
C
C*************************************************************************
C
      SUBROUTINE MG(N, X, V, RHO, RR, ERMR, ERMRM) 
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION X(3*N), RHO, V(3*N), RR(N,N),
     1                 DUMMY1, ERMR(N,N), ERMRM(N,N)
C
C  First calculate the gradient analytically.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY1=0.0D0
            DO J4=1,N
               DUMMY1=DUMMY1-2.0*RHO*RHO*(X(J3)-X(3*(J4-1)+J2))*
     1               ERMR(J4,J1)*ERMRM(J4,J1)*RR(J4,J1)
            ENDDO
            V(J3)=DUMMY1
         ENDDO
      ENDDO

      RETURN
      END
C
C*************************************************************************
C
C  Subroutine MS calculates the cartesian second derivative matrix 
C  analytically for the Morse potential.
C
C*************************************************************************
C
      SUBROUTINE MS(N, X, RHO, RR, ERMR, ERMRM, ERMRT) 
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6 
      DOUBLE PRECISION X(3*N), RHO, ERMR(N,N),
     1                 RR(N,N), DUMMY1, ERMRM(N,N),
     2                 ERMRT(N,N)
C
C  First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY1=0.0D0
            DO J4=1,N
               DUMMY1=DUMMY1+2.0*RHO*RHO*ERMR(J4,J1)
     1         *(  ERMRM(J4,J1) 
     2         * ( (RHO * (X(J3)-X(3*(J4-1)+J2) )*RR(J4,J1) )**2-1.0D0)  
     4         + (  RHO*(X(J3)-X(3*(J4-1)+J2)) )**2*RR(J4,J1)*
     5              ERMRT(J4,J1) ) * RR(J4,J1)
            ENDDO
            HESS(J3,J3)=DUMMY1
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
               DUMMY1=0.0D0
               DO J5=1,N
                 DUMMY1=DUMMY1
     1           +2.0*RHO**4*ERMR(J5,J1)
     2           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     3           *( ERMRM(J5,J1)*RR(J5,J1)+
     4           ERMRT(J5,J1)) * RR(J5,J1)**2
               ENDDO
               HESS(J3,3*(J1-1)+J4)=DUMMY1
               HESS(3*(J1-1)+J4,J3)=DUMMY1
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
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*ERMR(J4,J1)
     1      *(  (-RHO**2*(X(J3)-X(3*(J4-1)+J2))**2*RR(J4,J1)**2+1.0) 
     2           *   ERMRM(J4,J1)
     3           -   (X(J3)-X(3*(J4-1)+J2))**2 * RR(J4,J1) * RHO**2
     4           *   ERMRT(J4,J1)   )  * RR(J4,J1)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
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
                  HESS(J3,J6)=
     1           -2.0*RHO**4*ERMR(J4,J1)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( ERMRM(J4,J1)*RR(J4,J1)+
     4           ERMRT(J4,J1)) * RR(J4,J1)**2

                  HESS(J6,J3)=HESS(J3,J6)
               ENDDO
               DO J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1           -2.0*RHO**4*ERMR(J4,J1)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( ERMRM(J4,J1)*RR(J4,J1)+
     4           ERMRT(J4,J1)) * RR(J4,J1)**2
                  HESS(J6,J3)=HESS(J3,J6)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
