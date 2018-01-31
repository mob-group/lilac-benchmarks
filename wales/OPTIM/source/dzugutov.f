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
C-----------------------------------------------------------------------* 
C 
C  Energy and gradient for the Dzugutov potential.  
C 
C  16.02.2000
C  Sergei Simdyankin <ssim@nada.kth.se>
C
C correspondence with the parameters in the file "potential.f"
C V     = GRAD
C EDZ   = EREAL
C GTEST = GRADT
C STEST 
C
C meaning of the variables
C DIST - distance
C EDZ  - total interaction energy
C G    - (N,N) storage for gradient intermediate bits
C V    - (3*N) vector - gradient (force acting on a particle)
! 
C-----------------------------------------------------------------------*
      SUBROUTINE DZUGUTOV(NATOMS,X,V,EDZ,GTEST,STEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST, STEST
      INTEGER J1, J2, J3, J4, NATOMS, J5, J6
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY, F(NATOMS,NATOMS), 
     2                 DDUMMY, DDDUMMY,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6, R2(NATOMS,NATOMS),PARAM7
      
      EDZ=0.0D0
      IF (GTEST.OR.STEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            F(J1,J1)=0.0D0
            R2(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
                     DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               call derphi(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
               EDZ=EDZ+DUMMY
               R2(J2,J1)=1.0D0/DIST
               R2(J1,J2)=R2(J2,J1)
               G(J2,J1)=DDUMMY
               G(J1,J2)=G(J2,J1)
               IF (STEST) THEN
                  F(J2,J1)=DDDUMMY
                  F(J1,J2)=DDDUMMY
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               call derphi(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
               EDZ=EDZ+DUMMY
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.(GTEST.OR.STEST)) RETURN

      DO J1=1,NATOMS
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,NATOMS
            J2=3*J4
            XMUL2=G(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,NATOMS
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
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,NATOMS
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
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,NATOMS
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,NATOMS
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
      DO J1=1,3*NATOMS
         DO J2=J1+1,3*NATOMS
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END

C%%%%%%%
C values of the potential, first, and second derivatives
C if secder = .false. ddphi = 0 on output
C%%%%%%%
C-----------------------------------------------------------------------*
      subroutine derphi(r2,phi,dphi,ddphi,secder,m,A,c,aa,B,d,bb)
      implicit none
      DOUBLE PRECISION R2, PHI, DPHI, DDPHI
      logical secder
      integer m2

      DOUBLE PRECISION A, AA, B, BB, C, D, M
C     parameter(m=16, A=5.82d0, c=1.10d0, aa=1.87d0, B=1.28d0,  d=0.27d0, bb=1.94d0)
C     parameter(m=4,  A=3.00d0, c=0.52d0, aa=1.65d0, B=2.109d0, d=0.55d0, bb=1.94d0)
      DOUBLE PRECISION R, V1, DV1, DDV1, V2, DV2, DDV2, DUMMY
      DOUBLE PRECISION EXPO, EXPARG

      m2=m/2
      r = sqrt(r2)
C
C  Initialization needed for Sun!
C
      V1=0.0D0
      V2=0.0D0
      dV1=0.0D0
      ddV1=0.0D0
      dV2=0.0D0
      ddV2=0.0D0

      IF (R .LT. BB) THEN
         EXPARG = D/(R-BB)
C        IF (DABS(EXPARG).GT.708.0D0) THEN 
C           EXPO = 0.0D0
C        ELSE
            EXPO = DEXP(EXPARG)
C        ENDIF
         V2 =  B*EXPO
         DV2 = -V2*D/(R-BB)**2
         IF (SECDER) DDV2 = -DV2*D/(R-BB)**2 + 2*V2*D/(R-BB)**3
         IF (R.LT.AA) THEN 
            EXPARG = C/(R-AA)
C           IF (DABS(EXPARG).GT.708.0D0) THEN 
C              EXPO=0.0D0
C           ELSE
               EXPO = DEXP(EXPARG)
C           ENDIF
            DUMMY=1.0D0/R2**M2
            V1=A*(DUMMY-B)*EXPO
            DV1=-A*M*(DUMMY/R)*EXPO - V1*C/(R-AA)**2
            IF (SECDER) THEN
               DDV1 = A*M*(M+1)*(R**(-M-2))*EXPO 
     1             +  A*M*(R**(-M-1))*EXPO*C/((R-AA)**2) 
     2             -  DV1*C/((R-AA)**2) + 2*V1*C/((R-AA)**3) 
            ENDIF
         ENDIF
      ENDIF
  
      PHI=V1+V2
      DPHI=(DV1+DV2)/R
      IF (SECDER) DDPHI = DDV1 + DDV2

      RETURN
      end                       
