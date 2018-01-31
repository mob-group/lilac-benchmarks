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
C  Subroutine SW2 calculates the energy, Cartesian gradient and second
C  derivative matrix analytically for the SW Si potential.
C
C*************************************************************************
C
      SUBROUTINE SWTWO(N, X, V, ENERGY, E3, BOXLX, BOXLY, BOXLZ, GTEST, STEST, LAMBDA)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, ANV(N,N,3)
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, 
     1                 V(3*N), R(N,N), R4T, VEC(N,N,3),
     2                 G(N,N), BOXLX, BOXLY, BOXLZ, LAMBDA,
     3                 ER(N,N), F(N,N), R2(N,N), E3
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)
C     COMMON /ANVCOMMON/ ANV

C
C  Deal with any atoms that have left the box.
C
      IF (BULKT.AND.(.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
C 
C  Store distance matrices.
C
      ENERGY=0.0D0
C     IF (BULKT.AND.(.NOT.FIXIMAGE)) PRINT*,'recalculating ANV'
C     IF (GTEST) THEN
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
               IF (BULKT.AND.(.NOT.FIXIMAGE)) THEN
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
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R(J1,J1)=0.0D0
            ER(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=VEC(J2,J1,1)**2+VEC(J2,J1,2)**2+VEC(J2,J1,3)**2
               R(J2,J1)=SQRT(R2(J2,J1)) 
               IF (R(J2,J1).LT.SMALLA) THEN
                  R2(J2,J1)=1.0D0/R2(J2,J1)
                  R4T=R2(J2,J1)**2 
                  ER(J2,J1)=EXP(1.0D0/(-SMALLA + R(J2,J1)))
                  ENERGY=ENERGY+ER(J2,J1)*(-1.0D0 + BIGB*R4T) 
               ENDIF
               R(J1,J2)=R(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               ER(J1,J2)=ER(J2,J1)
            ENDDO 
         ENDDO
C
C  SWthree needs to know VEC
C
C     ELSE
C        DO J1=1,N
C           J3=3*(J1-1)
C           DO J2=J1+1,N
C              J4=3*(J2-1)
C              VEC1=X(J3+1)-X(J4+1)
C              VEC2=X(J3+2)-X(J4+2)
C              VEC3=X(J3+3)-X(J4+3)
C              IF (BULKT.AND.(.NOT.FIXIMAGE)) THEN
C                 ANV(J2,J1,1)=NINT(VEC1/BOXLX)
C                 ANV(J2,J1,2)=NINT(VEC2/BOXLY)
C                 ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
C              ENDIF
C              VEC1=VEC1-BOXLX*ANV(J2,J1,1)
C              VEC2=VEC2-BOXLY*ANV(J2,J1,2)
C              VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
C              R2T=VEC1**2+VEC2**2+VEC3**2
C              DIST=SQRT(R2T)
C              R4T=R2T**2
C              IF (DIST.LT.SMALLA) ENERGY=ENERGY+EXP(1.0D0/(-SMALLA + DIST))*(-1.0D0 + BIGB/R4T)
C           ENDDO
C        ENDDO
C     ENDIF

      ENERGY=BIGA*ENERGY

      IF (GTEST) CALL SWTWOG(G,R,R2,ER,V,X,N,VEC)
      
      IF (STEST) CALL SWTWOS(G,F,R,R2,ER,X,N,VEC)

      CALL SWTHREE(N, X, V, E3, GTEST, STEST, LAMBDA, R, R2, G, F, VEC)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE SWTWOG(G,R,R2,ER,V,X,N,VEC)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R(N,N), R2(N,N),
     1                 V(3*N), X(3*N), DUMMY, ER(N,N), VEC(N,N,3)
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)
C
C  Calculate the g tensor.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (R(J2,J1).LT.SMALLA) THEN
               G(J2,J1)=(BIGA*ER(J2,J1)*(-4*SMALLA**2*BIGB + 
     1                    (-1.0D0 + 8.0D0*SMALLA)*BIGB*R(J2,J1) - 4*BIGB/R2(J2,J1) + R(J2,J1)**5))*
     2                    R2(J2,J1)**3/(SMALLA - R(J2,J1))**2
            ELSE
               G(J2,J1)=0.0D0
            ENDIF
            G(J1,J2)=G(J2,J1)
         ENDDO
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
               IF (R(J4,J1).LT.SMALLA) DUMMY=DUMMY+G(J4,J1)*VEC(J4,J1,J2)
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE SWTWOS(G,F,R,R2,ER,X,N,VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R(N,N), VEC(N,N,3),
     1                 F(N,N), 
     2                 X(3*N),DUMMY, R2(N,N), ER(N,N)
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            IF (R(J2,J1).LT.SMALLA) THEN
               F(J2,J1)=BIGA*ER(J2,J1)*((-1.0D0 + 2*SMALLA - 2*R(J2,J1))/(SMALLA - R(J2,J1))**4 + 
     1                  R2(J2,J1)**3*(BIGB*(20.0D0 + (R(J2,J1)*(2.0D0*(4.0D0*SMALLA - 5.0D0*R(J2,J1))*
     1                  (SMALLA - R(J2,J1)) + R(J2,J1)))/(SMALLA - R(J2,J1))**4))) - G(J2,J1)
            ELSE
               F(J2,J1)=0.0D0
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
               IF (R(J4,J1).LT.SMALLA) DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)   
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
                  IF (R(J5,J1).LT.SMALLA) DUMMY=DUMMY + 
     1                F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
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
               IF (R(J4,J1).LT.SMALLA) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1) 
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
                  IF (R(J4,J1).LT.SMALLA) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
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

      SUBROUTINE SWLATMIN(N,XC,BOXLX,BOXLY,BOXLZ,V,LAMBDA)
      IMPLICIT NONE
      INTEGER N, NCOUNT, J1
      DOUBLE PRECISION XC(3*N),F1,F2,F3,GRAD,SECOND,XSAVE(3*N),
     1                 P2, P3,V(3*N),LAMBDA,BOXLX,BOXLY,BOXLZ,DIF,TEMP1
C
C  Value of DIF is the order of magnitude to which the lattice
C  constant can be optimised. Setting it smaller than 10^(-7)
C  causes numerical problems on the DEC.
C
      DIF=1.0D-7
      NCOUNT=1
      DO 20 J1=1,3*N
         XSAVE(J1)=XC(J1)
20    CONTINUE

10    TEMP1=BOXLX
      BOXLX=TEMP1+DIF
      BOXLY=TEMP1+DIF
      BOXLZ=TEMP1+DIF

      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F1=P2+P3

      BOXLX=TEMP1-DIF
      BOXLY=TEMP1-DIF
      BOXLZ=TEMP1-DIF
      DO 30 J1=1,3*N
         XC(J1)=XSAVE(J1)
30    CONTINUE
      CALL SWTWO  (N, XC,  V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC,  V, P3, .FALSE., .FALSE.,LAMBDA)
      F2=P2+P3

      GRAD=(F1-F2)/(2.0D0*DIF)

      BOXLX=TEMP1
      BOXLY=TEMP1
      BOXLZ=TEMP1
      DO 40 J1=1,3*N
         XC(J1)=XSAVE(J1)
40    CONTINUE
      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F1=P2+P3

      BOXLX=TEMP1+2.0D0*DIF
      BOXLY=TEMP1+2.0D0*DIF
      BOXLZ=TEMP1+2.0D0*DIF
      DO 50 J1=1,3*N
         XC(J1)=XSAVE(J1)
50    CONTINUE
      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F2=P2+P3

      BOXLX=TEMP1-2.0D0*DIF
      BOXLY=TEMP1-2.0D0*DIF
      BOXLZ=TEMP1-2.0D0*DIF
      DO 60 J1=1,3*N
         XC(J1)=XSAVE(J1)
60    CONTINUE
      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F3=P2+P3

      SECOND=(F3+F2-2.0D0*F1)/(4.0D0*DIF*DIF)
      PRINT*,'Energy for lattice cycle ',NCOUNT,' is ',F1
      PRINT*,'Gradient wrt box length=',GRAD
      PRINT*,'Second derivative wrt box length=',SECOND
      NCOUNT=NCOUNT+1
      IF (DABS(GRAD/SECOND).GT.1.0D-4) THEN
         BOXLX=BOXLX-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLY=BOXLY-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLZ=BOXLZ-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         PRINT*,'Step=',-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         GOTO 10
      ELSE
         BOXLX=BOXLX-GRAD/SECOND
         BOXLY=BOXLY-GRAD/SECOND
         BOXLZ=BOXLZ-GRAD/SECOND
         PRINT*,'Step=',-GRAD/SECOND
         IF (DABS(GRAD/SECOND).GT.1.0D-6) GOTO 10
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  Subroutine SW2 calculates the energy, Cartesian gradient and second
C  derivative matrix analytically for the SW Si potential.
C
C*************************************************************************
C
      SUBROUTINE SWTHREE(N, X, V, ENERGY, GTEST, STEST, LAMBDA, R, R2, G, F, VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, K1, K2, K3
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, GAMMA2, ER2131,
     1                 V(3*N), R(N,N), R122, R123, R132, R133, R12S, R13S,
     2                 G(N,N), r12DOTr13PT, r12DOTr13, R134, R124, ERLDOT,
     3                 ER(N,N), F(N,N), R2(N,N), THIRD,
     4                 ER21, ER31, R12, R13, R23, r13s3, r12s3, DOT1,
     5                 r13s4, r12s4, DOT2, ERL, 
     6                 VEC(N,N,3), V211, V212, V213, V311, V312, V313, V2112, V2122, V2132,
     7                 V3112, V3122, V3132
      DOUBLE PRECISION LAMBDA, GAMMA, SMALLA
      PARAMETER (GAMMA=1.2D0, SMALLA=1.8D0, THIRD=0.33333333333333333333D0)
      
C 
C  Store things
C
      DO J1=1,N
         ER(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (R(J2,J1).LT.SMALLA) THEN
               ER(J2,J1)=EXP(GAMMA/(-SMALLA + R(J2,J1)))
            ELSE
               ER(J2,J1)=0.0D0
            ENDIF
            ER(J1,J2)=ER(J2,J1)
C           WRITE(*,'(A,2I4,2E20.10)') 'J1,J2,R,ER=',J1,J2,R(J2,J1),ER(J2,J1)
         ENDDO
      ENDDO

      ENERGY=0.0D0
      DO J1=1,N
         K1=3*(J1-1)
         DO J2=1,N
            IF ((R(J2,J1).LT.SMALLA).AND.(J2.NE.J1)) THEN
               K2=3*(J2-1)
               ER21=ER(J2,J1)
               R12=R(J2,J1)
               V211=VEC(J2,J1,1)
               V212=VEC(J2,J1,2)
               V213=VEC(J2,J1,3)
               V2112=V211**2
               V2122=V212**2
               V2132=V213**2
               DO J3=1,N
                  IF (((R(J3,J1).LT.SMALLA).AND.(J3.NE.J1)).AND.(J3.NE.J2)) THEN
                     K3=3*(J3-1)
                     ER31=ER(J3,J1)
                     ER2131=ER21*ER31
                     R13=R(J3,J1)
                     R23=R(J3,J2)
                     V311=VEC(J3,J1,1)
                     V312=VEC(J3,J1,2)
                     V313=VEC(J3,J1,3)
                     V3112=V311**2
                     V3122=V312**2
                     V3132=V313**2
                     r12DOTr13=(V211*V311+V212*V312+V213*V313)/(r12*r13)
                     DOT1=r12DOTR13**2
                     r12DOTr13PT=(r12DOTr13+third)
                     DOT2=r12DOTr13PT**2

                     ENERGY=ENERGY+ER2131*DOT2

                     IF (GTEST) THEN
                        R122=R12**2
                        R132=R13**2
                        r12s=(r12 - smalla)**2
                        r13s=(r13 - smalla)**2
                        ERLDOT=ER2131*lambda*r12DOTr13PT

                       V(K1+1)=V(K1+1)+ (ERLDOT*
     1    ((2*r13*(r12 - r12DOTr13*r13)*V211 + 
     1         2*r12*(-(r12*r12DOTr13) + r13)*V311)/(r122*r132) + 
     1      gamma*r12DOTr13PT*(-(V211/(r12*r12s)) - 
     1         V311/(r13*r13s))))/2.0D0
                        V(K1+2)=V(K1+2)+ (ERLDOT*
     1    ((2*r13*(r12 - r12DOTr13*r13)*V212 + 
     1         2*r12*(-(r12*r12DOTr13) + r13)*V312)/(r122*r132) + 
     1      gamma*r12DOTr13PT*(-(V212/(r12*r12s)) - 
     1         V312/(r13*r13s))))/2.0D0
                        V(K1+3)=V(K1+3)+ (ERLDOT*
     1    ((2*r13*(r12 - r12DOTr13*r13)*V213 + 
     1         2*r12*(-(r12*r12DOTr13) + r13)*V313)/(r122*r132) + 
     1      gamma*r12DOTr13PT*(-(V213/(r12*r12s)) - 
     1         V313/(r13*r13s))))/2.0D0
                        V(K2+1)=V(K2+1)+ (ERLDOT*
     1    ((2*r12DOTr13 + (gamma*r12*r12DOTr13PT)/r12s)*
     1       V211 - (2*r12*V311)/r13))/(2.*r122)
                        V(K2+2)=V(K2+2)+ (ERLDOT*
     1    ((2*r12DOTr13 + (gamma*r12*r12DOTr13PT)/r12s)*
     1       V212 - (2*r12*V312)/r13))/(2.*r122)
                        V(K2+3)=V(K2+3)+ (ERLDOT*
     1    ((2*r12DOTr13 + (gamma*r12*r12DOTr13PT)/r12s)*
     1       V213 - (2*r12*V313)/r13))/(2.*r122)
                        V(K3+1)=V(K3+1)+ (ERLDOT*
     1    ((-2*r13*V211)/r12 + 
     1      (2*r12DOTr13 + (gamma*r12DOTr13PT*r13)/r13s)*
     1       V311))/(2.*r132)
                       V(K3+2)=V(K3+2)+ (ERLDOT*
     1    ((-2*r13*V212)/r12 + 
     1      (2*r12DOTr13 + (gamma*r12DOTr13PT*r13)/r13s)*
     1       V312))/(2.*r132)
                        V(K3+3)=V(K3+3)+ (ERLDOT*
     1    ((-2*r13*V213)/r12 + 
     1      (2*r12DOTr13 + (gamma*r12DOTr13PT*r13)/r13s)*
     1       V313))/(2.*r132)

                     ENDIF
                     IF (STEST) THEN
                        ERL=ER2131*lambda
                        GAMMA2=GAMMA**2
                        R122=R12**2
                        R123=R12*R122
                        R124=R12*R123
                        R132=R13**2
                        R133=R13*R132
                        R134=R13*R133
                        r12s=(r12 - smalla)**2
                        r13s=(r13 - smalla)**2
                        r12s3=r12s*(r12 - smalla)
                        r13s3=r13s*(r13 - smalla)
                        r12s4=r12s3*(r12 - smalla)
                        r13s4=r13s3*(r13 - smalla)

                        HESS(K1+1,K1+1)=HESS(K1+1,K1+1)+ (ERL*((2*(r13*(-r12 + r12DOTr13*r13)*V211 + 
     1            r12*(r12*r12DOTr13 - r13)*V311)**2)/(R124*R134)
     1       + (4*gamma*r12DOTr13PT*
     1         (r13*(r12 - r12DOTr13*r13)*V211 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V311)*
     1         (-(V211/(r12*r12s)) - 
     1           V311/(r13*r13s)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (V211/(r12*r12s) + 
     1          V311/(r13*r13s))**2 + 
     1      gamma*DOT2*
     1       (-(1/(r12*r12s)) - 1/(r13*r13s) + 
     1         ((3*r12 - smalla)*V2112)/
     1          (R123*r12s3) + 
     1         ((3*r13 - smalla)*V3112)/(R133*r13s3))
     1       - (2*r12DOTr13PT*(R133*(2*r12 - 3*r12DOTr13*r13)*
     1            V2112 + 
     1           2*r12*r13*(R122 - r12*r12DOTr13*r13 + R132)*V211*
     1            V311 + 
     1           R122*(R132*
     1               (R122*r12DOTr13 - 2*r12*r13 + r12DOTr13*R132) + 
     1              r12*(-3*r12*r12DOTr13 + 2*r13)*V3112)))/
     1       (R124*R134)))/2.
                        HESS(K1+1,K1+2)=HESS(K1+1,K1+2)+ (ERL*((2*(r13*(-r12 + r12DOTr13*r13)*V211 + 
     1           r12*(r12*r12DOTr13 - r13)*V311)*
     1         (r13*(-r12 + r12DOTr13*r13)*V212 + 
     1           r12*(r12*r12DOTr13 - r13)*V312))/(R124*R134) + 
     1      (2*gamma*r12DOTr13PT*
     1         (-(V211/(r12*r12s)) - 
     1           V311/(r13*r13s))*
     1         (r13*(r12 - r12DOTr13*r13)*V212 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V312))/(R122*R132)
     1       + (2*gamma*r12DOTr13PT*
     1         (r13*(r12 - r12DOTr13*r13)*V211 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V311)*
     1         (-(V212/(r12*r12s)) - 
     1           V312/(r13*r13s)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (-(V211/(r12*r12s)) - 
     1         V311/(r13*r13s))*
     1       (-(V212/(r12*r12s)) - 
     1         V312/(r13*r13s)) + 
     1      gamma*DOT2*
     1       (((3*r12 - smalla)*V211*V212)/
     1          (R123*r12s3) + 
     1         ((3*r13 - smalla)*V311*V312)/
     1          (R133*r13s3)) + 
     1      (2*r12DOTr13PT*(r12*V311*
     1            (-(r13*(R122 - r12*r12DOTr13*r13 + R132)*
     1                 V212) + 
     1              R122*(3*r12*r12DOTr13 - 2*r13)*V312) - 
     1           r13*V211*
     1            (R132*(2*r12 - 3*r12DOTr13*r13)*V212 + 
     1              r12*(R122 - r12*r12DOTr13*r13 + R132)*V312)))/
     1       (R124*R134)))/2.
                        HESS(K1+1,K1+3)=HESS(K1+1,K1+3)+ (ERL*((2*(r13*(-r12 + r12DOTr13*r13)*V211 + 
     1           r12*(r12*r12DOTr13 - r13)*V311)*
     1         (r13*(-r12 + r12DOTr13*r13)*V213 + 
     1           r12*(r12*r12DOTr13 - r13)*V313))/(R124*R134) + 
     1      (2*gamma*r12DOTr13PT*
     1         (-(V211/(r12*r12s)) - 
     1           V311/(r13*r13s))*
     1         (r13*(r12 - r12DOTr13*r13)*V213 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V313))/(R122*R132)
     1       + (2*gamma*r12DOTr13PT*
     1         (r13*(r12 - r12DOTr13*r13)*V211 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V311)*
     1         (-(V213/(r12*r12s)) - 
     1           V313/(r13*r13s)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (-(V211/(r12*r12s)) - 
     1         V311/(r13*r13s))*
     1       (-(V213/(r12*r12s)) - 
     1         V313/(r13*r13s)) + 
     1      gamma*DOT2*
     1       (((3*r12 - smalla)*V211*V213)/
     1          (R123*r12s3) + 
     1         ((3*r13 - smalla)*V311*V313)/
     1          (R133*r13s3)) + 
     1      (2*r12DOTr13PT*(r12*V311*
     1            (-(r13*(R122 - r12*r12DOTr13*r13 + R132)*
     1                 V213) + 
     1              R122*(3*r12*r12DOTr13 - 2*r13)*V313) - 
     1           r13*V211*
     1            (R132*(2*r12 - 3*r12DOTr13*r13)*V213 + 
     1              r12*(R122 - r12*r12DOTr13*r13 + R132)*V313)))/
     1       (R124*R134)))/2.
                        HESS(K1+1,K2+1)=HESS(K1+1,K2+1)+ (ERL*((-2*DOT1 + 
     1         r12DOTr13*(-6*r12DOTr13PT + (2*r12)/r13 - 
     1            (4*gamma*r12*r12DOTr13PT)/r12s) + 
     1         r12*r12DOTr13PT*(2/r13 - 
     1            (GAMMA2*r12*r12DOTr13PT)/r12s4 - 
     1            (3*gamma*r12*r12DOTr13PT)/r12s3 + 
     1            (2*gamma*r12)/(r13*r12s) + 
     1            (gamma*r12DOTr13PT*smalla)/r12s3))*
     1       V2112 + (r12*
     1         (4*(r12DOTr13 + r12DOTr13PT)*r13 + 
     1           r12*(-2 - 2*DOT1 - 2*r12DOTr13*r12DOTr13PT + 
     1              (4*gamma*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s) + 
     1           (gamma*R122*r12DOTr13PT*
     1              (-2*r12DOTr13 - (gamma*r12DOTr13PT*r13)/r13s)
     1              )/r12s)*V211*V311)/R132 + 
     1      R122*(r12DOTr13PT*
     1          (2*r12DOTr13 - (2*r12)/r13 + 
     1            (gamma*R122*r12DOTr13PT)/r12s3 - 
     1            (gamma*r12*r12DOTr13PT*smalla)/r12s3) + 
     1         (2*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V3112)/R133)))/(2.*R124)
                        HESS(K1+1,K2+2)=HESS(K1+1,K2+2)+ (ERL*(V211*
     1       ((-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12*r12DOTr13)/r13 + (2*r12*r12DOTr13PT)/r13 - 
     1            (GAMMA2*R122*DOT2)/r12s4 - 
     1            (3*gamma*R122*DOT2)/r12s3 - 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*R122*r12DOTr13PT)/(r13*r12s) + 
     1            (gamma*r12*DOT2*smalla)/r12s3)*
     1          V212 + (2*r12*
     1            (-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V312)/R132) + 
     1      (r12*V311*(r13*
     1            (-2*r12*DOT1 - 2*r12*r12DOTr13*r12DOTr13PT + 
     1              2*r12DOTr13*r13 + 2*r12DOTr13PT*r13 - 
     1              (2*gamma*R122*r12DOTr13*r12DOTr13PT)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r13s - 
     1              (GAMMA2*R122*DOT2*r13)/
     1               (r12s*r13s))*V212 + 
     1           2*r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V312))/R133))/(2.*R124)
                        HESS(K1+1,K2+3)=HESS(K1+1,K2+3)+ (ERL*(V211*
     1       ((-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12*r12DOTr13)/r13 + (2*r12*r12DOTr13PT)/r13 - 
     1            (GAMMA2*R122*DOT2)/r12s4 - 
     1            (3*gamma*R122*DOT2)/r12s3 - 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*R122*r12DOTr13PT)/(r13*r12s) + 
     1            (gamma*r12*DOT2*smalla)/r12s3)*
     1          V213 + (2*r12*
     1            (-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V313)/R132) + 
     1      (r12*V311*(r13*
     1            (-2*r12*DOT1 - 2*r12*r12DOTr13*r12DOTr13PT + 
     1              2*r12DOTr13*r13 + 2*r12DOTr13PT*r13 - 
     1              (2*gamma*R122*r12DOTr13*r12DOTr13PT)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r13s - 
     1              (GAMMA2*R122*DOT2*r13)/
     1               (r12s*r13s))*V213 + 
     1           2*r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V313))/R133))/(2.*R124)
                        HESS(K1+1,K3+1)=HESS(K1+1,K3+1)+ (ERL*(r12DOTr13PT*R132*
     1       (2*r12DOTr13 - (2*r13)/r12 + 
     1         (gamma*r12DOTr13PT*R132)/r13s3 - 
     1         (gamma*r12DOTr13PT*r13*smalla)/r13s3) + 
     1      (2*R132*((r12DOTr13 + r12DOTr13PT)*r13 + 
     1           r12*(-1 + (gamma*r12DOTr13PT*r13)/r12s))*
     1         V2112)/R123 + 
     1      (r13*(2*r13*(-1 - DOT1 - r12DOTr13*r12DOTr13PT - 
     1              (gamma*r12DOTr13*r12DOTr13PT*r13)/r13s) + 
     1           r12*(4*r12DOTr13 + 4*r12DOTr13PT - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r12s + 
     1              (4*gamma*r12DOTr13PT*r13)/r13s - 
     1              (GAMMA2*DOT2*R132)/
     1               (r12s*r13s)))*V211*
     1         V311)/R122 + 
     1      (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1         (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1         (GAMMA2*DOT2*R132)/r13s4 - 
     1         (3*gamma*DOT2*R132)/r13s3 - 
     1         (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1         (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1         (gamma*DOT2*r13*smalla)/r13s3)*
     1       V3112))/(2.*R134)
                        HESS(K1+1,K3+2)=HESS(K1+1,K3+2)+ (ERL*((r13*V211*
     1         (2*r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V212 + 
     1           r12*(2*r12*r12DOTr13 + 2*r12*r12DOTr13PT - 
     1              2*DOT1*r13 - 2*r12DOTr13*r12DOTr13PT*r13 - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r13s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*R132)/
     1               r13s - 
     1              (GAMMA2*r12*DOT2*R132)/
     1               (r12s*r13s))*V312))/
     1       R123 + V311*
     1       ((2*r13*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V212)/R122 + 
     1         (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1            (GAMMA2*DOT2*R132)/r13s4 - 
     1            (3*gamma*DOT2*R132)/r13s3 - 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1            (gamma*DOT2*r13*smalla)/r13s3)*
     1          V312)))/(2.*R134)
                        HESS(K1+1,K3+3)=HESS(K1+1,K3+3)+ (ERL*((r13*V211*
     1         (2*r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V213 + 
     1           r12*(2*r12*r12DOTr13 + 2*r12*r12DOTr13PT - 
     1              2*DOT1*r13 - 2*r12DOTr13*r12DOTr13PT*r13 - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r13s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*R132)/
     1               r13s - 
     1              (GAMMA2*r12*DOT2*R132)/
     1               (r12s*r13s))*V313))/
     1       R123 + V311*
     1       ((2*r13*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V213)/R122 + 
     1         (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1            (GAMMA2*DOT2*R132)/r13s4 - 
     1            (3*gamma*DOT2*R132)/r13s3 - 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1            (gamma*DOT2*r13*smalla)/r13s3)*
     1          V313)))/(2.*R134)
                        HESS(K1+2,K1+2)=HESS(K1+2,K1+2)+ (ERL*((2*(r13*(-r12 + r12DOTr13*r13)*V212 + 
     1            r12*(r12*r12DOTr13 - r13)*V312)**2)/(R124*R134)
     1       + (4*gamma*r12DOTr13PT*
     1         (r13*(r12 - r12DOTr13*r13)*V212 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V312)*
     1         (-(V212/(r12*r12s)) - 
     1           V312/(r13*r13s)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (V212/(r12*r12s) + 
     1          V312/(r13*r13s))**2 + 
     1      gamma*DOT2*
     1       (-(1/(r12*r12s)) - 1/(r13*r13s) + 
     1         ((3*r12 - smalla)*V2122)/
     1          (R123*r12s3) + 
     1         ((3*r13 - smalla)*V3122)/(R133*r13s3))
     1       - (2*r12DOTr13PT*(R133*(2*r12 - 3*r12DOTr13*r13)*
     1            V2122 + 
     1           2*r12*r13*(R122 - r12*r12DOTr13*r13 + R132)*V212*
     1            V312 + 
     1           R122*(R132*
     1               (R122*r12DOTr13 - 2*r12*r13 + r12DOTr13*R132) + 
     1              r12*(-3*r12*r12DOTr13 + 2*r13)*V3122)))/
     1       (R124*R134)))/2.
                        HESS(K1+2,K1+3)=HESS(K1+2,K1+3)+ (ERL*((2*(r13*(-r12 + r12DOTr13*r13)*V212 + 
     1           r12*(r12*r12DOTr13 - r13)*V312)*
     1         (r13*(-r12 + r12DOTr13*r13)*V213 + 
     1           r12*(r12*r12DOTr13 - r13)*V313))/(R124*R134) + 
     1      (2*gamma*r12DOTr13PT*
     1         (-(V212/(r12*r12s)) - 
     1           V312/(r13*r13s))*
     1         (r13*(r12 - r12DOTr13*r13)*V213 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V313))/(R122*R132)
     1       + (2*gamma*r12DOTr13PT*
     1         (r13*(r12 - r12DOTr13*r13)*V212 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V312)*
     1         (-(V213/(r12*r12s)) - 
     1           V313/(r13*r13s)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (-(V212/(r12*r12s)) - 
     1         V312/(r13*r13s))*
     1       (-(V213/(r12*r12s)) - 
     1         V313/(r13*r13s)) + 
     1      gamma*DOT2*
     1       (((3*r12 - smalla)*V212*V213)/
     1          (R123*r12s3) + 
     1         ((3*r13 - smalla)*V312*V313)/
     1          (R133*r13s3)) + 
     1      (2*r12DOTr13PT*(r12*V312*
     1            (-(r13*(R122 - r12*r12DOTr13*r13 + R132)*
     1                 V213) + 
     1              R122*(3*r12*r12DOTr13 - 2*r13)*V313) - 
     1           r13*V212*
     1            (R132*(2*r12 - 3*r12DOTr13*r13)*V213 + 
     1              r12*(R122 - r12*r12DOTr13*r13 + R132)*V313)))/
     1       (R124*R134)))/2.
                        HESS(K1+2,K2+1)=HESS(K1+2,K2+1)+ (ERL*((2*r12*V311*
     1         (r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V212 + 
     1           r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V312))/R133 + 
     1      V211*((-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12*r12DOTr13)/r13 + (2*r12*r12DOTr13PT)/r13 - 
     1            (GAMMA2*R122*DOT2)/r12s4 - 
     1            (3*gamma*R122*DOT2)/r12s3 - 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*R122*r12DOTr13PT)/(r13*r12s) + 
     1            (gamma*r12*DOT2*smalla)/r12s3)*
     1          V212 + (r12*
     1            (-2*r12*DOT1 - 2*r12*r12DOTr13*r12DOTr13PT + 
     1              2*r12DOTr13*r13 + 2*r12DOTr13PT*r13 - 
     1              (2*gamma*R122*r12DOTr13*r12DOTr13PT)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r13s - 
     1              (GAMMA2*R122*DOT2*r13)/
     1               (r12s*r13s))*V312)/
     1          R132)))/(2.*R124)
                        HESS(K1+2,K2+2)=HESS(K1+2,K2+2)+ (ERL*((-2*DOT1 + 
     1         r12DOTr13*(-6*r12DOTr13PT + (2*r12)/r13 - 
     1            (4*gamma*r12*r12DOTr13PT)/r12s) + 
     1         r12*r12DOTr13PT*(2/r13 - 
     1            (GAMMA2*r12*r12DOTr13PT)/r12s4 - 
     1            (3*gamma*r12*r12DOTr13PT)/r12s3 + 
     1            (2*gamma*r12)/(r13*r12s) + 
     1            (gamma*r12DOTr13PT*smalla)/r12s3))*
     1       V2122 + (r12*
     1         (4*(r12DOTr13 + r12DOTr13PT)*r13 + 
     1           r12*(-2 - 2*DOT1 - 2*r12DOTr13*r12DOTr13PT + 
     1              (4*gamma*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s) + 
     1           (gamma*R122*r12DOTr13PT*
     1              (-2*r12DOTr13 - (gamma*r12DOTr13PT*r13)/r13s)
     1              )/r12s)*V212*V312)/R132 + 
     1      R122*(r12DOTr13PT*
     1          (2*r12DOTr13 - (2*r12)/r13 + 
     1            (gamma*R122*r12DOTr13PT)/r12s3 - 
     1            (gamma*r12*r12DOTr13PT*smalla)/r12s3) + 
     1         (2*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V3122)/R133)))/(2.*R124)
                        HESS(K1+2,K2+3)=HESS(K1+2,K2+3)+ (ERL*(V212*
     1       ((-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12*r12DOTr13)/r13 + (2*r12*r12DOTr13PT)/r13 - 
     1            (GAMMA2*R122*DOT2)/r12s4 - 
     1            (3*gamma*R122*DOT2)/r12s3 - 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*R122*r12DOTr13PT)/(r13*r12s) + 
     1            (gamma*r12*DOT2*smalla)/r12s3)*
     1          V213 + (2*r12*
     1            (-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V313)/R132) + 
     1      (r12*V312*(r13*
     1            (-2*r12*DOT1 - 2*r12*r12DOTr13*r12DOTr13PT + 
     1              2*r12DOTr13*r13 + 2*r12DOTr13PT*r13 - 
     1              (2*gamma*R122*r12DOTr13*r12DOTr13PT)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r13s - 
     1              (GAMMA2*R122*DOT2*r13)/
     1               (r12s*r13s))*V213 + 
     1           2*r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V313))/R133))/(2.*R124)
                        HESS(K1+2,K3+1)=HESS(K1+2,K3+1)+ (ERL*((2*r13*V211*
     1         (r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V212 + 
     1           r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V312))/R123 + 
     1      V311*((r13*(2*r12*r12DOTr13 + 2*r12*r12DOTr13PT - 
     1              2*DOT1*r13 - 2*r12DOTr13*r12DOTr13PT*r13 - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r13s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*R132)/
     1               r13s - 
     1              (GAMMA2*r12*DOT2*R132)/
     1               (r12s*r13s))*V212)/
     1          R122 + (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1            (GAMMA2*DOT2*R132)/r13s4 - 
     1            (3*gamma*DOT2*R132)/r13s3 - 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1            (gamma*DOT2*r13*smalla)/r13s3)*
     1          V312)))/(2.*R134)
                        HESS(K1+2,K3+2)=HESS(K1+2,K3+2)+ (ERL*(r12DOTr13PT*R132*
     1       (2*r12DOTr13 - (2*r13)/r12 + 
     1         (gamma*r12DOTr13PT*R132)/r13s3 - 
     1         (gamma*r12DOTr13PT*r13*smalla)/r13s3) + 
     1      (2*R132*((r12DOTr13 + r12DOTr13PT)*r13 + 
     1           r12*(-1 + (gamma*r12DOTr13PT*r13)/r12s))*
     1         V2122)/R123 + 
     1      (r13*(2*r13*(-1 - DOT1 - r12DOTr13*r12DOTr13PT - 
     1              (gamma*r12DOTr13*r12DOTr13PT*r13)/r13s) + 
     1           r12*(4*r12DOTr13 + 4*r12DOTr13PT - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r12s + 
     1              (4*gamma*r12DOTr13PT*r13)/r13s - 
     1              (GAMMA2*DOT2*R132)/
     1               (r12s*r13s)))*V212*
     1         V312)/R122 + 
     1      (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1         (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1         (GAMMA2*DOT2*R132)/r13s4 - 
     1         (3*gamma*DOT2*R132)/r13s3 - 
     1         (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1         (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1         (gamma*DOT2*r13*smalla)/r13s3)*
     1       V3122))/(2.*R134)
                        HESS(K1+2,K3+3)=HESS(K1+2,K3+3)+ (ERL*((r13*V212*
     1         (2*r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V213 + 
     1           r12*(2*r12*r12DOTr13 + 2*r12*r12DOTr13PT - 
     1              2*DOT1*r13 - 2*r12DOTr13*r12DOTr13PT*r13 - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r13s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*R132)/
     1               r13s - 
     1              (GAMMA2*r12*DOT2*R132)/
     1               (r12s*r13s))*V313))/
     1       R123 + V312*
     1       ((2*r13*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V213)/R122 + 
     1         (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1            (GAMMA2*DOT2*R132)/r13s4 - 
     1            (3*gamma*DOT2*R132)/r13s3 - 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1            (gamma*DOT2*r13*smalla)/r13s3)*
     1          V313)))/(2.*R134)
                        HESS(K1+3,K1+3)=HESS(K1+3,K1+3)+ (ERL*((2*(r13*(-r12 + r12DOTr13*r13)*V213 + 
     1            r12*(r12*r12DOTr13 - r13)*V313)**2)/(R124*R134)
     1       + (4*gamma*r12DOTr13PT*
     1         (r13*(r12 - r12DOTr13*r13)*V213 + 
     1           r12*(-(r12*r12DOTr13) + r13)*V313)*
     1         (-(V213/(r12*r12s)) - 
     1           V313/(r13*r13s)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (V213/(r12*r12s) + 
     1          V313/(r13*r13s))**2 + 
     1      gamma*DOT2*
     1       (-(1/(r12*r12s)) - 1/(r13*r13s) + 
     1         ((3*r12 - smalla)*V2132)/
     1          (R123*r12s3) + 
     1         ((3*r13 - smalla)*V3132)/(R133*r13s3))
     1       - (2*r12DOTr13PT*(R133*(2*r12 - 3*r12DOTr13*r13)*
     1            V2132 + 
     1           2*r12*r13*(R122 - r12*r12DOTr13*r13 + R132)*V213*
     1            V313 + 
     1           R122*(R132*
     1               (R122*r12DOTr13 - 2*r12*r13 + r12DOTr13*R132) + 
     1              r12*(-3*r12*r12DOTr13 + 2*r13)*V3132)))/
     1       (R124*R134)))/2.
                        HESS(K1+3,K2+1)=HESS(K1+3,K2+1)+ (ERL*((2*r12*V311*
     1         (r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V213 + 
     1           r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V313))/R133 + 
     1      V211*((-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12*r12DOTr13)/r13 + (2*r12*r12DOTr13PT)/r13 - 
     1            (GAMMA2*R122*DOT2)/r12s4 - 
     1            (3*gamma*R122*DOT2)/r12s3 - 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*R122*r12DOTr13PT)/(r13*r12s) + 
     1            (gamma*r12*DOT2*smalla)/r12s3)*
     1          V213 + (r12*
     1            (-2*r12*DOT1 - 2*r12*r12DOTr13*r12DOTr13PT + 
     1              2*r12DOTr13*r13 + 2*r12DOTr13PT*r13 - 
     1              (2*gamma*R122*r12DOTr13*r12DOTr13PT)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r13s - 
     1              (GAMMA2*R122*DOT2*r13)/
     1               (r12s*r13s))*V313)/
     1          R132)))/(2.*R124)
                        HESS(K1+3,K2+2)=HESS(K1+3,K2+2)+ (ERL*((2*r12*V312*
     1         (r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V213 + 
     1           r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V313))/R133 + 
     1      V212*((-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12*r12DOTr13)/r13 + (2*r12*r12DOTr13PT)/r13 - 
     1            (GAMMA2*R122*DOT2)/r12s4 - 
     1            (3*gamma*R122*DOT2)/r12s3 - 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*R122*r12DOTr13PT)/(r13*r12s) + 
     1            (gamma*r12*DOT2*smalla)/r12s3)*
     1          V213 + (r12*
     1            (-2*r12*DOT1 - 2*r12*r12DOTr13*r12DOTr13PT + 
     1              2*r12DOTr13*r13 + 2*r12DOTr13PT*r13 - 
     1              (2*gamma*R122*r12DOTr13*r12DOTr13PT)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r13s - 
     1              (GAMMA2*R122*DOT2*r13)/
     1               (r12s*r13s))*V313)/
     1          R132)))/(2.*R124)
                        HESS(K1+3,K2+3)=HESS(K1+3,K2+3)+ (ERL*((-2*DOT1 + 
     1         r12DOTr13*(-6*r12DOTr13PT + (2*r12)/r13 - 
     1            (4*gamma*r12*r12DOTr13PT)/r12s) + 
     1         r12*r12DOTr13PT*(2/r13 - 
     1            (GAMMA2*r12*r12DOTr13PT)/r12s4 - 
     1            (3*gamma*r12*r12DOTr13PT)/r12s3 + 
     1            (2*gamma*r12)/(r13*r12s) + 
     1            (gamma*r12DOTr13PT*smalla)/r12s3))*
     1       V2132 + (r12*
     1         (4*(r12DOTr13 + r12DOTr13PT)*r13 + 
     1           r12*(-2 - 2*DOT1 - 2*r12DOTr13*r12DOTr13PT + 
     1              (4*gamma*r12DOTr13PT*r13)/r12s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s) + 
     1           (gamma*R122*r12DOTr13PT*
     1              (-2*r12DOTr13 - (gamma*r12DOTr13PT*r13)/r13s)
     1              )/r12s)*V213*V313)/R132 + 
     1      R122*(r12DOTr13PT*
     1          (2*r12DOTr13 - (2*r12)/r13 + 
     1            (gamma*R122*r12DOTr13PT)/r12s3 - 
     1            (gamma*r12*r12DOTr13PT*smalla)/r12s3) + 
     1         (2*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V3132)/R133)))/(2.*R124)
                        HESS(K1+3,K3+1)=HESS(K1+3,K3+1)+ (ERL*((2*r13*V211*
     1         (r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V213 + 
     1           r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V313))/R123 + 
     1      V311*((r13*(2*r12*r12DOTr13 + 2*r12*r12DOTr13PT - 
     1              2*DOT1*r13 - 2*r12DOTr13*r12DOTr13PT*r13 - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r13s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*R132)/
     1               r13s - 
     1              (GAMMA2*r12*DOT2*R132)/
     1               (r12s*r13s))*V213)/
     1          R122 + (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1            (GAMMA2*DOT2*R132)/r13s4 - 
     1            (3*gamma*DOT2*R132)/r13s3 - 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1            (gamma*DOT2*r13*smalla)/r13s3)*
     1          V313)))/(2.*R134)
                        HESS(K1+3,K3+2)=HESS(K1+3,K3+2)+ (ERL*((2*r13*V212*
     1         (r13*(-r12 + r12DOTr13*r13 + r12DOTr13PT*r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r12s)*
     1            V213 + 
     1           r12*(r12*r12DOTr13 + r12*r12DOTr13PT - r13 + 
     1              (gamma*r12*r12DOTr13PT*r13)/r13s)*
     1            V313))/R123 + 
     1      V312*((r13*(2*r12*r12DOTr13 + 2*r12*r12DOTr13PT - 
     1              2*DOT1*r13 - 2*r12DOTr13*r12DOTr13PT*r13 - 
     1              (2*gamma*r12*r12DOTr13*r12DOTr13PT*r13)/
     1               r12s + 
     1              (2*gamma*r12*r12DOTr13PT*r13)/r13s - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*R132)/
     1               r13s - 
     1              (GAMMA2*r12*DOT2*R132)/
     1               (r12s*r13s))*V213)/
     1          R122 + (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1            (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1            (GAMMA2*DOT2*R132)/r13s4 - 
     1            (3*gamma*DOT2*R132)/r13s3 - 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1            (gamma*DOT2*r13*smalla)/r13s3)*
     1          V313)))/(2.*R134)
                        HESS(K1+3,K3+3)=HESS(K1+3,K3+3)+ (ERL*(r12DOTr13PT*R132*
     1       (2*r12DOTr13 - (2*r13)/r12 + 
     1         (gamma*r12DOTr13PT*R132)/r13s3 - 
     1         (gamma*r12DOTr13PT*r13*smalla)/r13s3) + 
     1      (2*R132*((r12DOTr13 + r12DOTr13PT)*r13 + 
     1           r12*(-1 + (gamma*r12DOTr13PT*r13)/r12s))*
     1         V2132)/R123 + 
     1      (r13*(2*r13*(-1 - DOT1 - r12DOTr13*r12DOTr13PT - 
     1              (gamma*r12DOTr13*r12DOTr13PT*r13)/r13s) + 
     1           r12*(4*r12DOTr13 + 4*r12DOTr13PT - 
     1              (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r12s + 
     1              (4*gamma*r12DOTr13PT*r13)/r13s - 
     1              (GAMMA2*DOT2*R132)/
     1               (r12s*r13s)))*V213*
     1         V313)/R122 + 
     1      (-2*DOT1 - 6*r12DOTr13*r12DOTr13PT + 
     1         (2*r12DOTr13*r13)/r12 + (2*r12DOTr13PT*r13)/r12 - 
     1         (GAMMA2*DOT2*R132)/r13s4 - 
     1         (3*gamma*DOT2*R132)/r13s3 - 
     1         (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1         (2*gamma*r12DOTr13PT*R132)/(r12*r13s) + 
     1         (gamma*DOT2*r13*smalla)/r13s3)*
     1       V3132))/(2.*R134)
                        HESS(K2+1,K2+1)=HESS(K2+1,K2+1)+ (ERL*((2*DOT1 + 
     1         r12DOTr13*(6*r12DOTr13PT + 
     1            (4*gamma*r12*r12DOTr13PT)/r12s) + 
     1         (gamma*r12*DOT2*
     1            (gamma*r12 + 3*R122 - 4*r12*smalla + smalla**2))/
     1          r12s4)*V2112 - 
     1      (4*r12*(R122*(r12DOTr13 + r12DOTr13PT) + 
     1           (r12DOTr13 + r12DOTr13PT)*smalla**2 + 
     1           r12*(gamma*r12DOTr13PT - 2*r12DOTr13*smalla - 
     1              2*r12DOTr13PT*smalla))*V211*V311)/
     1       (r13*r12s) + 
     1      R122*(r12DOTr13PT*
     1          (-2*r12DOTr13 - (gamma*r12*r12DOTr13PT)/r12s) + 
     1         (2*V3112)/R132)))/(2.*R124)
                        HESS(K2+1,K2+2)=HESS(K2+1,K2+2)+ (ERL*((2*r12*V311*
     1         (r13*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12*r12DOTr13PT)/r12s)*V212 + 
     1           r12*V312))/R132 + 
     1      V211*((2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1            (GAMMA2*R122*DOT2)/r12s4 + 
     1            (2*gamma*R122*DOT2)/r12s3 + 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (gamma*r12*DOT2)/r12s)*V212
     1          + (2*r12*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12*r12DOTr13PT)/r12s)*V312)/
     1          r13)))/(2.*R124)
                        HESS(K2+1,K2+3)=HESS(K2+1,K2+3)+ (ERL*((2*r12*V311*
     1         (r13*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12*r12DOTr13PT)/r12s)*V213 + 
     1           r12*V313))/R132 + 
     1      V211*((2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1            (GAMMA2*R122*DOT2)/r12s4 + 
     1            (2*gamma*R122*DOT2)/r12s3 + 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (gamma*r12*DOT2)/r12s)*V213
     1          + (2*r12*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12*r12DOTr13PT)/r12s)*V313)/
     1          r13)))/(2.*R124)
                        HESS(K2+1,K3+1)=HESS(K2+1,K3+1)+ (ERL*(2*R132*
     1       (-r12DOTr13 + r12DOTr13PT*(-1 - (gamma*r12)/r12s))*
     1       V2112 + r12*r13*
     1       (2 + 2*DOT1 + 
     1         2*r12DOTr13*r12DOTr13PT*
     1          (1 + (gamma*r12)/r12s + 
     1            (gamma*r13)/r13s) + 
     1         (GAMMA2*r12*DOT2*r13)/
     1          (r12s*r13s))*V211*
     1       V311 + 2*R122*
     1       (r12DOTr13PT*R132 + 
     1         (-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V3112))
     1    )/(2.*R123*R133)
                        HESS(K2+1,K3+2)=HESS(K2+1,K3+2)+ (ERL*(2*r12*V311*
     1       (r13*V212 + 
     1         r12*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V312) + 
     1      r13*V211*((-2*r12DOTr13*r13 - 2*r12DOTr13PT*r13 - 
     1            (2*gamma*r12*r12DOTr13PT*r13)/r12s)*
     1          V212 + r12*
     1          (2*DOT1 + 2*r12DOTr13*r12DOTr13PT + 
     1            (2*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (GAMMA2*r12*DOT2*r13)/
     1             (r12s*r13s))*V312)))/
     1  (2.*R123*R133)
                        HESS(K2+1,K3+3)=HESS(K2+1,K3+3)+ (ERL*(2*r12*V311*
     1       (r13*V213 + 
     1         r12*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V313) + 
     1      r13*V211*((-2*r12DOTr13*r13 - 2*r12DOTr13PT*r13 - 
     1            (2*gamma*r12*r12DOTr13PT*r13)/r12s)*
     1          V213 + r12*
     1          (2*DOT1 + 2*r12DOTr13*r12DOTr13PT + 
     1            (2*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (GAMMA2*r12*DOT2*r13)/
     1             (r12s*r13s))*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+2,K2+2)=HESS(K2+2,K2+2)+ (ERL*((2*DOT1 + 
     1         r12DOTr13*(6*r12DOTr13PT + 
     1            (4*gamma*r12*r12DOTr13PT)/r12s) + 
     1         (gamma*r12*DOT2*
     1            (gamma*r12 + 3*R122 - 4*r12*smalla + smalla**2))/
     1          r12s4)*V2122 - 
     1      (4*r12*(R122*(r12DOTr13 + r12DOTr13PT) + 
     1           (r12DOTr13 + r12DOTr13PT)*smalla**2 + 
     1           r12*(gamma*r12DOTr13PT - 2*r12DOTr13*smalla - 
     1              2*r12DOTr13PT*smalla))*V212*V312)/
     1       (r13*r12s) + 
     1      R122*(r12DOTr13PT*
     1          (-2*r12DOTr13 - (gamma*r12*r12DOTr13PT)/r12s) + 
     1         (2*V3122)/R132)))/(2.*R124)
                        HESS(K2+2,K2+3)=HESS(K2+2,K2+3)+ (ERL*((2*r12*V312*
     1         (r13*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12*r12DOTr13PT)/r12s)*V213 + 
     1           r12*V313))/R132 + 
     1      V212*((2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1            (GAMMA2*R122*DOT2)/r12s4 + 
     1            (2*gamma*R122*DOT2)/r12s3 + 
     1            (4*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (gamma*r12*DOT2)/r12s)*V213
     1          + (2*r12*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12*r12DOTr13PT)/r12s)*V313)/
     1          r13)))/(2.*R124)
                        HESS(K2+2,K3+1)=HESS(K2+2,K3+1)+ (ERL*(2*r13*V211*
     1       (r13*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12*r12DOTr13PT)/r12s)*V212 + 
     1         r12*V312) + 
     1      r12*V311*(r13*
     1          (2*DOT1 + 2*r12DOTr13*r12DOTr13PT + 
     1            (2*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (GAMMA2*r12*DOT2*r13)/
     1             (r12s*r13s))*V212 + 
     1         2*r12*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V312)))/
     1  (2.*R123*R133)
                        HESS(K2+2,K3+2)=HESS(K2+2,K3+2)+ (ERL*(2*R132*
     1       (-r12DOTr13 + r12DOTr13PT*(-1 - (gamma*r12)/r12s))*
     1       V2122 + r12*r13*
     1       (2 + 2*DOT1 + 
     1         2*r12DOTr13*r12DOTr13PT*
     1          (1 + (gamma*r12)/r12s + 
     1            (gamma*r13)/r13s) + 
     1         (GAMMA2*r12*DOT2*r13)/
     1          (r12s*r13s))*V212*
     1       V312 + 2*R122*
     1       (r12DOTr13PT*R132 + 
     1         (-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V3122))
     1    )/(2.*R123*R133)
                        HESS(K2+2,K3+3)=HESS(K2+2,K3+3)+ (ERL*(2*r12*V312*
     1       (r13*V213 + 
     1         r12*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V313) + 
     1      r13*V212*((-2*r12DOTr13*r13 - 2*r12DOTr13PT*r13 - 
     1            (2*gamma*r12*r12DOTr13PT*r13)/r12s)*
     1          V213 + r12*
     1          (2*DOT1 + 2*r12DOTr13*r12DOTr13PT + 
     1            (2*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (GAMMA2*r12*DOT2*r13)/
     1             (r12s*r13s))*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+3,K2+3)=HESS(K2+3,K2+3)+ (ERL*((2*DOT1 + 
     1         r12DOTr13*(6*r12DOTr13PT + 
     1            (4*gamma*r12*r12DOTr13PT)/r12s) + 
     1         (gamma*r12*DOT2*
     1            (gamma*r12 + 3*R122 - 4*r12*smalla + smalla**2))/
     1          r12s4)*V2132 - 
     1      (4*r12*(R122*(r12DOTr13 + r12DOTr13PT) + 
     1           (r12DOTr13 + r12DOTr13PT)*smalla**2 + 
     1           r12*(gamma*r12DOTr13PT - 2*r12DOTr13*smalla - 
     1              2*r12DOTr13PT*smalla))*V213*V313)/
     1       (r13*r12s) + 
     1      R122*(r12DOTr13PT*
     1          (-2*r12DOTr13 - (gamma*r12*r12DOTr13PT)/r12s) + 
     1         (2*V3132)/R132)))/(2.*R124)
                        HESS(K2+3,K3+1)=HESS(K2+3,K3+1)+ (ERL*(2*r13*V211*
     1       (r13*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12*r12DOTr13PT)/r12s)*V213 + 
     1         r12*V313) + 
     1      r12*V311*(r13*
     1          (2*DOT1 + 2*r12DOTr13*r12DOTr13PT + 
     1            (2*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (GAMMA2*r12*DOT2*r13)/
     1             (r12s*r13s))*V213 + 
     1         2*r12*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+3,K3+2)=HESS(K2+3,K3+2)+ (ERL*(2*r13*V212*
     1       (r13*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12*r12DOTr13PT)/r12s)*V213 + 
     1         r12*V313) + 
     1      r12*V312*(r13*
     1          (2*DOT1 + 2*r12DOTr13*r12DOTr13PT + 
     1            (2*gamma*r12*r12DOTr13*r12DOTr13PT)/r12s + 
     1            (2*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (GAMMA2*r12*DOT2*r13)/
     1             (r12s*r13s))*V213 + 
     1         2*r12*(-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+3,K3+3)=HESS(K2+3,K3+3)+ (ERL*(2*R132*
     1       (-r12DOTr13 + r12DOTr13PT*(-1 - (gamma*r12)/r12s))*
     1       V2132 + r12*r13*
     1       (2 + 2*DOT1 + 
     1         2*r12DOTr13*r12DOTr13PT*
     1          (1 + (gamma*r12)/r12s + 
     1            (gamma*r13)/r13s) + 
     1         (GAMMA2*r12*DOT2*r13)/
     1          (r12s*r13s))*V213*
     1       V313 + 2*R122*
     1       (r12DOTr13PT*R132 + 
     1         (-r12DOTr13 - r12DOTr13PT - 
     1            (gamma*r12DOTr13PT*r13)/r13s)*V3132))
     1    )/(2.*R123*R133)
                        HESS(K3+1,K3+1)=HESS(K3+1,K3+1)+ (ERL*(r12DOTr13PT*R132*
     1       (-2*r12DOTr13 - (gamma*r12DOTr13PT*r13)/r13s) + 
     1      (2*R132*V2112)/R122 - 
     1      (4*r13*(gamma*r12DOTr13PT*r13 + 
     1           (r12DOTr13 + r12DOTr13PT)*r13s)*V211*
     1         V311)/(r12*r13s) + 
     1      (2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1         (GAMMA2*DOT2*R132)/r13s4 + 
     1         (2*gamma*DOT2*R132)/r13s3 + 
     1         (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1         (gamma*DOT2*r13)/r13s)*V3112))
     1   /(2.*R134)
                        HESS(K3+1,K3+2)=HESS(K3+1,K3+2)+ (ERL*((2*r13*V211*
     1         (r13*V212 + 
     1           r12*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12DOTr13PT*r13)/r13s)*V312))/
     1       R122 + V311*
     1       ((2*r13*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12DOTr13PT*r13)/r13s)*V212)/
     1          r12 + (2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1            (GAMMA2*DOT2*R132)/r13s4 + 
     1            (2*gamma*DOT2*R132)/r13s3 + 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (gamma*DOT2*r13)/r13s)*V312))
     1    )/(2.*R134)
                        HESS(K3+1,K3+3)=HESS(K3+1,K3+3)+ (ERL*((2*r13*V211*
     1         (r13*V213 + 
     1           r12*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12DOTr13PT*r13)/r13s)*V313))/
     1       R122 + V311*
     1       ((2*r13*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12DOTr13PT*r13)/r13s)*V213)/
     1          r12 + (2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1            (GAMMA2*DOT2*R132)/r13s4 + 
     1            (2*gamma*DOT2*R132)/r13s3 + 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (gamma*DOT2*r13)/r13s)*V313))
     1    )/(2.*R134)
                        HESS(K3+2,K3+2)=HESS(K3+2,K3+2)+ (ERL*(r12DOTr13PT*R132*
     1       (-2*r12DOTr13 - (gamma*r12DOTr13PT*r13)/r13s) + 
     1      (2*R132*V2122)/R122 - 
     1      (4*r13*(gamma*r12DOTr13PT*r13 + 
     1           (r12DOTr13 + r12DOTr13PT)*r13s)*V212*
     1         V312)/(r12*r13s) + 
     1      (2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1         (GAMMA2*DOT2*R132)/r13s4 + 
     1         (2*gamma*DOT2*R132)/r13s3 + 
     1         (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1         (gamma*DOT2*r13)/r13s)*V3122))
     1   /(2.*R134)
                        HESS(K3+2,K3+3)=HESS(K3+2,K3+3)+ (ERL*((2*r13*V212*
     1         (r13*V213 + 
     1           r12*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12DOTr13PT*r13)/r13s)*V313))/
     1       R122 + V312*
     1       ((2*r13*(-r12DOTr13 - r12DOTr13PT - 
     1              (gamma*r12DOTr13PT*r13)/r13s)*V213)/
     1          r12 + (2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1            (GAMMA2*DOT2*R132)/r13s4 + 
     1            (2*gamma*DOT2*R132)/r13s3 + 
     1            (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1            (gamma*DOT2*r13)/r13s)*V313))
     1    )/(2.*R134)
                        HESS(K3+3,K3+3)=HESS(K3+3,K3+3)+ (ERL*(r12DOTr13PT*R132*
     1       (-2*r12DOTr13 - (gamma*r12DOTr13PT*r13)/r13s) + 
     1      (2*R132*V2132)/R122 - 
     1      (4*r13*(gamma*r12DOTr13PT*r13 + 
     1           (r12DOTr13 + r12DOTr13PT)*r13s)*V213*
     1         V313)/(r12*r13s) + 
     1      (2*DOT1 + 6*r12DOTr13*r12DOTr13PT + 
     1         (GAMMA2*DOT2*R132)/r13s4 + 
     1         (2*gamma*DOT2*R132)/r13s3 + 
     1         (4*gamma*r12DOTr13*r12DOTr13PT*r13)/r13s + 
     1         (gamma*DOT2*r13)/r13s)*V3132))
     1   /(2.*R134)

                        HESS(K1+2,K1+1)=HESS(K1+1,K1+2) 
                        HESS(K1+3,K1+1)=HESS(K1+1,K1+3) 
                        HESS(K2+1,K1+1)=HESS(K1+1,K2+1)
                        HESS(K2+2,K1+1)=HESS(K1+1,K2+2)
                        HESS(K2+3,K1+1)=HESS(K1+1,K2+3)
                        HESS(K3+1,K1+1)=HESS(K1+1,K3+1)
                        HESS(K3+2,K1+1)=HESS(K1+1,K3+2)
                        HESS(K3+3,K1+1)=HESS(K1+1,K3+3)

                        HESS(K1+3,K1+2)=HESS(K1+2,K1+3)
                        HESS(K2+1,K1+2)=HESS(K1+2,K2+1)
                        HESS(K2+2,K1+2)=HESS(K1+2,K2+2)
                        HESS(K2+3,K1+2)=HESS(K1+2,K2+3)
                        HESS(K3+1,K1+2)=HESS(K1+2,K3+1)
                        HESS(K3+2,K1+2)=HESS(K1+2,K3+2)
                        HESS(K3+3,K1+2)=HESS(K1+2,K3+3)

                        HESS(K2+1,K1+3)=HESS(K1+3,K2+1)
                        HESS(K2+2,K1+3)=HESS(K1+3,K2+2)
                        HESS(K2+3,K1+3)=HESS(K1+3,K2+3)
                        HESS(K3+1,K1+3)=HESS(K1+3,K3+1)
                        HESS(K3+2,K1+3)=HESS(K1+3,K3+2)
                        HESS(K3+3,K1+3)=HESS(K1+3,K3+3)

                        HESS(K2+2,K2+1)=HESS(K2+1,K2+2)
                        HESS(K2+3,K2+1)=HESS(K2+1,K2+3)
                        HESS(K3+1,K2+1)=HESS(K2+1,K3+1)
                        HESS(K3+2,K2+1)=HESS(K2+1,K3+2)
                        HESS(K3+3,K2+1)=HESS(K2+1,K3+3)

                        HESS(K2+3,K2+2)=HESS(K2+2,K2+3)
                        HESS(K3+1,K2+2)=HESS(K2+2,K3+1)
                        HESS(K3+2,K2+2)=HESS(K2+2,K3+2)
                        HESS(K3+3,K2+2)=HESS(K2+2,K3+3)

                        HESS(K3+1,K2+3)=HESS(K2+3,K3+1)
                        HESS(K3+2,K2+3)=HESS(K2+3,K3+2)
                        HESS(K3+3,K2+3)=HESS(K2+3,K3+3)

                        HESS(K3+2,K3+1)=HESS(K3+1,K3+2)
                        HESS(K3+3,K3+1)=HESS(K3+1,K3+3)

                        HESS(K3+3,K3+2)=HESS(K3+2,K3+3)

                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      ENERGY=ENERGY*LAMBDA/2.0D0

      RETURN
      END
