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
C  Energy and derivatives of a general Born-Mayer potential with
C  Tosi-Fumi parameters specified after the keyword TOSI in the
C  odata file. Energy is in hartree, length in Bohr.
C
      SUBROUTINE TOSIFUMI(N, X, V, POTEL, APP, AMM, APM, RHO, ZSYM, GTEST, STEST)
      use porfuncs
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, I, J, J1, J2, K, L, J3
      DOUBLE PRECISION X(3*N), V(3*N), POTEL, 
     1                 APP, AMM, APM, RHO, DUM
      DOUBLE PRECISION AC(N,N),Q(N),RD(N,N)
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST

      DO I=1,N
         IF (ZSYM(I).EQ.'MI') THEN
            Q(I)=-1.0D0
         ELSE IF (ZSYM(I).EQ.'PL') THEN
            Q(I)=1.0D0
         ELSE
            WRITE(*,'(A)') ' All atoms must be type PL or MI for Tosi-Fumi'
            STOP
         ENDIF
      ENDDO

       DO I=1,N
          DO J=I,N
             IF (ZSYM(I).EQ.'PL') THEN 
                IF (ZSYM(J).EQ.'PL') THEN 
                   AC(I,J)=APP
                   AC(J,I)=APP
                ELSE IF (ZSYM(J).EQ.'MI') THEN 
                   AC(I,J)=APM
                   AC(J,I)=APM
                ENDIF
             ELSE IF (ZSYM(I).EQ.'MI') THEN
                IF (ZSYM(J).EQ.'MI') THEN 
                   AC(I,J)=AMM
                   AC(J,I)=AMM
                ELSE IF (ZSYM(J).EQ.'PL') THEN 
                   AC(I,J)=APM
                   AC(J,I)=APM
                ENDIF
             ENDIF 
         ENDDO
      ENDDO
C
C  Calculate the interparticle distances.
C 
      DO I=1,N
         K=1+3*(I-1)         
         DO J=I+1,N
            L=K+3*(J-I)
            RD(I,J)= DSQRT((X(K)-X(L))*(X(K)-X(L))+
     1                     (X(K+1)-X(L+1))*(X(K+1)-X(L+1))+
     2                     (X(K+2)-X(L+2))*(X(K+2)-X(L+2))) 
            RD(J,I)=RD(I,J)
         ENDDO
      ENDDO

      POTEL=0.0D0
      DO I=1,N
         DO J=I+1,N
            POTEL=POTEL+Q(I)*Q(J)/RD(J,I)+AC(J,I)*DEXP(-RD(J,I)/RHO)
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
C
C  Gradient.
C
      DO J1=1,N
         DO J2=1,3
            DUM=0.0D0
            DO J3=1,N
               IF (J1.NE.J3) THEN
                  DUM=DUM-(Q(J3)*Q(J1)/RD(J3,J1)**3+AC(J3,J1)*DEXP(-RD(J3,J1)/RHO)/(RD(J3,J1)*RHO))
     1                       *(X(3*(J1-1)+J2)-X(3*(J3-1)+J2))
               ENDIF
            ENDDO
            V(3*(J1-1)+J2)=DUM
         ENDDO
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  Hessian.
C
      DO I=1,3*N
         DO J=1,3*N
            HESS(I,J)=0.0D0
         ENDDO
      ENDDO
      DO I=1,3*N
         DO J=I,3*N
C 
C  Determine the cartesian coordinate of X(I) and X(J).
C
            K=MOD((I-1),3)
            L=MOD((J-1),3)
C 
C  This IF statement calls SHESS if I and J are the same
C  Cartesian coordinate and calls DHESS if the are different
C  Cartesian coordinates
C
            IF(K .EQ. L)THEN
               CALL NSHESS(I,J,K,L,X,N,RHO,AC,Q,RD)
            ELSE
               CALL NDHESS(I,J,K,L,X,N,RHO,AC,Q,RD)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
C....................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE THE SAME
 
      SUBROUTINE NSHESS(I,J,K,L,CFG,NUM,RHO,AC,Q,RD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM
      DOUBLE PRECISION CFG(3*NUM), RW, RHO, R, QQ, AB, 
     1                 EP, DR, DF
      DOUBLE PRECISION AC(NUM,NUM),Q(NUM),RD(NUM,NUM)

      RW=RHO
 
c.....I1 and J1 are the particle numbers for X(I) and X(J)

      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
c.....if the particles are the same do the next two loops
c.....there are two loops so that we can skip particle I1
c.....without resorting to an IF statement

      IF(I1 .EQ. J1)THEN
      
      DO 10 M = 1, I1-1
        
      N = 3*(M-1) + K + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
      DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1     + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2     - (AB*EP/(RW*R))
 
      HESS(I,I) = HESS(I,I) + DF
      
   10 CONTINUE
 
      DO 20 M = I1+1, NUM
         
      N = 3*(M-1) + K + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
      DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1     + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2     - (AB*EP/(RW*R))
 
      HESS(I,I) = HESS(I,I) + DF
      
   20 CONTINUE
 
      ELSE
 
c.....do this if the particles are not the same
      R = RD(I1,J1)
      QQ = Q(I1)*Q(J1)
      AB = AC(I1,J1)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(J))*(CFG(I)-CFG(J))
      
      HESS(I,J) = (QQ/(R*R*R)) - (3.0D0*QQ*DR/(R*R*R*R*R))
     1     - (AB*DR*EP/(RW*R*R*R)) - (AB*DR*EP/(RW*RW*R*R))
     2     + (AB*EP/(RW*R))
 
      HESS(J,I) = HESS(I,J)
 
      ENDIF
 
      RETURN
      END
C..........................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE NOT THE SAME
 
      SUBROUTINE NDHESS(I,J,K,L,CFG,NUM,RHO,AC,Q,RD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM,N1
      DOUBLE PRECISION CFG(3*NUM),RHO,RW,R,QQ,AB,
     1                 EP,DR,DF
      DOUBLE PRECISION AC(NUM,NUM),Q(NUM),RD(NUM,NUM)
      RW=RHO
      
      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
c.....do the next to loops if the particles are the same

      IF(I1 .EQ. J1)THEN
 
      DO 10 M = 1, I1-1
 
      N = 3*(M-1) + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N+K))*(CFG(J)-CFG(N+L))
      
      DF = (3.0D0*QQ*DR/(R*R*R*R*R))
     1    + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
 
      HESS(I,J) = HESS(I,J) + DF
      
   10 CONTINUE
 
      DO 20 M = I1+1 ,NUM
         
      N = 3*(M-1) + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N+K))*(CFG(J)-CFG(N+L))
      
      DF = (3.0D0*QQ*DR/(R*R*R*R*R))
     1    + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
 
      HESS(I,J) = HESS(I,J) + DF
      
   20 CONTINUE
      HESS(J,I) = HESS(I,J)
 
      ELSE
 
c.....do this if the particles are not the same
      N = K + 3*(J1-1) + 1
      N1 = L + 3*(I1-1) + 1
 
      R = RD(I1,J1)
      QQ = Q(I1)*Q(J1)
      AB = AC(I1,J1)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(N1)-CFG(J))
      
      HESS(I,J) = -(3.0D0*QQ*DR/(R*R*R*R*R))
     1         - (AB*DR*EP/(RW*R*R*R)) - (AB*DR*EP/(RW*RW*R*R))
 
      HESS(J,I) = HESS(I,J)
 
      ENDIF
 
      RETURN
      END
C
C  Energy and analytic derivatives of possible dispersion terms.
C
      SUBROUTINE TOSIFUMIC6(N, X, V, EDISP, C6PP, C6MM, C6PM, ZSYM, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J6, J4, J5
      DOUBLE PRECISION X(3*N), R6,
     1                 V(3*N), R2(N,N), DUMMY3, DUMMY4,
     2                 R8(N,N), G(N,N), 
     3                 F(N,N), DUMMY, EDISP, C6PP, C6MM, C6PM
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST
C 
C  Store distance matrices.
C
      EDISP=0.0D0
      DO J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         DO J2=J1+1,N
            R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R6=R2(J2,J1)**3
            IF (ZSYM(J2).EQ.'PL') THEN
               IF (ZSYM(J1).EQ.'PL') THEN
                  EDISP=EDISP-C6PP*R6
               ELSE
                  EDISP=EDISP-C6PM*R6
               ENDIF
            ELSE
               IF (ZSYM(J1).EQ.'PL') THEN
                  EDISP=EDISP-C6PM*R6
               ELSE
                  EDISP=EDISP-C6MM*R6
               ENDIF
            ENDIF
            R8(J2,J1)=R2(J2,J1)**4
            R2(J1,J2)=R2(J2,J1)
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
C
C  Calculate the g tensor.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (ZSYM(J2).EQ.'PL') THEN
               IF (ZSYM(J1).EQ.'PL') THEN
                  G(J2,J1)=6.0D0*C6PP*R8(J2,J1)
               ELSE
                  G(J2,J1)=6.0D0*C6PM*R8(J2,J1)
               ENDIF
            ELSE
               IF (ZSYM(J1).EQ.'PL') THEN
                  G(J2,J1)=6.0D0*C6PM*R8(J2,J1)
               ELSE
                  G(J2,J1)=6.0D0*C6MM*R8(J2,J1)
               ENDIF
            ENDIF
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
C  First calculate the gradient analytically.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=V(J3)+DUMMY
         ENDDO
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  f tensor
C
      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            IF (ZSYM(J2).EQ.'PL') THEN
               IF (ZSYM(J1).EQ.'PL') THEN
                  F(J2,J1)=-42.0D0*C6PP*R8(J2,J1)-G(J2,J1)
               ELSE
                  F(J2,J1)=-42.0D0*C6PM*R8(J2,J1)-G(J2,J1)
               ENDIF
            ELSE
               IF (ZSYM(J1).EQ.'PL') THEN
                  F(J2,J1)=-42.0D0*C6PM*R8(J2,J1)-G(J2,J1)
               ELSE
                  F(J2,J1)=-42.0D0*C6MM*R8(J2,J1)-G(J2,J1)
               ENDIF
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
            HESS(J3,J3)=HESS(J3,J3)+DUMMY
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
                  DUMMY3=F(J5,J1)*R2(J5,J1)
                  DUMMY4=(X(J3)-X(3*J5-3+J2))*(X(3*J1-3+J4)-X(3*J5-3+J4))
                  DUMMY=DUMMY + DUMMY3*DUMMY4
               ENDDO
               HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+DUMMY
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
               HESS(3*(J4-1)+J2,J3)=HESS(3*(J4-1)+J2,J3)-F(J4,J1)*R2(J4,J1)*
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
                  HESS(J6,J3)=HESS(J6,J3)-F(J4,J1)*R2(J4,J1)
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

C**************************************************************************
C
C  Energy and analytic derivatives of possible first order induction.
C
      SUBROUTINE TOSIFUMIPOL(N, X, V, EIND, ALPHAP, ALPHAM, ZSYM, DAMP, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2
      DOUBLE PRECISION X(3*N), EIND, DAMP,
     1                 V(3*N), 
     2                 ALPHAP, ALPHAM, VEC1(3*N),
     3                 VEC2(3*N), DIF, TEMP1
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST

      CALL TFIND(N, X, EIND, ALPHAP, ALPHAM, ZSYM, DAMP)

      IF (.NOT.GTEST) RETURN
C
C  Analytic gradient.
C
      CALL TFGRAD(N, X, VEC1, ALPHAP, ALPHAM, ZSYM, DAMP)
      DO J1=1,3*N
         V(J1)=V(J1)+VEC1(J1)
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  Analytic second derivatives. Not done for damping!
C
C     PRINT*,'Analytic second derivatives:'
C     CALL TFSEC(N, X, ALPHAP, ALPHAM, ZSYM, DAMP, NPROD)
C     DO J1=1,3*N
C        DO J2=1,3*N
C           WRITE(*,'(A,2I4,F20.10)') 'J1,J2,A=',J1,J2,HESS(J2,J1)
C        ENDDO
C     ENDDO
C
C  Numerical second derivatives.
C
C     PRINT*,'Numerical second derivatives:'
      DIF=1.0D-4
      DO J1=1,3*N
         TEMP1=X(J1)
         X(J1)=X(J1)+DIF
         CALL TFGRAD(N, X, VEC1, ALPHAP, ALPHAM, ZSYM, DAMP)
         X(J1)=X(J1)-2.0D0*DIF
         CALL TFGRAD(N, X, VEC2, ALPHAP, ALPHAM, ZSYM, DAMP)
         DO J2=1,3*N
            HESS(J2,J1)=HESS(J2,J1)+(VEC1(J2)-VEC2(J2))/(2.0D0*DIF)
C           WRITE(*,'(A,2I3,F20.10)') 'J1,J2,A=',J1,J2,(VEC1(J2)-VEC2(J2))/(2.0D0*DIF)
         ENDDO
         X(J1)=TEMP1
      ENDDO

      RETURN
      END

C*******************************************************************************
C
C  First order induction energy
C
      SUBROUTINE TFIND(N, X, EIND, ALPHAP, ALPHAM, ZSYM, DAMP)
      IMPLICIT NONE
      INTEGER N, J1, J2, NPROD(N,N)
      DOUBLE PRECISION X(3*N), EIND, DAMP, DUMMY3,
     1                 R(N,N), 
     2                 DUMMY, ALPHAP, ALPHAM, VEC(3), 
     3                 F2, RR(N,N)
      CHARACTER(LEN=5) ZSYM(N)
C
C Statement functions
C
      F2(DUMMY)=1.0D0-DEXP(-DUMMY/DAMP)*(1.0D0+DUMMY/DAMP+DUMMY**2/(2.0D0*DAMP**2))
C     F2(DUMMY)=1.0D0
C 
C  Store distance matrices.
C
      DO J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         NPROD(J1,J1)=1
         DO J2=J1+1,N
            R(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1              +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2              +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R(J2,J1)=1.0D0/DSQRT(R(J2,J1))
            RR(J2,J1)=1.0D0/R(J2,J1)
            R(J1,J2)=R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            IF (ZSYM(J2).EQ.ZSYM(J1)) THEN
               NPROD(J2,J1)=1
               NPROD(J1,J2)=1
            ELSE
               NPROD(J2,J1)=-1
               NPROD(J1,J2)=-1
            ENDIF
         ENDDO
      ENDDO
C
C  First order induction energy.
C
      EIND=0.0D0
      DO J1=1,N
         DUMMY=0.0D0
         VEC(1)=0.0D0
         VEC(2)=0.0D0
         VEC(3)=0.0D0
         DO J2=1,N
            DUMMY3=F2(RR(J2,J1))*R(J2,J1)**3
            IF (ZSYM(J2).EQ.'PL') THEN
               VEC(1)=VEC(1)+(X(3*J1-2)-X(3*J2-2))*DUMMY3
               VEC(2)=VEC(2)+(X(3*J1-1)-X(3*J2-1))*DUMMY3
               VEC(3)=VEC(3)+(X(3*J1)-X(3*J2))*DUMMY3
            ELSE
               VEC(1)=VEC(1)-(X(3*J1-2)-X(3*J2-2))*DUMMY3
               VEC(2)=VEC(2)-(X(3*J1-1)-X(3*J2-1))*DUMMY3
               VEC(3)=VEC(3)-(X(3*J1)-X(3*J2))*DUMMY3
            ENDIF
         ENDDO
         IF (ZSYM(J1).EQ.'PL') THEN
            EIND=EIND-ALPHAP*(VEC(1)**2+VEC(2)**2+VEC(3)**2)/2.0D0
         ELSE
            EIND=EIND-ALPHAM*(VEC(1)**2+VEC(2)**2+VEC(3)**2)/2.0D0
         ENDIF
      ENDDO

      RETURN
      END

C*************************************************************************
C
C  Analytic first derivatives of possible first order induction.
C
      SUBROUTINE TFGRAD(N, X, V, ALPHAP, ALPHAM, ZSYM, DAMP)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, NPROD(N,N)
      DOUBLE PRECISION X(3*N), DAMP,
     1                 V(3*N), R(N,N), DUMMY3, DUMMY4,
     2                 DUMMY, ALPHAP, ALPHAM, DUMMY5, DUMMY6,
     3                 DOT(N,N), R2(N,N), DUMMY7, DUMMY8,
     4                 R3(N,N), XJ3, RR(N,N),
     5                 DAMP2, DAMP3, RJ4J1, DSAVE(N,N)
      CHARACTER(LEN=5) ZSYM(N)
C
C Statement functions
C
C     F2(DUMMY)=1.0D0-DEXP(-DUMMY/DAMP)*(1.0D0+DUMMY/DAMP+DUMMY**2/(2.0D0*DAMP**2))
C     F2P(DUMMY)=DEXP(-DUMMY/DAMP)*DUMMY**2/(2.0D0*DAMP**3)

      DAMP2=DAMP**2
      DAMP3=DAMP**3
C 
C  Store distance matrices.
C
      DO J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         R3(J1,J1)=0.0D0
         DOT(J1,J1)=X(3*J1-2)**2 + X(3*J1-1)**2 + X(3*J1)**2
         NPROD(J1,J1)=1
         DO J2=J1+1,N
            DUMMY=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1           +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2           +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RR(J2,J1)=DSQRT(DUMMY)
            DSAVE(J2,J1)=DEXP(-RR(J2,J1)/DAMP)
            DUMMY=1.0D0/DUMMY
            R2(J2,J1)=DUMMY
            R(J2,J1)=DSQRT(DUMMY)
            R3(J2,J1)=R(J2,J1)*DUMMY
            R(J1,J2)=R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R3(J1,J2)=R3(J2,J1)
            DSAVE(J1,J2)=DSAVE(J2,J1)
            DOT(J2,J1)=X(3*J2-2)*X(3*J1-2) + X(3*J2-1)*X(3*J1-1) + X(3*J2)*X(3*J1)
            DOT(J1,J2)=DOT(J2,J1)
            IF (ZSYM(J2).EQ.ZSYM(J1)) THEN
               NPROD(J2,J1)=1
               NPROD(J1,J2)=1
            ELSE
               NPROD(J2,J1)=-1
               NPROD(J1,J2)=-1
            ENDIF
         ENDDO
      ENDDO

      DO J1=1,3*N
         V(J1)=0.0D0
      ENDDO
C
C  Gradient.
C
      DO J1=1,N           !   J1 = l
         DO J2=1,3        !   J2 = alpha
            J3=3*J1-3+J2 
            XJ3=X(J3)
            DO J4=1,N        !  J4 = i
               DUMMY=0.0D0
               DUMMY3=DOT(J4,J4)-DOT(J4,J1)
               DUMMY4=(XJ3-X(3*J4-3+J2))*R2(J4,J1)
               DUMMY5=1.0D0-DSAVE(J4,J1)*(1.0D0+RR(J4,J1)/DAMP+RR(J4,J1)**2/(2.0D0*DAMP2))
               DUMMY6=DSAVE(J4,J1)*RR(J4,J1)**2/(2.0D0*DAMP**3)
               DUMMY7=X(3*J4-3+J2)
               DUMMY8=R(J4,J1)
               DO J5=1,N     !  J5 = k
                  DUMMY=DUMMY-(NPROD(J5,J1)*R3(J5,J4)*(1.0D0 + 
     1      DSAVE(J5,J4)*(-1.0D0 - RR(J5,J4)/DAMP - (0.5D0*RR(J5,J4)**2)/DAMP2))*
     2    (DUMMY6*DUMMY8*(DUMMY7 - XJ3)*(DUMMY3 + DOT(J5,J1) - DOT(J5,J4)) + 
     3      DUMMY5*(DUMMY7 + 3.0D0*DUMMY4*(DUMMY3 + DOT(J5,J1) - DOT(J5,J4)) - X(-3 + J2 + 3*J5))) )
               ENDDO
               IF (ZSYM(J4).EQ.'PL') THEN
                  V(J3)=V(J3)-ALPHAP*R3(J4,J1)*DUMMY
               ELSE
                  V(J3)=V(J3)-ALPHAM*R3(J4,J1)*DUMMY
               ENDIF
            ENDDO

            DUMMY=0.0D0
            DO J4=1,N            ! J4 = j
               DUMMY3=DOT(J1,J1)-DOT(J4,J1)
               DUMMY4=(XJ3-X(3*J4-3+J2))*R2(J4,J1)
               DUMMY5=R3(J4,J1)
               DUMMY6=X(3*J4-3+J2)
               DUMMY7=1.0D0-DSAVE(J4,J1)*(1.0D0+RR(J4,J1)/DAMP+RR(J4,J1)**2/(2.0D0*DAMP2))
               DUMMY8=DSAVE(J4,J1)*RR(J4,J1)**2/(2.0D0*DAMP**3)
               RJ4J1=R(J4,J1)
               DO J5=1,N         ! J5 = k
                  DUMMY=DUMMY+DUMMY5*NPROD(J5,J4)*R3(J5,J1)*
     1  ( (1.0D0-DSAVE(J5,J1)*(1.0D0+RR(J5,J1)/DAMP+RR(J5,J1)**2/(2.0D0*DAMP2)))
     2       *(DUMMY8*(-DUMMY6 + XJ3)*(DUMMY3 - DOT(J5,J1) + DOT(J5,J4))*
     2        RJ4J1 + DUMMY7*(-DUMMY6 + 2.0D0*XJ3 - 
     3          3.0D0*(DUMMY3 - DOT(J5,J1) + DOT(J5,J4))*
     4           (DUMMY4 + R2(J5,J1)*(XJ3 - X(-3 + J2 + 3*J5))) - X(-3 + J2 + 3*J5))
     5       ) + DUMMY7*(DUMMY3 - DOT(J5,J1) + DOT(J5,J4))*
     6           DSAVE(J5,J1)*RR(J5,J1)/(2.0D0*DAMP3)*(XJ3 - X(-3 + J2 + 3*J5)))
               ENDDO
            ENDDO
            IF (ZSYM(J1).EQ.'PL') THEN
               V(J3)=V(J3)-ALPHAP*DUMMY/2.0D0
            ELSE
               V(J3)=V(J3)-ALPHAM*DUMMY/2.0D0
            ENDIF
         ENDDO
      ENDDO

      RETURN 
      END
C
C  Analytic second derivatives of possible first order induction.
C
      SUBROUTINE TFSEC(N, X, ALPHAP, ALPHAM, ZSYM, DAMP, NPROD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NPROD(N,N), J7
      DOUBLE PRECISION X(3*N), DAMP,
     1                 R(N,N), 
     2                 DUMMY, ALPHAP, ALPHAM, 
     3                 DOT(N,N), R2(N,N),
     4                 R3(N,N), XJ3, XJ6, DOT11, DOT44, DOT41,
     5                 XJ1J5, XJ4J2, R2J4J1, R3J4J1, EX, AL(N), XJ5, R3J6J1,
     6                 XJ6J2, XJ6J4, DOT66, DOT61, R2J6J1, RR(N,N)
      CHARACTER(LEN=5) ZSYM(N)
C
C  Statement functions.
C
C     F2(DUMMY)=1.0D0-DEXP(-DUMMY/DAMP)*(1.0D0+DUMMY/DAMP+DUMMY**2/(2.0D0*DAMP**2))
C     F2P(DUMMY)=DEXP(-DUMMY/DAMP)*DUMMY**2/(2.0D0*DAMP**3)
C     F2PP(DUMMY)=DEXP(-DUMMY/DAMP)*DUMMY*(1.0D0-DUMMY/(2.0D0*DAMP))/DAMP**3
C 
C  Store distance matrices.
C
      DO J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DOT(J1,J1)=X(3*J1-2)**2 + X(3*J1-1)**2 + X(3*J1)**2
         IF (ZSYM(J1).EQ.'PL') THEN
            AL(J1)=ALPHAP
         ELSE
            AL(J1)=ALPHAM
         ENDIF
         DO J2=J1+1,N
            DUMMY=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1           +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2           +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RR(J2,J1)=DSQRT(DUMMY)
            DUMMY=1.0D0/DUMMY
            R2(J2,J1)=DUMMY
            R(J2,J1)=DSQRT(DUMMY)
            R3(J2,J1)=R(J2,J1)*DUMMY
            RR(J1,J2)=RR(J2,J1)
            R(J1,J2)=R(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R3(J1,J2)=R3(J2,J1)
            DOT(J2,J1)=X(3*J2-2)*X(3*J1-2) + X(3*J2-1)*X(3*J1-1) + X(3*J2)*X(3*J1)
            DOT(J1,J2)=DOT(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,3*N
         DO J2=1,3*N
            HESS(J2,J1)=0.0D0
         ENDDO
      ENDDO
C
C  Second derivatives - coordinates on the same atom.
C
      DO J1=1,N                         ! J1 = l     
         DOT11=DOT(J1,J1)
         DO J2=1,3                      ! J2 = alpha
            J3=3*J1-3+J2
            XJ3=X(J3)
            DO J4=J2,3                  ! J4 = beta
               J5=3*J1-3+J4
               XJ5=X(J5)
               EX=0.0D0
               IF (J4.EQ.J2) EX=1.0D0
               DO J6=1,N                ! J6 = i
                  DOT66=DOT(J6,J6)
                  DOT61=DOT(J6,J1)
                  R2J6J1=R2(J6,J1)
                  R3J6J1=R3(J6,J1)
                  XJ6J2=X(3*J6-3+J2)
                  XJ6J4=X(3*J6-3+J4)
                  DUMMY=2.0D0*EX*R3J6J1**2 - 12.0D0*R2J6J1**4*(XJ3 - XJ6J2)*(XJ5 - XJ6J4) + 
     1         18.0D0*(DOT11 - 2.0D0*DOT61 + DOT66)*R2J6J1**5*(XJ3 - XJ6J2)*(XJ5 - XJ6J4)
                  DO J7=1,N             ! J7 = k
                     DUMMY=DUMMY+2.0D0*R2J6J1*R3J6J1*NPROD(J7,J1)*R3(J7,J6)*
     1                    (3.0D0*EX*(DOT61 - DOT66 - DOT(J7,J1) + DOT(J7,J6)) - 
     2                  15.0D0*R2J6J1*(XJ3 - XJ6J2)*(XJ5 - XJ6J4)*(DOT61 - DOT66 - DOT(J7,J1) + DOT(J7,J6)) + 
     3                    3.0D0*(XJ5 - XJ6J4)*(XJ6J2 - X(-3 + J2 + 3*J7)) + 
     4                    3.0D0*(XJ3 - XJ6J2)*(XJ6J4 - X(-3 + J4 + 3*J7)))
                  ENDDO
                  HESS(J5,J3)=HESS(J5,J3)-AL(J6)*DUMMY/2.0D0
               ENDDO
C
C  Case l = i.
C
               DUMMY=0.0D0
               DO J6=1,N              ! J6 = j
                  R3J6J1=R3(J6,J1)
                  R2J6J1=R2(J6,J1)
                  XJ6J2=X(3*J6-3+J2)
                  XJ6J4=X(3*J6-3+J4)
                  DOT61=DOT(J6,J1)
                  DO J7=1,N           ! J7 = k
                     DUMMY=DUMMY+R3J6J1*NPROD(J7,J6)*R3(J7,J1)*(EX*
     1     (2.0D0 - 3.0D0*(DOT11 - DOT61 - DOT(J7,J1) + DOT(J7,J6))*(R2J6J1 + R2(J7,J1))) - 
     2    3.0D0*(2.0D0*XJ3 - XJ6J2 - 3.0D0*(DOT11 - DOT61 - DOT(J7,J1) + DOT(J7,J6))*
     3        (R2J6J1*(XJ3 - XJ6J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7))) - X(-3 + J2 + 3*J7))*
     4     (R2J6J1*(XJ5 - XJ6J4) + R2(J7,J1)*(XJ5 - X(-3 + J4 + 3*J7))) + 
     5    6.0D0*(DOT11 - DOT61 - DOT(J7,J1) + DOT(J7,J6))*
     6     (R2J6J1**2*(XJ3 - XJ6J2)*(XJ5 - XJ6J4) + 
     7       R2(J7,J1)**2*(XJ3 - X(-3 + J2 + 3*J7))*(XJ5 - X(-3 + J4 + 3*J7))) - 
     8    3.0D0*(R2J6J1*(XJ3 - XJ6J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7)))*
     9     (2.0D0*XJ5 - XJ6J4 - X(-3 + J4 + 3*J7)))
                  ENDDO 
               ENDDO
               HESS(J5,J3)=HESS(J5,J3)-AL(J1)*DUMMY/2.0D0
               HESS(J3,J5)=HESS(J5,J3)
            ENDDO
         ENDDO
      ENDDO
C
C  Second derivatives - different atoms.
C
      DO J1=1,N             ! J1 = l     
         DOT11=DOT(J1,J1)
         DO J2=1,3          ! J2 = alpha
            J3=3*J1-3+J2 
            XJ3=X(J3)
            DO J4=J1+1,N    ! J4 = m
               DOT44=DOT(J4,J4)
               DOT41=DOT(J4,J1)
               XJ4J2=X(3*J4-3+J2)
               R2J4J1=R2(J4,J1)
               R3J4J1=R3(J4,J1)
               DO J5=1,3    ! J5 = beta
                  J6=3*J4-3+J5
                  XJ1J5=X(3*J1-3+J5)
                  XJ6=X(J6)
                  EX=0.0D0
                  IF (J5.EQ.J2) EX=1.0D0
                  DO J7=1,N   ! J7 = i
                     DUMMY=2.0D0*NPROD(J4,J1)*R3(J7,J1)*R3(J7,J4)*
     1  (EX - 3.0D0*(R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7))*(XJ1J5 - X(-3 + J5 + 3*J7)) + 
     2       R2(J7,J4)*(XJ4J2 - X(-3 + J2 + 3*J7))*(XJ6 - X(-3 + J5 + 3*J7))) + 
     3    9.0D0*(DOT41 - DOT(J7,J1) - DOT(J7,J4) + DOT(J7,J7))*R2(J7,J1)*R2(J7,J4)*
     -     (XJ3 - X(-3 + J2 + 3*J7))*(XJ6 - X(-3 + J5 + 3*J7)))
                     HESS(J6,J3)=HESS(J6,J3)-AL(J7)*DUMMY/2.0D0
                  ENDDO
C
C  i = l term
C
                  DUMMY=0.0D0
                  DO J7=1,N    ! J7 = j
                     DUMMY=DUMMY+2.0D0*R3J4J1*NPROD(J7,J4)*R3(J7,J1)*
     1  (-EX + 3.0D0*EX*R2J4J1*(DOT11 - DOT41 - DOT(J7,J1) + DOT(J7,J4)) - 
     2    3.0D0*R2J4J1*(-XJ1J5 + XJ6)*(2.0D0*XJ3 - XJ4J2 - 
     3       3.0D0*(DOT11 - DOT41 - DOT(J7,J1) + DOT(J7,J4))*
     4        (R2J4J1*(XJ3 - XJ4J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7))) - X(-3 + J2 + 3*J7))
     5     - 3.0D0*(2.0D0*R2J4J1**2*(XJ3 - XJ4J2)*(XJ1J5 - XJ6)*
     6        (DOT11 - DOT41 - DOT(J7,J1) + DOT(J7,J4)) + 
     7       (R2J4J1*(XJ3 - XJ4J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7)))*
     8        (-XJ1J5 + X(-3 + J5 + 3*J7))))
                  ENDDO
                  HESS(J6,J3)=HESS(J6,J3)-AL(J1)*DUMMY/2.0D0
C
C  i = m term
C
                  DUMMY=0.0D0
                  DO J7=1,N    ! J7 = j
                     DUMMY=DUMMY+2.0D0*R3J4J1*NPROD(J7,J1)*R3(J7,J4)*
     1  (-EX - 3.0D0*EX*R2J4J1*(DOT41 - DOT44 - DOT(J7,J1) + DOT(J7,J4)) - 
     2    3.0D0*(-2.0D0*R2J4J1**2*(XJ3 - XJ4J2)*(XJ1J5 - XJ6)*
     3        (DOT41 - DOT44 - DOT(J7,J1) + DOT(J7,J4)) + 
     4       (-XJ4J2 + X(-3 + J2 + 3*J7))*
     5        (R2J4J1*(-XJ1J5 + XJ6) + R2(J7,J4)*(XJ6 - X(-3 + J5 + 3*J7)))) - 
     6    3.0D0*R2J4J1*(XJ3 - XJ4J2)*(-XJ1J5 + 2.0D0*XJ6 - 
     7       3.0D0*(-DOT41 + DOT44 + DOT(J7,J1) - DOT(J7,J4))*
     8        (R2J4J1*(-XJ1J5 + XJ6) + R2(J7,J4)*(XJ6 - X(-3 + J5 + 3*J7))) - X(-3 + J5 + 3*J7)))
                  ENDDO
                  HESS(J6,J3)=HESS(J6,J3)-AL(J4)*DUMMY/2.0D0
                  HESS(J3,J6)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN 
      END

